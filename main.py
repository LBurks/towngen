#! /usr/bin/python

import numpy
from numpy import median
import numpy as np
import scipy.spatial as spatial
import cairo
import math
import random
import time
from shapely.geometry import Polygon, LineString
import noise
import sys
from draw import drawOutlines, drawPolys

def voronoi_def_polys_2d(diagram, radius=None):
    box = Polygon([[0, 0], [0, 1], [1, 1], [1, 0]])
    if diagram.points.shape[1] != 2:
        raise ValueError("Requires 2D input")
    new_regions = []
    new_vertices = diagram.vertices.tolist()
    center = diagram.points.mean(axis=0)
    if radius is None:
        radius = diagram.points.ptp().max()
    all_ridges = {}
    for (p1, p2), (v1,v2) in zip(diagram.ridge_points, diagram.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))
    for p1, region in enumerate(diagram.point_region):
        vertices = diagram.regions[region]
        if all(v >= 0 for v in vertices):
            polygon = Polygon([diagram.vertices[v] for v in vertices])
            polygon = polygon.intersection(box)
            resultcoord = [[p[0], p[1]] for p in polygon.exterior.coords]
            final_region = []
            for v in resultcoord:
                if v in new_vertices:
                    final_region.append(new_vertices.index(v))
                else:
                    final_region.append(len(new_vertices))
                    new_vertices.append(v)
            new_regions.append(final_region)
            continue
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]
        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                continue
            t = diagram.points[p2] - diagram.points[p1]
            t /= numpy.linalg.norm(t)
            n = numpy.array([-t[1], t[0]])
            midpoint = diagram.points[[p1,p2]].mean(axis=0)
            direction = numpy.sign(numpy.dot(midpoint - center, n)) * n
            far_point = diagram.vertices[v2] + direction * radius
            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())
            all_ridges[p1][ridges.index((p2, v1, v2))] = (p2, len(new_vertices) - 1, v2)
        vs = numpy.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = numpy.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = numpy.array(new_region)[numpy.argsort(angles)]
        polygon = Polygon([new_vertices[v] for v in new_region])
        polygon = polygon.intersection(box)
        resultcoord = [[p[0], p[1]] for p in polygon.exterior.coords]
        final_region = []
        for v in resultcoord:
            if v in new_vertices:
                final_region.append(new_vertices.index(v))
            else:
                final_region.append(len(new_vertices))
                new_vertices.append(v)
        new_regions.append(final_region)
    for p1 in all_ridges:
        targets = []
        for p2, v1, v2 in all_ridges[p1]:
            ridge = LineString([new_vertices[v1], new_vertices[v2]])
            if len(box.intersection(ridge).bounds) <= 0:
                targets.append(all_ridges[p1].index((p2, v1, v2)) - len(targets))
        for i in targets:
            all_ridges[p1].pop(i)
    return new_regions, numpy.asarray(new_vertices), all_ridges

def find_centroid(vertices):
    area = 0
    center_x = 0
    center_y = 0
    polygon = vertices.tolist()
    polygon.append(polygon[0]) #Because p(n) and p(0) need to be part of the calculations as well
    for i in range(len(polygon) - 1):
        vert0 = polygon[i]
        vert1 = polygon[i+1]
        step =  (vert0[0]*vert1[1]) - (vert1[0]*vert0[1])
        area += step
        center_x += (vert0[0] + vert1[0]) * step
        center_y += (vert0[1] + vert1[1]) * step
    area /= 2
    return [center_x / (6.0*area), center_y / (6.0*area)]

class point():
    def __init__(self, co):
        self.coordinates = co
        self.land = False
        self.water = False
        self.ocean = False
        self.elevation = 0
        self.border = False
        self.visible = True

class voronoiDiagram():
    def __init__(self):
        # The points used to generate the voronoi diagram. This should be the point class above
        self.points = []
        # The centers of the polygons that make up the voronoi diagram. Should be the point class above
        self.centers = []
        # arrays of the vertices that makes up each polygon of the diagram. Should be the point class above
        self.regions = []
        # the coordinates for each vertex
        self.vertices = []
        # A set to define each ridge, the key is the point with an array of tuples (p2, v1, v2)
        self.ridges = {}

# A flood fill algorithm that starts by looking for the polygon in the upper-left corner.
# It then assumes that is an ocean tile and looks for neighbors
#  that are both water and don't already have the ocean flag raised
def floodFill(diagram):
    lens = [math.sqrt((pow(p.coordinates[0], 2) + pow(p.coordinates[1], 2))) for p in diagram.points]
    corner = lens.index(min(lens))
    diagram.points[corner].ocean = True
    oceans = [corner]
    changes = [corner]
    for p1 in oceans:
        for p2, v1, v2 in diagram.ridges[p1]:
            if diagram.points[p2].water == True and diagram.points[p2].ocean == False:
                diagram.points[p2].ocean = True
                oceans.append(p2)

# An attempt to emulate the
def getElevations(diagram):
    queue = []
    for v in diagram.vertices:
        if v.border:
            v.elevation = 0
            queue.append(diagram.vertices.index(v))
    for q in queue:
        for p2, v1, v2 in diagram.ridges[q]:
            newElevation = 0.01 + diagram.points[q].elevation
            if (not diagram.points[q].water and not diagram.points[p2].water):
                newElevation += 1
            if newElevation > diagram.points[p2].elevation:
                diagram.points[p2].elevation = newElevation
                queue.append(p2)
    return 0;

def getVertexNeighbors(diagram):
    point = random.randint(0, len(diagram.vertices))
    segments = []
    for p1 in diagram.ridges:
        for p2, v1, v2 in diagram.ridges[p1]:
            if (v1 == point or v2 == point) and ([v1, v2] not in segments or [v2, v1] not in segments):
                segments.append([v1, v2])
    print(point)
    print(segments)

def gen_map(starting_points, iterations):
    box = Polygon([[0, 0], [0, 1], [1, 1], [1, 0]])
    points = starting_points
    centers = []
    for i in range(iterations):
        print("Generation", i)
        if i > 0:
            points = centers
        diagram = spatial.Voronoi(points)
        regions, vertices, ridges = voronoi_def_polys_2d(diagram)
        centers=[]
        for r in regions:
            poly = Polygon(vertices[r])
            poly = poly.intersection(box)
            polygon = [p for p in poly.exterior.coords]
            center = find_centroid(vertices[r])
            centers.append(center)
    testDiagram = voronoiDiagram()
    testDiagram.regions = regions
    testDiagram.ridges = ridges
    for p in vertices:
        target = point(p)
        target.land = perlinForm(target.coordinates)
        target.water = not target.land
        target.ocean = False
        target.border = (p[0] <= 0 or p[0] >= 1) ^ (p[1] <= 0 or p[1] >= 1)
        testDiagram.vertices.append(target)
    for p in points:
        target = point(p)
        region = testDiagram.regions[len(testDiagram.points)]
        vert_land = [testDiagram.vertices[v].land for  v in region]
        target.land = (vert_land.count(True) >= vert_land.count(False))
        target.water = not target.land
        target.ocean = False
        testDiagram.points.append(target)
    floodFill(testDiagram)
    getVertexNeighbors(testDiagram)
    # getElevations(testDiagram)
    print("Done generating diagram")
    return testDiagram


# This generates random points for the map.
def randomdata(max_points):
    random.seed(int(time.time()))
    xset = [random.uniform(0, 1) for x in range(max_points)]
    yset = [random.uniform(0, 1) for y in range(max_points)]
    return np.array([[x, y] for x, y in zip(xset,yset)])

# This gets sample data from some text files.
# The purpose is to get a consistent map to test results of various features i'm working on.
# This will go away in the future
def sampledata(max_points):
    xfile = open('x.sample', 'r')
    yfile = open('y.sample', 'r')
    limit = min(max_points, 100)
    return np.array([[float(x), float(y)] for x, y, i in zip(xfile, yfile, range(limit))])



#This is a blind translation from https://github.com/amitp/mapgen2/blob/4394df0e04101dbbdc36ee1e61ad7d62446bb3f1/Map.as
#The comments are left as is as well
def radialForm(point):
    #This is to treat the center as (0, 0) and change the range to properly create a centered island
    x = (point[0] - 0.5) * 3
    y = (point[1] - 0.5) * 3
    ISLAND_FACTOR = 1.06 #1 leads to no small islands, 2 leads to a lot
    bumps = random.randint(1, 6)
    startAngle = random.uniform(0, 2*numpy.pi)
    dipAngle = random.uniform(0, 2*numpy.pi)
    dipWidth = random.uniform(0.2, 0.7)
    angle = math.atan2(y, x)
    length = 0.5 * (max([math.fabs(point[0]), math.fabs(point[1])]) + math.sqrt(math.pow(x, 2) + math.pow(y, 2)))
    r1 = 0.5 + 0.4*math.sin(startAngle + bumps + math.cos((bumps+3) * angle))
    r2 = 0.7 - 0.2*math.sin(startAngle + bumps - math.sin((bumps+2) * angle))
    if math.fabs(angle - dipAngle) < dipWidth or math.fabs(angle - dipAngle + 2*numpy.pi) < dipWidth or math.fabs(angle - dipAngle - 2*numpy.pi) < dipWidth:
        r1 = r2 = 0.2
    return (length < r1 or (length > r1 * ISLAND_FACTOR and length < r2))

def perlinForm(point):
    random.seed(int(time.time()))
    x = (point[0] - 0.5) * 3
    y = (point[1] - 0.5) * 3
    seed = random.uniform(0, 2048)
    # seed = 771.048769494
    c = noise.pnoise3(x, y, seed, 8)
    length = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
    return (math.fabs(c) < (0.3 - 0.3 * math.pow(length, 2)))

if len(sys.argv) > 1:
    MAX_POINTS = int(sys.argv[1])
else:
    MAX_POINTS = 200
WIDTH = 1000.0
HEIGHT = 1000.0
yratio = 1 #float(HEIGHT) / 200
xratio = 1 #float(WIDTH) / 200

random.seed(int(time.time()))

points = sampledata(MAX_POINTS)
sorted_points = points[points[:, 0].argsort()]
print(sorted_points)
diagram  = gen_map(points, 1)

drawPolys(diagram, "map")
drawOutlines(diagram, "outlines")
