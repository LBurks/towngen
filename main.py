#! /usr/bin/python

import numpy
from numpy import median
import scipy.spatial as spatial
import cairo
import math
import random
from shapely.geometry import Polygon, LineString


def calcIntersection(p1, p2):
    for x in [0, 1]:
        t = (x - p1[0]) / (p2[0] - p1[0])
        print(t)
        if 0 <= t <= 1:
            print("x is equal to", x)
            return [x, p1[1] + t * (p2[1] - p2[1])]
    for y in [0, 1]:
        t = (y - p1[1]) / (p2[1] - p1[0])
        if 0 <= t <= 1:
            print("y is equal to", y)
            return [p1[0] + t * (p2[0] - p2[0]), y]
    print(p1, p2, "Does not intersect the bounding box")
    return p2


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
        print("---")
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
        print(ridges)
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

def gen_map(starting_points, iterations):
    box = Polygon([[0, 0], [0, 1], [1, 1], [1, 0]])
    points = starting_points
    centers = []
    for i in range(iterations):
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
    result = {'points': points, 'centers': centers, 'regions': regions, 'vertices': vertices, 'ridges': ridges }#, 'ridge_points': ridge_points, 'ridge_vertices': ridge_vertices}
    return diagram, result


# This generates random points for the map.
def randomdata(max_points):
    xset = [random.uniform(0, 1) for x in range(max_points)]
    yset = [random.uniform(0, 1) for y in range(max_points)]
    return numpy.array([[x, y] for x, y in zip(xset,yset)])

# This gets sample data from some text files.
# The purpose is to get a consistent map to test results of various features i'm working on.
# This will go away in the future
def sampledata(max_points):
    xfile = open('x.sample', 'r')
    yfile = open('y.sample', 'r')
    limit = min(max_points, 100)
    return numpy.array([[float(x), float(y)] for x, y, i in zip(xfile, yfile, range(limit))])

MAX_POINTS = 100
WIDTH = 100.0
HEIGHT = 100.0
yratio = 1 #float(HEIGHT) / 200
xratio = 1 #float(WIDTH) / 200

old_diagram, diagram  = gen_map(sampledata(MAX_POINTS), 3)

surface = cairo.SVGSurface("example.svg", WIDTH, HEIGHT)
context = cairo.Context(surface)
context.scale(WIDTH, HEIGHT)
context.set_line_width(0.001)
context.set_font_size(0.005)

box = Polygon([[0, 0], [0, 1], [1, 1], [1, 0]])

for r in diagram['regions']:
    context.set_source_rgba(0, 0, 0, 1)
    context.move_to(diagram['vertices'][r[0]][0], diagram['vertices'][r[0]][1])
    for v in diagram['vertices'][r[1:]]:
        context.line_to(v[0], v[1])
    context.close_path()
    context.stroke()

for r in diagram['ridges']:
    ridges = diagram['ridges'][r]
    for p2, v1, v2 in ridges:
        context.move_to(diagram['points'][r][0], diagram['points'][r][1])
        ridge = LineString([diagram['vertices'][v1], diagram['vertices'][v2]])
        context.line_to(diagram['points'][p2][0], diagram['points'][p2][1])
        context.stroke()

# for p1, p2 in diagram['ridge_points']:
#     context.set_line_width(0.005)
#     context.move_to(diagram['points'][p1][0], diagram['points'][p1][1])
#     context.line_to(diagram['points'][p2][0], diagram['points'][p2][1])
#     context.stroke()
#
# for v1, v2 in diagram['ridge_vertices']:
#     context.set_line_width(0.001)
#     context.move_to(diagram['vertices'][v1][0], diagram['vertices'][v1][1])
#     context.line_to(diagram['vertices'][v2][0], diagram['vertices'][v2][1])
#     context.stroke()


context.set_line_width(0.001)

context.set_source_rgb(0.7, 0.2, 0.2)

for p in diagram['points']:
    context.move_to(0, 0)
    context.arc(p[0], p[1], 0.004, 0, 2*numpy.pi)
    context.fill()
    context.stroke()

context.set_source_rgb(0.2, 0.2, 0.7)

for p in diagram['vertices']:
    context.move_to(0, 0)
    context.arc(p[0], p[1], 0.004, 0, 2*numpy.pi)
    context.fill()
    context.stroke()

# for v in vertices:#point in points:
#     context.move_to(v[0], v[1])
#     context.arc(v[0], v[1], 0.02, 0, 2*numpy.pi)
#     context.fill()
#     context.stroke()

# for r in regions:
#     print("poly:", r)
#     for i in range(len(r)-1):
#         print(intersect(vertices[r[i]], vertices[r[i+1]]))
#     print(intersect(vertices[r[-1]], vertices[r[0]]))



# context.set_source_rgb(0.2, 0.2, 0.7)
#
#
# for i in range(len(regions)):#r in regions:
#     polyverts = [vertices[v] for v in regions[i]]
#     centroid = find_centroid(polyverts)
#     context.move_to(centroid[0], centroid[1])
#     context.arc(centroid[0], centroid[1], 0.004, 0, 2*numpy.pi)
#     context.fill()
# context.stroke()


# print("***DIAGRAM***")
# for o in diagram:
#     print(o, diagram[o])
