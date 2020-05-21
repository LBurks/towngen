import random
import time
from math import sqrt

import numpy as np
import scipy.spatial as spatial

from fortune import fortune
from draw import drawFortune

def sampledata(max_points = 100):
    xfile = open('x.sample', 'r')
    yfile = open('y.sample', 'r')
    limit = max_points
    return [[float(x), float(y)] for x, y, i in zip(xfile, yfile, range(limit))]

def randomdata(max_points):
    random.seed(int(time.time()))
    xset = [random.uniform(.2, .8) for x in range(max_points)]
    yset = [random.uniform(.2, .8) for y in range(max_points)]
    return np.array([[x, y] for x, y in zip(xset,yset)])

def voronoi_finite_polygons_2d(vor, radius=None):
    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

#def fortunetest():
#   points = sampledata(15)
#    print(points)
#    diagram = spatial.Voronoi(points) 
#    regions, vertices = voronoi_finite_polygons_2d(diagram)
#    vertices, regions, center = fortune(points)
#   drawFortune(points, vertices, regions, center)

def fortunetest1():
    points = sampledata(10)
    #print(points)
    diagram = spatial.Voronoi(points) 
    regions, vertices = voronoi_finite_polygons_2d(diagram)
    fortune(points, [0, 0, 1, 1])

fortunetest1()
