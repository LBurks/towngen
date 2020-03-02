import numpy as np
import math
import random
from draw import drawFortune
from collections import deque

class arc:
    def __init__(self, point, a = None, b = None):
        self.focus = point
        self.prev = a
        self.next = b
        self.seg1 = self.seg2 = None
        self.event = None

class segment:
    def __init__(self, point):
        self.start = point
        self.end = None
        self.done = False
        output.append(self)

    def finish(self, point):
        if self.done: return
        self.end = point
        self.done = True

class event:
    def __init__(self, x, p, a):
        self.x = x
        self.point = p
        self.arc = a
        self.valid = True

points = [[0.25, 0.5], [0.3, 0.2]]
X0 = 0
X1 = 1
Y0 = 0
Y1 = 1
points.sort()
qpoints = deque(points)
events = deque([])
output = []
root = None

def siteEvent(point):
    global root
    if root is None:
        root = arc(point)
        return
    i = root
    while i is not None:
        meet, z = intersect(point, i)
        if (meet):
            if (i.next is not None and not intersect(point, i.next)[0]):
                i.next.prev = arc(i.focus, i, i.next)
                i.next = i.next.prev
            else: i.next = arc(i.focus, i)
            i.next.seg1 = i.seg1
            i.next.prev = arc(point,i,i.next)
            i.next = i.next.prev
            i = i.next
            i.prev.seg2 = i.seg1 = segment(z)
            i.next.seg1 = i.seg2 = segment(z)

            checkCircleEvent(i, point[0])
            checkCircleEvent(i.prev, point[0])
            checkCircleEvent(i.next, point[0])

            return
        i = i.next
    i = root
    while i.next is not None:
        i = i.next
    i.next = arc(point, i)
    start = [0, (i.next.point.y + i.focus.y) / 2]
    i.seg2 = i.next.seg1 = segment(start)

def intersect(point, i):
    res = [0, 0]
    if (i.focus[0] == point[0]): return False
    a = None
    b = None
    if i.prev is not None:
        a = intersection(i.prev.focus, i.focus, point[0])[1]
    if i.next is not None:
        b = intersection(i.focus, i.next.focus, point[0])[1]
    if ((i.prev is None or a <= point[1]) and (i.next is None or point[1] <= b)):
        res[1] = point[1]
        res[0] = ((pow(i.focus[0], 2) + pow(i.focus[1]-res[1], 2)) - pow(point[1], 2)) / (2 * i.focus[0] - 2 * point[0])
        return True, res
    return False, res

def checkCircleEvent(i, x0):
    if i.event is not None and i.event.point[0] != x0:
        i.event.valid = False
    i.event = None
    if (i.prev is None) or (i.next is None):
        return;
    res, x, o = circle(i.prev.focus, i.focus, i.next.focus)
    if res and x > x0:
        i.event = event(x, o, i)
        events.append(i.event)

def intersection(p0, p1, l):
    res = [0, 0]
    p = p0
    if p0[0] == p1[0]:
        res[1] = (p0[1] + p1[1]) / 2
    elif p1[0] == l:
        res[1] = p1[1]
    elif p0[0] == l:
        res[1] = p0[1]
        p = p1
    else:
        z0 = 2*(p0[0]- l)
        z1 = 2*(p1[0]- l)
        a = 1/z0 - 1/z1
        b = -2 * (p0[1]/z0 - p1[1]/z1)
        c = (pow(p0[1], 2) + pow(p0[0], 2) - pow(l, 2))/z0 - (pow(p1[1], 2) + pow(p1[0], 2) - pow(l, 2))/z1
    res[0] = (pow(p[0], 2) + pow(p[1] - res[1], 2) - pow(l, 2))/(2*p[0]-2*l)
    return res

def circle(a, b, c):
    x = [0, 0]
    o = [0, 0]
    if ((b[0]-a[0]) * (a[1]-a[1]) - (c[0]-a[0]) * (b[1]-a[1])) > 0:
        return False, x, o
    A = b[0] - a[0]
    B = b[1] - a[1]
    C = c[0] - a[0]
    D = c[1] - a[1]
    E = A*(a[0]+b[0]) + B*(a[1]+b[1])
    F = C*(a[0]+c[0]) + D*(a[1]+c[1])
    G = 2*(A*(c[1]-b[1])) - 2*(B*(c[0]-b[0]))

    if G == 0: return False, x, o
    o[0] = (D*E-B*F)/G
    o[1] = (A*F-C*E)/G

    x = o[0] + math.sqrt(pow(a[0] - o[0], 2) + pow(a[1] - o[1], 2))
    return True, x, o

def vertexEvent(event):
    if event.valid:
        seg = segment(event.point)
        a = event.arc
        if a.prev is not None:
            a.prev.next = a.next
            a.prev.seg2 = seg
        if a.next is not None:
            a.next.prev = a.prev
            a.prev.seg1 = seg
        if a.seg1 is not None: a.seg1.finish(event.point)
        if a.seg2 is not None: a.seg2.finish(event.point)
        if a.prev is not None: checkCircleEvent(a.prev, event.x)
        if a.next is not None: checkCircleEvent(a.next, event.x)

def finish_edges():
    l = X0 + (X1-X0) + (Y1-Y0)
    i = root
    while i is not None:
        if i.seg2 is not None:
            i.seg2.finish(intersection(i.focus, i.next.focus, l*2))
        i = i.next

def print_output():
    print(X0, X1, Y0, Y1)
    for s in output:
        p0 = s.start
        p1 = s.end
        print(p0, p1)


class diagram():
    def __init__(self, points):
        self.points = np.copy(points)
        self.vertices = np.array([])
        self.regions = np.array([])
        #Neighboring regions and the vertices that form their boundary
        self.neighbor_regions = {}
        #Neighboring vertices and the points that form their boundary
        self.vertex_neighbors = {}


def fortune():
    while len(qpoints) > 0:
        if (len(events) > 0 and events[0].point[0] < qpoints[0][0]):
            event = events.popleft()
            vertexEvent(event)
        else:
            point = qpoints.popleft()
            siteEvent(point)
    while len(events) > 0:
        event = events.popleft()
        vertexEvent(event)
    print("here")
    finish_edges()
    print_output()
    drawFortune(points, output)

fortune()
