import numpy as np
import math
import random
from draw import drawFortune
from collections import deque

X0 = 0
X1 = 1
Y0 = 0
Y1 = 1
events = deque([])
beachline = []
vertices = []
regions = []
mypoints = np.array([])

class arc:
    def __init__(self, point, l = None, r = None):
        self.focus = point
        self.lpoint = l
        self.rpoint = r
        self.event = None

#center is the point that will be added to vertices
#arcs are the arcs that form the event, as a list.
#x is the x coordinate where the event will take place
#valid is weather is is a valid event or not. I might not even need it since I can just drop the item from the queue
class event:
    def __init__(self, c, a, x):
        self.center = c
        self.arc = a
        self.x = x

def newEvent(e):
    for q in range(len(events)):
        if e.x < events[q].x: 
            events.insert(q, e)
            return
    events.append(e)


def printBeachline(points):
    print("printing current beachline")
    for b in beachline:
        print("address:", b)
        print("Focus:", points.index(b.focus))
        print("Event:", b.event)
        print("Left point:", b.lpoint)
        print("Right point:", b.rpoint)
        print("---")

def siteEvent(point):
    if len(beachline) == 0:
        beachline.append(arc(point))
        return
    for i in range(len(beachline)):
        ise, z = intersect(point, i)
        if ise:
            if (i < len(beachline) - 1 & !intersect(point, i+1)[0]):
                beachline[i+1:i+1] = [arc(point), arc(beachline[i].focus, r=beachline[i].rpoint)]
                beachline[i].rpoint = None
                if beachline[i-1].event is not None:
                    events.remove(beachline[i-1].event)
                    beachline[i-1].event = None
                if beachline[i].event is not None:
                    events.remove(beachline[i].event)
                    beachline[i].event = None
                if beachline[i+1].event is not None:
                    events.remove(beachline[i+1].event)
                    beachline[i+1].event = None
            else : 
            r, c, d = circle(i-1)
            if r:
                e = event(c, beachline[i-1], c[0] + d)
                beachline[i-1].event = e
                newEvent(e)
            r, c, d = circle(i)
            if r:
                e = event(c, beachline[i], c[0] + d)
                beachline[i].event = e
                newEvent(e)
            r, c, d = circle(i+1)
            if r:
                e = event(c, beachline[i+1], c[0] + d)
                beachline[i+1].event = e
                newEvent(e)
            return

def vertEvent(ev):
    global events 
    if (ev.center[0] > 1 or ev.center[0] < 0 or ev.center[1] > 1 or ev.center[1] < 0):
        print("truncate ")
    vertices.append(ev.center)
    i = beachline.index(ev.arc)
    larc = beachline[i-1]
    carc = ev.arc
    rarc = beachline[i+1]
    if rarc.event is not None:
        events.remove(rarc.event)
        rarc.event = None
    if larc.event is not None:
        events.remove(larc.event)
        larc.event = None
    del beachline[i]
    events = deque(filter(lambda e: e.arc != ev.arc, events))
    larc.rpoint = len(vertices) - 1
    rarc.lpoint = len(vertices) - 1
    r, c, d = circle(i)
    if r:
        e = event(c, beachline[i], c[0] + d)
        beachline[i].event = e
        newEvent(e)
    r, c, d = circle(i-1)
    if r:
        e = event(c, beachline[i-1], c[0] + d)
        beachline[i-1].event = e
        newEvent(e)
    return i

#where does a new parabola with focus p intersect with arc i?
def intersect(p, i):
    larc = beachline[i-1] if i > 0 else None
    carc = beachline[i]
    rarc = beachline[i+1] if i < len(beachline)-1 else None 
    res = [-1, -1]
    if (beachline[i].focus[0] == p[0]): return False, res
    a = None
    b = None
    if (larc is not None):
        a = intersection(larc.focus, carc.focus, p[0])[1]
    if (rarc is not None):
        b  = intersection(carc.focus, rarc.focus, p[0])[1]
    if (larc is None or a <= p[1]) and (rarc is None or p[1] <= b):
        res[1] = p[1]
        res[0] = (pow(carc.focus[0],2) + pow(carc.focus[1] - res[1], 2) - pow(p[0], 2))/(2*carc.focus[0] - 2*p[0])
        return True, res
    return False, res

#asking where two parabolas with foci p0/p1 and directrix x=l intersect?
def intersection(p0, p1, l):
    res = p = p0.copy()
    #p0/p1 sit on the same x coordinate, so y will be the middle of them
    if (p0[0] == p1[0]):
        res[1] = (p0[1] + p1[1]) / 2
    #if p1 rests on the line    
    elif (p1[0] == l):
        res[1] = p1[1]
    #if p0 rests on the line
    elif (p0[0] == l):
        res[1] = p0[1]
        p = p1.copy()
    else:
        z0 = 2*(p0[0] - l)
        z1 = 2*(p1[0] - l)
        a = 1/z0 - 1/z1
        b = -2*(p0[1]/z0 - p1[1]/z1)
        c = (pow(p0[1], 2) + pow(p0[0], 2) - pow(l, 2))/z0 - (pow(p1[1], 2) + pow(p1[0], 2) - pow(l, 2))/z1
        res[1] = ( -b - math.sqrt(pow(b, 2) - 4 * a * c)) / (2 * a)
    res[0] = (pow(p[0], 2) + pow(p[1] - res[1], 2) - pow(l, 2))/(2*p[0]-2*l)
    return res

#i is the index of the arc on the beachline we're creating a circle event for
def circle(i):
    #a, b, c are the arcs' foci. left, center, and right respectively
    a = beachline[i-1].focus
    b = beachline[i].focus
    c = beachline[i+1].focus
    print(mypoints[mypoints.index(a)])
    print(mypoints[mypoints.index(b)])
    print(mypoints[mypoints.index(c)])
    if (b[0]-a[0]) * (c[1]-a[1]) - (c[0]-a[0]) * (b[1]-a[1]) > 0:
        print("This won't form a useful circle.")
        return False, [], -1

    A = b[0] - a[0]
    B = b[1] - a[1]
    C = c[0] - a[0]
    D = c[1] - a[1]
    E = A*(a[0] + b[0]) + B*(a[1] + b[1])
    F = C*(a[0] + c[0]) + D*(a[1] + c[1])
    G = 2*(A*(c[1]-b[1]) - B*(c[0]-b[0]))
    if G == 0: 
        print("points are colinear")
        return False, [], -1

    xc = (D*E-B*F)/G
    yc = (A*F-C*E)/G

    distance = math.sqrt(pow((xc - b[0]), 2) + pow((yc - b[1]), 2))
    print("returns an event")
    return True, [xc, yc], distance



#    mid1 = np.array([a[0] + b[0], a[1] + b[1]])/2
#    mid2 = np.array([b[0] + c[0], b[1] + c[1]])/2
#
#    m1 = (b[0] - a[0])/(a[1] - b[1])
#    m2 = (b[0] - c[0])/(c[1] - b[1])
#
#    p1 = mid1[1] - (m1 * mid1[0])
#    p2 = mid2[1] - (m2 * mid2[0])
#
#    xc = (p1-p2)/(m2-m1)
#    yc = p1 + (m1 * xc)
#    distance = math.sqrt(pow((xc - b[0]), 2) + pow((yc - b[1]), 2))
#    return True, [xc, yc], distance

#i = arc.
#i = an x coordiante
def checkCircleEvent(i, x0):
    if (i.event is not None and i.event.x != x0):
        print("hi")

def wrapUp():
    print(mypoints)
    center = [sum([p[0] for p in mypoints]) / len(mypoints), sum([p[1] for p in mypoints]) / len(mypoints)]
    for i in range(len(beachline) - 1):
        n1 = mypoints.index(beachline[i].focus)
        n2 = mypoints.index(beachline[i+1].focus)
        if beachline[i].rpoint is None and beachline[i+1].lpoint is None:
            print("find a point they both share, then stretch line to bounding box")
            print("n1:", regions[n1])
            print("n2:", regions[n2])
            print("target:", list(set(regions[n1]).intersection(regions[n2])))
            target = list(set(regions[n1]).intersection(regions[n2]))[0]
            beachline[i].rpoint = target
            beachline[i+1].lpoint = target
        if beachline[i].rpoint == beachline[i+1].lpoint:
            point = vertices[beachline[i].rpoint]
            if (0 < point[0] < 1 and 0 < point[1] < 1):
                print("Stretch line to bounding box")
                print("1. get the perpendicular line between the two points")
                tangent = np.array([mypoints[n2][0] - mypoints[n1][0], mypoints[n2][1] - mypoints[n1][1]])
                tangent /= np.linalg.norm(tangent)
                print(tangent)
                normal = np.array([-tangent[1], tangent[0]])
                midpoint = np.array([mypoints[n1], mypoints[n2]]).mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, normal)) * normal
                print(direction)
                print("2. draw perpendicular line from target vertex to nearest edge of bounding box")
            else:
                print("ignore this point, it won't lead to anythiing inside the bounding box'")

def wrapUp1():
    print(mypoints)
    #1. Find the center of the points. 
    center = [sum([p[0] for p in mypoints]) / len(mypoints), sum([p[1] for p in mypoints]) / len(mypoints)]
    #2. Find the slope of the line created by the two nodes. 
    #Get the inverse to find the slope of the ridge we're creating
    #3. Find which line this should intersect. 
    #3a. Compare the line between the center and the vertex (or the two nodes)
    #3b. 
    return center



def fortune(points):
    global mypoints
    global regions
    mypoints = points
    qpoints = deque(sorted(points, key=lambda p: p[0]))
    regions = [[] for x in range(len(points)) ]
    finalregions = []
    while len(qpoints) > 0:
        if (len(events) > 0 and events[0].x < qpoints[0][0]):
            e = events.popleft()
            i = vertEvent(e)
            regions[points.index(e.arc.focus)].append(len(vertices) - 1)
            regions[points.index(beachline[i-1].focus)].append(len(vertices) - 1)
            regions[points.index(beachline[i].focus)].append(len(vertices) - 1)
        else:
            point = qpoints.popleft()
            siteEvent(point)
        printBeachline(points)
    for i in range(1,len(beachline)-2):
        print(i)
        r, c, d = circle(i)
        if r:
            e = event(c, beachline[i], c[0] + d)
            beachline[i].event = e
            newEvent(e)
    while len(events) > 0:
        e = events.popleft()
        i = vertEvent(e)
        regions[points.index(e.arc.focus)].append(len(vertices) - 1)
        regions[points.index(beachline[i-1].focus)].append(len(vertices) - 1)
        regions[points.index(beachline[i].focus)].append(len(vertices) - 1)
        printBeachline(points)
    center = wrapUp1()
    for b in beachline:
        print(b.focus)
    for r in regions:
        c = np.array([vertices[v] for v in r]).mean(axis=0)
        finalregions.append(sorted(r, key=lambda p: np.arctan2(vertices[p][1] - c[1], vertices[p][0] - c[0])))
    print(finalregions)
    print(vertices)
    print(center)
    return vertices, finalregions, center
