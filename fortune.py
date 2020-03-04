import numpy as np
import math
import random
from draw import drawFortune
from collections import deque

X0 = 0
X1 = 1
Y0 = 0
Y1 = 1
#points.sort()
#qpoints = deque(points)
events = deque([])
root = [None]

class arc1:
    def __init__(self, point, a = None, b = None):
        self.focus = point
        self.prev = a
        self.next = b
        self.event = None

#center is the point that will be added to vertices
#arcs are the arcs that form the event, as a list.
#x is the x coordinate where the vent will take place
#valid is weather is is a valid event or not. I might not even need it since I can just drop the item from the queue
class event1:
    def __init__(self):
        self.center = c
        self.arcs = a
        self.x = x
        self.valid = True

def siteEvent1(point):
    print("---")
    global root
    if root is None:
        root = arc1(point)
        return
    else:
        i = root
        while i is not None:
            print(i)
            ise, z = intersect1(point, i)
            if ise:
                if (i.next is not None and not intersect1(point, i.next)[0]):
                    i.next.prev = arc1(i.focus, i, i.next)
                    i.next = i.next.prev
                else: i.next = arc1(i.focus, i)
                i.next.prev = arc(point, i, i.next)
                i.next = i.next.prev
                return
                checkCircleEvent1(i, point[0])
                checkCircleEvent1(i.next, point[0])
                checkCircleEvent1(i.prev, point[0])
            i = i.next
    i = root
    while i.next is not None:
        i = i.next
    i.next = arc1(point)

#where does a new parabola with focus p intersect with arc i?
def intersect1(p, i):
    res = [-1, -1]
    if (i.focus[0] == p[0]): return False, res
    a = None
    b = None
    if (i.prev is not None):
        a = intersection1(i.prev.focus, i.focus, p[0])[1]
    if (i.next is not None):
        b  = intersection1(i.focus, i.next.focus, p[0])[1]
    if (i.prev is None or a <= p[1]) and (i.next is None or p[1] <= b):
        res[1] = p[1]
        res[0] = (pow(i.focus[0],2) + pow(i.focus[1] - res[1], 2) - pow(p[0], 2))/(2*i.focus[0] - 2*p[0])
        return True, res
    return False, res

#asking where two parabolas with foci p0/p1 and directrix y=l intersect?
def intersection1(p0, p1, l):
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

def checkCircleEvent1(i, x0):
    if (i.event is not None and i.event.x != x0):
        print("hi")


def fortune1(points):
    qpoints = deque(sorted(points, key=lambda p: p[0]))
    while len(qpoints) > 0:
        if (len(events) > 0 and events[0].point[0] < qpoints[0][0]):
            print("event")
#            event = events.popleft()
#            vertexEvent(event)
        else:
            point = qpoints.popleft()
            siteEvent1(point)
#    while len(events) > 0:
#        event = events.popleft()
#        vertexEvent(event)
