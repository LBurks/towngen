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

class arc:
    def __init__(self, point, a = None, b = None):
        self.focus = point
        self.prev = a
        self.next = b
        self.lpoint = None
        self.rpoint = None
        self.event = None

#center is the point that will be added to vertices
#arcs are the arcs that form the event, as a list.
#x is the x coordinate where the event will take place
#valid is weather is is a valid event or not. I might not even need it since I can just drop the item from the queue
class event:
    def __init__(self):
        self.center = c
        self.arcs = a
        self.x = x

def siteEvent(point):
    if len(beachline) == 0:
        beachline.append(arc(point))
    else:
        for i in range(len(beachline)):
            ise, z = intersect(point, i)
            if ise:
                beachline[i+1:i+1] = [arc(point), arc(beachline[i].focus)]
                if beachline[i].event is not None:
                    events.remove(events.index(beachline[i].event))
                    beachline[i].event = None
                print("Check for events that will remove the new arcs")
                circle(beachline[i-1], beachline[i], beachline[i+1])

#where does a new parabola with focus p intersect with arc i?
def intersect(p, i):
    larc = beachline[i-1] if i > 0 else None
    carc = beachline[i]
    rarc = beachline[i+1] if i > len(beachline)-1 else None
    
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

def circle(larc, carc, rarc):
    if carc.focus[0] >= larc.focus[0] and carc.focus[0] >= rarc.focus[0]:
        print("center arc is the furthest point out. This circle won't form a vertex")
        return
    print("Get the perpendicular bisectors for both lines.")
    print("If the lines are parallel, there is not a circle, return some failure")
    print("If they aren't find the intersect, this is your circle's midpoint.")
    print("Add the result of the distance formula to the x value of the circle's center. That is your site for the event.")

def checkCircleEvent(i, x0):
    if (i.event is not None and i.event.x != x0):
        print("hi")


def fortune(points):
    qpoints = deque(sorted(points, key=lambda p: p[0]))
    while len(qpoints) > 0:
        if (len(events) > 0 and events[0].point[0] < qpoints[0][0]):
            print("event")
        else:
            point = qpoints.popleft()
            siteEvent(point)
        print(beachline)
