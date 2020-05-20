import numpy as np
import random
from math import sqrt
from draw import drawFortune
from collections import deque

bounding_box = np.array([])
points = np.array([])
beachline = []

#focus: The node that is the focus of the arc
#event: The event that will remove this arc
class arc:
    def __init__(self, f, e = None):
        self.focus = f
        self.event = e

#center: The vertex this event will represent
#arc: The arc this event will remove
#x: The x coordinate the event will take place at
#valid: Weather the event is valid or not.
class event:
    def __init__(self, c, a, x):
        self.center = c
        self.arc = a
        self.x = x
        self.valid = True

def vertEvent(event):
    print("Vertex Event")
    if not event.valid: return

def siteEvent(point):
    #1. Add the root of the beachline if it doesn't exist. We're done
    if len(beachline) == 0:
        beachline.append(arc(point))
        return
    #2. Find where a line drawn from the new site intersects the beachline
    for i in range(len(beachline)):
        ise, z = intersect(i, point)
        if ise:
            # new parabola intersects arc i. duplicate i if needed.
            if (i < len(beachline) - 1 and not intersect(i+1, point)[0]):
                beachline[i+1:i+1] = [arc(point), arc(beachline[i].focus)]
            else: beachline.append(arc(point))
            i += 1 #i now points to the new arc.

            if i in range(1, len(beachline) - 1):
                check_event(i-1, point)
                check_event(i, point)
                check_event(i+1, point)
            return

def intersect(i, p):
    focus = points[beachline[i].focus]
    point = points[p]
    res = [-1, -1]

    #if the focus and point both have the same x value, they won't intersect
    if focus[0] == point[0]: return False, res

    #test where the bounds of the center arc is.
    a, b = None, None
    if i > 0:
        a = intersection(i-1, i, point[0])
    if i < len(beachline) - 1:
        b = intersection(i, i+1, point[0])

    if (a is None or a[1] <= point[1]) and (b is None or point[1] <= b[1]):
        #Actually find the intersect, y is the y-value for our site. 
        #x-value is derived from parabola equation
        res[1] = point[1]    
        res[0] = (pow(focus[0], 2) + pow(focus[1] - res[1], 2) - pow(point[0], 2)) / (2 * focus[0] - 2 * point[0])

        return True, res
    return False, res

def intersection(a0, a1, l):
    f0 = points[beachline[a0].focus]
    f1 = points[beachline[a1].focus]
    p = f0
    res =[0, 0]
    #The points are on the same vertical line. The y value will be in the middle of them
    if f0[0] == f1[0]: res[1] = (f0[1] + f1[1]) / 2
    #The lower point is on the directrix. the y value will be the same
    elif f1[0] == l: res[1] = f1[1]
    #the upper point is on the directrix. The y value will be the same. 
    #The upper point also will not form a parabola, so our parabola focus needs to change.
    elif f0[0] == l: 
        res[1] = f0[1]
        p = f1
    else:
        #Apply the quadratic fomula
        z0 = 2 * (f0[0] - l)
        z1 = 2 * (f1[0] - l)

        a = (1/z0 - 1/z1)
        b = -2 * (f0[1]/z0 - f1[1]/z1 )
        c = (pow(f0[0], 2) + pow(f0[1], 2) - pow(l, 2)) / z0  - (pow(f1[0], 2) + pow(f1[1], 2) - pow(l, 2)) / z1

        res[1] = (-b - sqrt(pow(b, 2) - 4 * a * c)) / (2*a)

    res[0] = (pow(p[0], 2) + pow(p[1]-res[1], 2) - pow(l, 2)) /(2 * p[0] - 2 * l)
    return res

def check_event(i, p):
    a = beachline[i]
    x0 = points[p][0]
    if (a.event is not None and a.event.x != x0): 
        a.event.valid = False
    a.event = None
    if i == 0 or i == len(beachline) - 1:
        res, c, x = circle(i)
        if res:
            e = event(c, a, x)
            a.event = e
            events.push(e)

def circle(i):
    #a, b, c are the arcs' foci. left, center, and right respectively
    a = points[beachline[i-1].focus]
    b = points[beachline[i].focus]
    c = points[beachline[i+1].focus]
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
    x = distance + xc
    print("returns an event")
    return True, [xc, yc], x

def circle1(i):
    a = points[beachline[i-1].focus]
    b = points[beachline[i].focus]
    c = points[beachline[i+1].focus]
    #center = the intersection of the two perpendicular bisectors
    center = [0, 0]
    #x = the x value the event will happen at.
    x = 0
    #Should return this if all three points are on the same line (no perpendicular bisector) or the midpont is the outermost point (The center doesn't actually point where the center arc terminates)
    det = a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])
    if det <= 0 or (a[0] < b[0] and c[0] < b[0]): return
    #find the midpoints
    mid1 = [(a[0] + b[0])/2, (a[1]+b[1])/2]
    mid2 = [(b[0] + c[0])/2, (b[1]+c[1])/2]
    #slope of the perpendicular bisector
    m1 = (b[0] - a[0])/(b[1] - a[1])
    m2 = (c[0] - b[0])/(c[1] - b[1])
    #I believe this will get me my center
    center[0] = (m1*mid1[0] - m2*mid2[0] - mid1[1] + mid2[1]) / (m1 - m2)
    center[1] = m1 * (center[0] - mid1[0]) + mid1[1]
    distance = sqrt(pow(center[0] - a[0], 2) + pow(center[1] - a[1], 2))
    x = center[0] + distance
    print('here')
    e = event(center, beachline[i], x)
    return

def fortune(p, bb):
    global bounding_box, points
    bounding_box = bb
    points = p
    queue = deque(sorted([a for a in range(0, len(p))], key=lambda a: p[a][0]))
    events = deque([])
    #1. process sites and vertices concurrently.
    while len(queue) > 0:
        if (len(events) > 0 and events[0].x < points[queue[0]][0]): vertEvent(events.popleft())
        else: siteEvent(queue.popleft())
    #2. Determine events for beachline after all sites have been processed.
    #3. Process each event.
    while len(events) > 0:
        vertEvent()
