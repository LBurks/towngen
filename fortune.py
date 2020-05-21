from collections import deque

beachline = []

#focus: The node that is the focus of the arc
#event: The event that will remove this arc
class arc:
    def __init__(self, f, e = None):
        self.focus = f
        self.event = e

def siteEvent(point):
    #Add a root node
    if len(beachline) == 0:
        beachline.append(arc(point))
        return
    #find where on our beachline the new node fits
    #find intersection of line and parabola
    #sqrt(pow(i.x - p.x,2)) = sqrt((pow(f.x-i.x,2) + pow(f.y-1.y,2))
    #attach arc to the end, since there is no intersection
    beachline.append(arc(point))

def fortune(points, bounding_box):
    queue = deque(sorted([a for a in range(0, len(points))], key=lambda a: points[a][0]))
    events = deque([])
    while len(queue) > 0:
        print("Add a node to the beachline")
        queue.popleft()

