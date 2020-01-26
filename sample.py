import math

def area_for_polygon(polygon):
    result = 0
    imax = len(polygon) - 1
    for i in range(0,imax):
        result += (polygon[i]['x'] * polygon[i+1]['y']) - (polygon[i+1]['x'] * polygon[i]['y'])
    result += (polygon[imax]['x'] * polygon[0]['y']) - (polygon[0]['x'] * polygon[imax]['y'])
    return result / 2.

def centroid_for_polygon(polygon):
    area = area_for_polygon(polygon)
    imax = len(polygon) - 1

    result_x = 0
    result_y = 0
    for i in range(0,imax):
        result_x += (polygon[i]['x'] + polygon[i+1]['x']) * ((polygon[i]['x'] * polygon[i+1]['y']) - (polygon[i+1]['x'] * polygon[i]['y']))
        result_y += (polygon[i]['y'] + polygon[i+1]['y']) * ((polygon[i]['x'] * polygon[i+1]['y']) - (polygon[i+1]['x'] * polygon[i]['y']))
    result_x += (polygon[imax]['x'] + polygon[0]['x']) * ((polygon[imax]['x'] * polygon[0]['y']) - (polygon[0]['x'] * polygon[imax]['y']))
    result_y += (polygon[imax]['y'] + polygon[0]['y']) * ((polygon[imax]['x'] * polygon[0]['y']) - (polygon[0]['x'] * polygon[imax]['y']))
    result_x /= (area * 6.0)
    result_y /= (area * 6.0)

    return {'x': result_x, 'y': result_y}

def bottommost_index_for_polygon(polygon):
    bottommost_index = 0
    for index, point in enumerate(polygon):
        if (point['y'] < polygon[bottommost_index]['y']):
            bottommost_index = index
    return bottommost_index

def angle_for_vector(start_point, end_point):
    y = end_point['y'] - start_point['y']
    x = end_point['x'] - start_point['x']
    angle = 0

    if (x == 0):
        if (y > 0):
            angle = 90.0
        else:
            angle = 270.0
    elif (y == 0):
        if (x > 0):
            angle = 0.0
        else:
            angle = 180.0
    else:
        angle = math.degrees(math.atan((y+0.0)/x))
        if (x < 0):
            angle += 180
        elif (y < 0):
            angle += 360

    return angle

def convex_hull_for_polygon(polygon):
    starting_point_index = bottommost_index_for_polygon(polygon)
    convex_hull = [polygon[starting_point_index]]
    polygon_length = len(polygon)

    hull_index_candidate = 0 #arbitrary
    previous_hull_index_candidate = starting_point_index
    previous_angle = 0
    while True:
        smallest_angle = 360

        for j in range(0,polygon_length):
            if (previous_hull_index_candidate == j):
                continue
            current_angle = angle_for_vector(polygon[previous_hull_index_candidate], polygon[j])
            if (current_angle < smallest_angle and current_angle > previous_angle):
                hull_index_candidate = j
                smallest_angle = current_angle

        if (hull_index_candidate == starting_point_index): # we've wrapped all the way around
            break
        else:
            convex_hull.append(polygon[hull_index_candidate])
            previous_angle = smallest_angle
            previous_hull_index_candidate = hull_index_candidate

    return convex_hull
