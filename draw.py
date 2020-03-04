import cairo
import numpy

#Draws the Voronoi diagram and Delauny triangulation. Unused and kept as reference
def drawOutlines(diagram, outfile = "outline", width = 1000.0, height = 1000.0):
    surface = cairo.SVGSurface(outfile + ".svg", width, height)
    context = cairo.Context(surface)
    context.scale(width, height)
    context.set_line_width(0.001)
    context.set_font_size(0.005)
    context.set_source_rgba(1, 1, 1, 1)
    context.rectangle(0, 0, width, height)
    context.fill()
    context.stroke()

    for r in diagram.regions:
        context.set_source_rgba(0, 0, 0, 1)
        context.move_to(diagram.vertices[r[0]].coordinates[0], diagram.vertices[r[0]].coordinates[1])
        for v in [diagram.vertices[x] for x in r[1:]]:
            context.line_to(v.coordinates[0], v.coordinates[1])
        context.close_path()
        context.stroke()
    for r in diagram.ridges:
        ridges = diagram.ridges[r]
        for p2, v1, v2 in ridges:
            context.move_to(diagram.points[r].coordinates[0], diagram.points[r].coordinates[1])
            context.line_to(diagram.points[p2].coordinates[0], diagram.points[p2].coordinates[1])
            context.stroke()
    context.set_source_rgb(0.7, 0.2, 0.2)
    for p in diagram.points:
        context.move_to(0, 0)
        context.arc(p.coordinates[0], p.coordinates[1], 0.004, 0, 2*numpy.pi)
        context.fill()
        context.stroke()
    context.set_source_rgb(0.2, 0.2, 0.7)
    for p in diagram.vertices:
        context.move_to(0, 0)
        context.arc(p.coordinates[0], p.coordinates[1], 0.004, 0, 2*numpy.pi)
        context.fill()
        context.stroke()
    surface.write_to_png(outfile + ".png")


#Draws each polygon and colors them green for land, blue for water.
def drawPolys(diagram, outfile = "polys", width = 1000.0, height = 1000.0):
    print("Drawing map")
    surface = cairo.SVGSurface(outfile + ".svg", width, height)
    context = cairo.Context(surface)
    context.scale(width, height)
    context.set_line_width(0.001)
    context.set_source_rgba(1, 0, 0, 0)
    res = []
    for r, p in zip(diagram.regions, diagram.points):
        context.move_to(diagram.vertices[r[0]].coordinates[0], diagram.vertices[r[0]].coordinates[1])
        # for v in diagram.vertices[r[1:]]:
        for v in [diagram.vertices[x] for x in r[1:]]:
            context.line_to(v.coordinates[0], v.coordinates[1])
        context.close_path()
        if p.land:
            context.set_source_rgba(0, 0.8, 0, 1)
        else:
            if p.ocean:
                context.set_source_rgba(0, 0, 0.5, 1)
            else:
                context.set_source_rgba(0, 0, 1, 1)
        context.fill_preserve()
        context.stroke()
        for y in [diagram.vertices[x] for x in r if diagram.vertices[x].border]:
            context.set_source_rgba(1, 0, 0, 1)
            context.move_to(0, 0)
            context.arc(y.coordinates[0], y.coordinates[1], 0.004, 0, 2*numpy.pi)
            context.fill()
            context.stroke()
    surface.write_to_png(outfile + ".png")


def drawFortune(points, lines, outfile = "fortune", width = 1000.0, height = 1000.0):
#    print("Drawing map")
    surface = cairo.SVGSurface(outfile + ".svg", width, height)
    context = cairo.Context(surface)
    context.scale(width, height)
    context.set_line_width(0.001)
    context.set_font_size(0.005)
    context.set_source_rgba(1, 1, 1, 1)
    context.rectangle(0, 0, width, height)
    context.fill()
    context.stroke()
    context.set_source_rgb(0.7, 0.2, 0.2)
    for p in points:
        context.move_to(0, 0)
        context.arc(p[0], p[1], 0.004, 0, 2*numpy.pi)
        context.fill()
        context.stroke()
    for l in lines:
        p1 = l.start
        p2 = l.end
        context.set_source_rgb(0.2, 0.2, 0.7)
        if (p1 is not None):
            context.move_to(0, 0)
            context.arc(p1[0], p1[1], 0.004, 0, 2*numpy.pi)
            context.fill()
        if (p2 is not None):
            context.move_to(0, 0)
            context.arc(p2[0], p2[1], 0.004, 0, 2*numpy.pi)
            context.fill()
        context.stroke()
        context.set_source_rgb(0, 0, 0)
        if (p1 is not None and p2 is not None):
            context.move_to(p1[1], p1[0])
            context.line_to(p2[1], p2[0])
            context.stroke()
    surface.write_to_png(outfile + ".png")
