#!/usr/bin/python

import cairo
from math import *

def pyt(*args):
        result = 0
        for a in args:
                result += a ** 2
        return sqrt(result)

# Size of the images in pixels
size = 128

def intersect(a, b, c, d):
    # Vectors AB and CD
    ux = b[0] - a[0]
    uy = b[1] - a[1]
    vx = c[0] - d[0]
    vy = c[1] - d[1]
    det = ux * vy - uy * vx
    if det == 0:
        return None
    kx = c[0] - a[0]
    ky = c[1] - a[1]
    t1 = (kx * vy - ky * vx) / det
    t2 = (ux * ky - uy * kx) / det
    e = (a[0] + t1 * ux, a[1] + t1 * uy)
    #f = (c[0] - t2 * vx, c[1] - t2 * vy)
    #print(e, f)
    return e

def point_on(a, b, f):
    u = b[0] - a[0]
    v = b[1] - a[1]
    return (a[0] + f * u, a[1] + f * v)

def draw(filename, alpha, beta):
    with cairo.ImageSurface(cairo.Format.ARGB32, size, size) as surface:
        ctx = cairo.Context(surface)

        corners = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1), (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]

        dist = 2.7
        scale = 0.3 * size

        # Rotate cube
        for i, c in enumerate(corners):
            x = c[0] * cos(alpha) - c[2] * sin(alpha)
            y = c[1]
            z = c[0] * sin(alpha) + c[2] * cos(alpha)
            corners[i] = (x, y, z)
        # Rotate cube
        for i, c in enumerate(corners):
            x = c[0]
            y = c[1] * cos(beta) - c[2] * sin(beta)
            z = c[1] * sin(beta) + c[2] * cos(beta)
            corners[i] = (x, y, z)

        # Project cube
        corners_p = []
        for c in corners:
            #x = c[0] * dist / pyt(c[0], c[1]) / (c[2] - dist)
            #y = c[1] * dist / pyt(c[0], c[1]) / (c[2] - dist)
            x = c[0] / (c[2] - dist)
            y = c[1] / (c[2] - dist)
            corners_p.append((x * scale + size / 2, y * scale + size / 2))

        #ctx.set_source_rgb(0, 0, 0)
        #ctx.move_to(0, 0)
        #ctx.line_to(size, 0)
        #ctx.line_to(size, size)
        #ctx.line_to(0, size)
        #ctx.close_path()
        #ctx.fill()

        ctx.set_line_width(1.5)
        #ctx.set_source_rgb(86/255, 117/255, 167/255)
        ctx.set_source_rgb(0.7, 0.7, 0.7)
        i1 = intersect(corners_p[0], corners_p[1], corners_p[3], corners_p[2])
        if i1:
            ctx.set_source_rgb(0.7, 0.7, 0.7)
            ctx.move_to(*i1)
            ctx.line_to(*(point_on(corners_p[0], corners_p[1], 1.5)))
            ctx.stroke()
            ctx.move_to(*i1)
            ctx.line_to(*(point_on(corners_p[3], corners_p[2], 1.5)))
            ctx.stroke()
            ctx.move_to(*i1)
            ctx.line_to(*(point_on(corners_p[4], corners_p[5], 1.5)))
            ctx.stroke()
            ctx.move_to(*i1)
            ctx.line_to(*(point_on(corners_p[7], corners_p[6], 1.5)))
            ctx.stroke()
            #ctx.set_source_rgb(191/255, 129/255, 80/255)
            ctx.set_source_rgb(1, 1, 1)
            ctx.arc(*i1, 3, 0, 2 * pi)
            ctx.fill()
        i2 = intersect(corners_p[0], corners_p[3], corners_p[1], corners_p[2])
        if i2:
            ctx.set_source_rgb(0.7, 0.7, 0.7)
            ctx.move_to(*i2)
            ctx.line_to(*(point_on(corners_p[0], corners_p[3], 1.5)))
            ctx.stroke()
            ctx.move_to(*i2)
            ctx.line_to(*(point_on(corners_p[1], corners_p[2], 1.5)))
            ctx.stroke()
            ctx.move_to(*i2)
            ctx.line_to(*(point_on(corners_p[4], corners_p[7], 1.5)))
            ctx.stroke()
            ctx.move_to(*i2)
            ctx.line_to(*(point_on(corners_p[5], corners_p[6], 1.5)))
            ctx.stroke()
            ctx.set_source_rgb(191/255, 129/255, 80/255)
            ctx.set_source_rgb(1, 1, 1)
            ctx.arc(*i2, 3, 0, 2 * pi)
            ctx.fill()
        i3 = intersect(corners_p[0], corners_p[4], corners_p[1], corners_p[5])
        if i3:
            ctx.set_source_rgb(0.7, 0.7, 0.7)
            ctx.move_to(*i3)
            ctx.line_to(*(point_on(corners_p[0], corners_p[4], 1.5)))
            ctx.stroke()
            ctx.move_to(*i3)
            ctx.line_to(*(point_on(corners_p[1], corners_p[5], 1.5)))
            ctx.stroke()
            ctx.move_to(*i3)
            ctx.line_to(*(point_on(corners_p[2], corners_p[6], 1.5)))
            ctx.stroke()
            ctx.move_to(*i3)
            ctx.line_to(*(point_on(corners_p[3], corners_p[7], 1.5)))
            ctx.stroke()
            ctx.set_source_rgb(191/255, 129/255, 80/255)
            ctx.set_source_rgb(1, 1, 1)
            ctx.arc(*i3, 3, 0, 2 * pi)
            ctx.fill()

        ctx.set_line_width(2.5)
        ctx.set_source_rgb(1, 1, 1)
        ctx.move_to(*corners_p[0])
        ctx.line_to(*corners_p[1])
        ctx.line_to(*corners_p[2])
        ctx.line_to(*corners_p[3])
        ctx.close_path()
        ctx.stroke()

        ctx.move_to(*corners_p[4])
        ctx.line_to(*corners_p[5])
        ctx.line_to(*corners_p[6])
        ctx.line_to(*corners_p[7])
        ctx.close_path()
        ctx.stroke()

        ctx.move_to(*corners_p[0])
        ctx.line_to(*corners_p[4])
        ctx.stroke()
        ctx.move_to(*corners_p[1])
        ctx.line_to(*corners_p[5])
        ctx.stroke()
        ctx.move_to(*corners_p[2])
        ctx.line_to(*corners_p[6])
        ctx.stroke()
        ctx.move_to(*corners_p[3])
        ctx.line_to(*corners_p[7])
        ctx.stroke()

        surface.write_to_png(filename + ".png")

draw("1-point", 0, 0)
draw("2-point", radians(37), 0)
draw("3-point", radians(37), radians(35))
