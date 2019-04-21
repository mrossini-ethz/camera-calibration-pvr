# Blender Plugin: Camera Calibration with Perspective Views of Rectangles
# Copyright (C) 2017  Marco Rossini
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# version 2 as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# This Blender plugin is based on the research paper "Recovery of Intrinsic
# and Extrinsic Camera Parameters Using Perspective Views of Rectangles" by
# T. N. Tan, G. D. Sullivan and K. D. Baker, Department of Computer Science,
# The University of Reading, Berkshire RG6 6AY, UK, Email: T.Tan@reading.ac.uk,
# from the Proceedings of the British Machine Vision Conference, published by
# the BMVA Press.

import mathutils
from math import pi, sin, cos, atan2
from . import cameraplane

def calibrate_camera(pa, pb, pc, pd, scale, focal, W, L):
    # Ensure that AB and CD are parallel
    if not cameraplane.is_collinear(pb - pa, pc - pd):
        # Change the order of the vertices so that AB and CD are parallel
        pa, pb, pc, pd = pb, pc, pd, pa

    # Ensure that AB is longer than CD
    AB = (pb - pa).length
    CD = (pd - pc).length
    if CD > AB:
        # Change the order of the vertices to make AB longer than CD
        pa, pb, pc, pd = pc, pd, pa, pb
        AB, CD = CD, AB
    vAB = pb - pa

    # Calculate the angle of the parallels with respect to the image
    alpha = atan2(vAB[1], vAB[0])
    # Ensure the angle is between -180 and 180 degrees
    if alpha > pi/2:
        alpha -= pi
    elif alpha < -pi/2:
        alpha += pi

    # Decide plane orientation
    is_horizontal = abs(alpha) <= pi/4

    # Determine the camera rotation angle
    if is_horizontal:
        cam_rot = alpha
    elif alpha > pi/4:
        cam_rot = alpha - pi/2
    elif alpha < -pi/4:
        cam_rot = alpha + pi/2

    # Calculate the vanishing point
    v1 = cameraplane.get_vanishing_point(pa, pd, pb, pc)

    # Midpoints of trapezoid base and top
    mA = (pa + pb) / 2
    mC = (pd + pc) / 2

    # "Up" direction
    vup = mathutils.Vector((-sin(alpha), cos(alpha)))
    # Determine the z-offset direction
    if vup.dot(pc - pa) > 0:
        # Bottom situation
        zsign = 1
    else:
        # Top situation
        zsign = -1
    # Special situation
    if not is_horizontal and alpha > 0:
        zsign = -zsign

    # Trapezoid height
    vn = mathutils.Vector((vAB[1], -vAB[0]))
    H = abs(vn.dot(pc - pa) / vn.length)

    # Set camera shift parameters
    shift_x = -v1[0] / scale
    shift_y = -v1[1] / scale

    if focal:
        # Focal length is given. Calculate rectangle length.
        y = focal * W * scale / AB / 32
        L = y * (AB / CD - 1)
    else:
        # Rectangle length is given. Calculate focal length (assume sensor size 32 mm).
        y = L / (AB / CD - 1)
        focal = AB * 32 * y / W / scale

    # Calculate camera x- and z-positions
    z = zsign * H * W / AB * (y + L) / L
    x = -abs(vAB.normalized().dot(mA - v1)) / scale * y * 32 / focal

    if is_horizontal:
        coords = (mathutils.Vector((-W/2, -L/2, 0)), mathutils.Vector((W/2, -L/2, 0)), mathutils.Vector((W/2, L/2, 0)), mathutils.Vector((-W/2, L/2, 0)))
    else:
        x, z = z, x
        coords = (mathutils.Vector((0, -L/2, -W/2)), mathutils.Vector((0, -L/2, W/2)), mathutils.Vector((0, L/2, W/2)), mathutils.Vector((0, L/2, -W/2)))

    # Reconstruct the rectangle using the focal length and return the results, together with the shift value
    return (focal, mathutils.Vector((x, -L / 2 - y, z)), cam_rot, shift_x, shift_y, coords, L)
