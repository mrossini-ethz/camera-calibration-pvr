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

from math import atan2
import mathutils

def get_transformation(ra, rb, rc, rd):
    """Average the vectors AD, BC and AB, DC and normalize them"""
    ex = (rb - ra + rc - rd).normalized()
    ey = (rd - ra + rc - rb).normalized()
    # Get the unit vector in z-direction by using the cross product
    # Normalize, because rx and ry may not be perfectly perpendicular
    ez = ex.cross(ey).normalized()
    return [ex, ey, ez, (ra + rb + rc + rd) / 4.0]

def get_rot_angles(ex, ey, ez):
    """Get the x- and y-rotation from the ez unit vector"""
    rx = atan2(ez[1], ez[2])
    rx_matrix = mathutils.Euler((rx, 0.0, 0.0), "XYZ")
    # Rotate the ez vector by the previously found angle
    ez.rotate(rx_matrix)
    # Negative value because of right handed rotation
    ry = - atan2(ez[0], ez[2])
    # Rotate the ex vector by the previously found angles
    rxy_matrix = mathutils.Euler((rx, ry, 0.0), "XYZ")
    ex.rotate(rxy_matrix)
    # Negative value because of right handed rotation
    rz = - atan2(ex[1], ex[0])
    return [rx, ry, rz]

def apply_transformation(vertices, translation, rotation):
    n = len(vertices)
    result = []
    for i in range(len(vertices)):
        result.append(vertices[i] - translation)
        result[-1].rotate(rotation)
    return result

def vertex_apply_transformation(p, scale, rotation, translation):
    # Make a copy of the vertex
    p = p.copy()
    # Apply the scale
    for i in range(3):
        p[i] *= scale[i]
    # Apply rotation
    p.rotate(rotation)
    # Apply translation and project to x-y-plane
    p = p + translation
    return p
