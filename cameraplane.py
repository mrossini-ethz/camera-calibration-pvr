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

from math import pi
import mathutils
from . import algebra

def intersect_2d(pa, pb, pc, pd):
    """Find the intersection point of the lines AB and DC (2 dimensions)."""
    # Helper vectors
    ad = pd - pa
    ab = pb - pa
    cd = pd - pc
    # Solve linear system of equations s * ab + t * cd = ad
    tmp = algebra.solve_linear_system_2d(ab[0], cd[0], ad[0], ab[1], cd[1], ad[1])
    # Check for parallel lines
    if not tmp:
        return None
    s, t = tmp
    # Return the intersection point
    return pa + s * ab

def get_vanishing_point(pa, pb, pc, pd):
    """Get the vanishing point of the lines AB and DC."""
    return intersect_2d(pa, pb, pc, pd)

def get_vanishing_points(pa, pb, pc, pd):
    """Get the two vanishing points of the rectangle defined by the corners pa pb pc pd"""
    return (get_vanishing_point(pa, pb, pd, pc), get_vanishing_point(pa, pd, pb, pc))

def get_camera_plane_vector(p, scale, focal_length = 1.0):
    """Convert a 2d point in the camera plane into a 3d vector from the camera onto the camera plane"""
    # field_of_view = 2 * atan(sensor_size / 2 * focal_length), assume sensor_size = 32
    s = (16.0 / focal_length) / (scale / 2.0)
    return mathutils.Vector((p[0] * s, p[1] * s, -1.0))

def is_collinear(v1, v2):
    """Determines whether the two given vectors are collinear using a limit of 0.1 degrees"""
    limit = 0.1 * pi / 180
    # Test the angle
    return abs(v1.angle(v2)) < limit

def is_trapezoid(pa, pb, pc, pd):
    """Determines whether the polygon with the vertices pa, pb, pc, pd is a trapezoid"""
    return is_collinear(pb - pa, pc - pd) or is_collinear(pd - pa, pc - pb)

def is_trapezoid_but_not_rectangle(pa, pb, pc, pd):
    """Determines whether the polygon with the vertices pa, pb, pc, pd is a trapezoid"""
    a = is_collinear(pb - pa, pc - pd)
    b = is_collinear(pd - pa, pc - pb)
    # Exclusive OR
    return a != b

def is_to_the_right(a, b, c):
    """Checks whether the rotation angle from vector AB to vector BC is between 0 and 180 degrees when rotating to the right. Returns a number."""
    # Vector from a to b
    ab = b - a
    # Vector from b to c
    bc = c - b
    # This is a simple dot product with bc and the vector perpendicular to ab (rotated clockwise)
    return - ab[0] * bc[1] + ab[1] * bc[0]

def is_convex(pa, pb, pc, pd):
    """Checks whether the given quadrilateral corners form a convex quadrilateral."""
    # Check, which side each point is on
    to_the_right = []
    to_the_right.append(is_to_the_right(pa, pb, pc))
    to_the_right.append(is_to_the_right(pb, pc, pd))
    to_the_right.append(is_to_the_right(pc, pd, pa))
    to_the_right.append(is_to_the_right(pd, pa, pb))
    # Check whether all are on the same side
    a = True
    b = True
    for ttr in to_the_right:
        a = a and ttr > 0
        b = b and ttr < 0
    return a or b
