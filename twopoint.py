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
from math import sqrt

from . import cameraplane
from . import solverectangle

def calculate_focal_length(pa, pb, pc, pd, pe, pf, scale):
    # Determine which two edges of the polygon ABCD are parallel, reorder if necessary
    if cameraplane.is_collinear(pb - pa, pc - pd):
        # rename vertices to make AD and BC parallel
        tmp = [pa, pb, pc, pd]
        pa = tmp[1]
        pb = tmp[2]
        pc = tmp[3]
        pd = tmp[0]
    # Get the horizon direction vector
    vertical = pd - pa + pc - pb
    horizon = mathutils.Vector((-vertical[1], vertical[0]))
    # FIXME: remove debugging printouts in this function
    print("horizon", horizon)
    # Determine the vanishing point of the polygon ABCD
    vanish1 = cameraplane.get_vanishing_point(pa, pb, pc, pd)
    print("vanish1", vanish1)
    # Intersect the dangling edge with the horizon to find the second vanishing point
    vanish2 = cameraplane.get_vanishing_point(pe, pf, vanish1, vanish1 + horizon)
    print("vanish2", vanish2)
    # Find the rotation point
    # FIXME: don't use the x-coordinate directly
    t = -vanish1[0] / horizon[0]
    optical_centre = vanish1 + t * horizon
    # Get the camera shift
    shift = -optical_centre[1] / scale
    print("shift", shift)
    # Find the focal length
    dist = sqrt((vanish1 - optical_centre).length * (vanish2 - optical_centre).length)
    # Assume sensor size of 32
    focal = dist / (scale / 2.) * 16
    print("focal", focal)
    return (focal, optical_centre, shift)

def calibrate_camera(pa, pb, pc, pd, pe, pf, scale):
    # Get the focal length, the optical centre and the vertical shift of the camera
    focal, optical_centre, shift = calculate_focal_length(pa, pb, pc, pd, pe, pf, scale)
    # Correct for the camera shift
    pa = pa - optical_centre
    pb = pb - optical_centre
    pc = pc - optical_centre
    pd = pd - optical_centre
    # Reconstruct the rectangle using the focal length and return the results, together with the shift value
    return (focal,) + solverectangle.reconstruct_rectangle(pa, pb, pc, pd, scale, focal) + (shift,)
