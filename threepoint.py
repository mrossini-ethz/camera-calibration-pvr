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

import bpy
import mathutils
from math import sqrt, pi, degrees

from . import algebra
from . import cameraplane
from . import solverectangle

def calculate_focal_length(pa, pb, pc, pd, scale):
    """Get the vanishing points of the rectangle as defined by pa, pb, pc and pd"""
    pm, pn = cameraplane.get_vanishing_points(pa, pb, pc, pd)
    # Calculate the vectors from camera to the camera plane where the vanishing points are located
    vm = cameraplane.get_camera_plane_vector(pm, scale)
    vn = cameraplane.get_camera_plane_vector(pn, scale)
    # Calculate the focal length
    return sqrt(abs(vm.dot(vn)))

def calculate_focal_length_shifted(vertices, attached_vertices, dangling_vertices, scale):
    # Get the 3 vanishing points
    v1 = cameraplane.get_vanishing_point(vertices[0], vertices[3], vertices[1], vertices[2])
    v2 = cameraplane.get_vanishing_point(attached_vertices[0], dangling_vertices[0], attached_vertices[1], dangling_vertices[1])
    v3 = cameraplane.get_vanishing_point(vertices[0], vertices[1], vertices[3], vertices[2])
    print("Vanishing points:", v1, v2, v3)
    # Use that v1, v2, and v3 are all perpendicular to each other to solve for the lens shift
    # x*(v2x - v3x) + y*(v2y - v3y) = v1y*v2y - v1y*v3y + v1x*v2x - v1x*v3x
    # x*(v1x - v3x) + y*(v1y - v3y) = v1y*v2y - v2y*v3y + v1x*v2x - v2x*v3x
    A1 = v2[0] - v3[0]
    B1 = v2[1] - v3[1]
    C1 = v1[1] * v2[1] - v1[1] * v3[1] + v1[0] * v2[0] - v1[0] * v3[0]
    A2 = v1[0] - v3[0]
    B2 = v1[1] - v3[1]
    C2 = v1[1] * v2[1] - v2[1] * v3[1] + v1[0] * v2[0] - v2[0] * v3[0]
    # Solve A1 * x + B1 * y = C1, A2 * x + B2 * y = C2
    shift_x, shift_y = algebra.solve_linear_system_2d(A1, B1, C1, A2, B2, C2)
    optical_centre = mathutils.Vector((shift_x, shift_y))
    shift_x /= -scale
    shift_y /= -scale
    print("Shift:", shift_x, shift_y)
    # Get the focal length
    vm = cameraplane.get_camera_plane_vector(v1 - optical_centre, scale)
    vn = cameraplane.get_camera_plane_vector(v3 - optical_centre, scale)
    # Calculate the focal length
    focal = sqrt(abs(vm.dot(vn)))
    print("Focal:", focal)
    return focal, shift_x, shift_y, optical_centre

def calibrate_camera(pa, pb, pc, pd, scale):
    # Calculate the focal length of the camera
    focal = calculate_focal_length(pa, pb, pc, pd, scale)
    # Reconstruct the rectangle using this focal length
    return (focal,) + solverectangle.reconstruct_rectangle(pa, pb, pc, pd, scale, focal)

def calibrate_camera_shifted(vertices, attached_vertices, dangling_vertices, scale):
    # Get the focal length, the optical centre and the vertical shift of the camera
    focal, shift_x, shift_y, optical_centre = calculate_focal_length_shifted(vertices, attached_vertices, dangling_vertices, scale)
    # Correct for the camera shift
    for i in range(len(vertices)):
        vertices[i] -= optical_centre
    # Reconstruct the rectangle using the focal length and return the results, together with the shift values
    return (focal,) + solverectangle.reconstruct_rectangle(vertices[0], vertices[1], vertices[2], vertices[3], scale, focal) + (shift_x, shift_y)
