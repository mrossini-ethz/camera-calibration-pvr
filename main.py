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
from math import sqrt, pi, atan2, degrees, sin, cos
from sys import float_info

from . import polynomial
from . import rootfinder
from . import algebra
from . import cameraplane
from . import transformation
from . import reference
from . import scene

### Algorithms for intrinsic camera parameters ###################################

# The algorithm solves the following intrinsic parameters of the camera:
# - (F) focal length
# - (X) lens shift in x
# - (Y) lens shift in y
# - (A) pixel aspect ratio
# There are the following solutions for the extrinsic parameters:
# - (P) camera position
# - (R) camera rotation
# The following inputs are accepted:
# - (S) single rectangle
# - (V) single rectangle with dangling vertex
# - (D) dual rectangle
# Naming convention for the different solvers (examples):
# - solve_F_S(): solves for the focal length using a single rectangle as input
# - solve_FY_V(): solves for the focal length and y-shift from a single rectangle with dangling vertex
# - solve_FXY_VV(): solves for the focal length and x- and y- shift from a single rectangle with two dangling vertices
# - calibrate_camera_F_PR_S(): calibrates focal length, position and rotation from a single rectangle
# - calibrate_camera_FX_PR_V(): calibrates focal length, x-shift, position and rotation from a single rectangle with dangling vertex
# - calibrate_camera_FXY_PR_VV(): calibrates focal length, x- and y-shift, position and rotation from a single rectangle with two dangling vertices
# - calibrate_camera_FXY_P_S(): calibrates focal length, x- and y-shift and position from a single rectangle in 1-point perspective

def solve_F_S(pa, pb, pc, pd, scale):
    """Get the vanishing points of the rectangle as defined by pa, pb, pc and pd"""
    pm, pn = cameraplane.get_vanishing_points(pa, pb, pc, pd)
    # Calculate the vectors from camera to the camera plane where the vanishing points are located
    vm = cameraplane.get_camera_plane_vector(pm, scale)
    vn = cameraplane.get_camera_plane_vector(pn, scale)
    # Calculate the focal length
    return sqrt(abs(vm.dot(vn)))

def solve_FY_V(pa, pb, pc, pd, pe, pf, scale):
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

def solve_FXY_VV(vertices, attached_vertices, dangling_vertices, scale):
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

### Helper functions for extrinsic camera parameter algorithms ###################

def get_lambda_d_poly_a(qab, qac, qad, qbc, qbd, qcd):
    """Equation A (see paper)"""
    d4 = qac * qbd ** 2 - qad * qbc * qbd
    d3 = qab * qad * qbc + qad ** 2 * qbc * qbd + qad ** 2 * qcd + qbc * qbd - 2 * qab * qac * qbd - qab * qad * qbd * qcd - qac * qad * qbd ** 2
    d2 = qab ** 2 * qac + qab ** 2 * qad * qcd + 3 * qab * qac * qad * qbd + qab * qbd * qcd - qab * qad ** 2 * qbc - qab * qbc - qac * qad ** 2 - qad * qbc * qbd - 2 * qad * qcd
    d1 = qab * qad * qbc + 2 * qac * qad + qcd - 2 * qab ** 2 * qac * qad - qab ** 2 * qcd - qab * qac * qbd
    d0 = qab ** 2 * qac - qac
    return polynomial.make_poly([d0, d1, d2, d3, d4])

def get_lambda_d_poly_b(qab, qac, qad, qbc, qbd, qcd):
    """Equation B (see paper)"""
    d4 = qbd - qbd * qcd ** 2
    d3 = qab * qcd ** 2 + qac * qbd * qcd + 2 * qad * qbd * qcd ** 2 - qab - 2 * qad * qbd - qad * qbc * qcd
    d2 = 2 * qab * qad + qac * qad * qbc + qad ** 2 * qbc * qcd + qad **2 * qbd + qbc * qcd - qab * qac * qcd - qab * qad * qcd ** 2 - 3 * qac * qad * qbd * qcd - qbd * qcd ** 2
    d1 = qab * qac * qad * qcd + qac ** 2 * qad * qbd + 2 * qac * qbd * qcd - qab * qad ** 2 - qac * qad ** 2 * qbc - qac * qbc - qad * qbc * qcd
    d0 = qac * qad * qbc - qac ** 2 * qbd
    return polynomial.make_poly([d0, d1, d2, d3, d4])

def get_lambda_d(pa, pb, pc, pd, scale, focal_length):
    """Calculate the vectors from camera to the camera plane where the rectangle corners are located"""
    va = cameraplane.get_camera_plane_vector(pa, scale, focal_length).normalized()
    vb = cameraplane.get_camera_plane_vector(pb, scale, focal_length).normalized()
    vc = cameraplane.get_camera_plane_vector(pc, scale, focal_length).normalized()
    vd = cameraplane.get_camera_plane_vector(pd, scale, focal_length).normalized()
    # Calculate dot products
    qab = va.dot(vb)
    qac = va.dot(vc)
    qad = va.dot(vd)
    qbc = vb.dot(vc)
    qbd = vb.dot(vd)
    qcd = vc.dot(vd)
    # Determine the equation that needs to be solved
    pa = polynomial.norm(get_lambda_d_poly_a(qab, qac, qad, qbc, qbd, qcd))
    pb = polynomial.norm(get_lambda_d_poly_b(qab, qac, qad, qbc, qbd, qcd))
    print("A:", pa)
    print("B:", pb)
    p = polynomial.reduce(polynomial.sub(pa, pb))
    print("P:", p)
    # Solve the equation
    roots = rootfinder.find_poly_roots(p)
    print("Solutions:")
    # Iterate over all roots
    solutions = []
    for ld in roots:
        # Calculate the other parameters
        #ld = 1.10201
        lb = (qad * ld - 1) / (qbd * ld - qab)
        lc = (qad * ld - ld ** 2) / (qac - qcd * ld)
        # Scale the vectors pointing to the corners from the camera plane to 3d space
        ra = va
        rb = vb * lb
        rc = vc * lc
        rd = vd * ld
        # Printout for debugging
        print("x:", ld)
        # Corner angles
        angles = [degrees((rb - ra).angle(rd - ra)), degrees((ra - rb).angle(rc - rb)), degrees((rb - rc).angle(rd - rc)), degrees((rc - rd).angle(ra - rd))]
        print("Corner angles:", angles)
        # Rectangle size
        width = (rb - ra).length
        height = (rd - ra).length
        # Flatness (normal distance of point rd to plane defined by ra, rb, rc
        n = (ra - rb).cross(rc - rb)
        d = n.dot(ra)
        dist = abs(n.dot(rd) - d) / n.length
        print("Flatness:", dist, "=", dist / max(width, height) * 100, "%")
        # Calculate badness
        badness = 0.0
        # FIXME: angle badness and flatness badness should be weighted somehow
        for ang in angles:
            badness += abs(ang - 90)
        badness += abs(dist / max(width, height) * 100)
        print("Badness:", badness)
        solutions.append((badness, [ra, rb, rc, rd]))
    # Chose solution with best score
    best_badness = solutions[0][0]
    best_index = 0
    for i in range(1, len(solutions)):
        if best_badness > solutions[i][0]:
            best_index = i
            best_badness = solutions[i][0]
    # Return the best solution
    return solutions[best_index][1]

### Algorithms for extrinsic camera parameters ###################################

def reconstruct_rectangle(pa, pb, pc, pd, scale, focal):
    # Calculate the coordinates of the rectangle in 3d
    coords = get_lambda_d(pa, pb, pc, pd, scale, focal)
    # Calculate the transformation of the rectangle
    trafo = transformation.get_transformation(coords[0], coords[1], coords[2], coords[3])
    # Reconstruct the rotation angles of the transformation
    angles = transformation.get_rot_angles(trafo[0], trafo[1], trafo[2])
    xyz_matrix = mathutils.Euler((angles[0], angles[1], angles[2]), "XYZ")
    # Reconstruct the camera position and the corners of the rectangle in 3d such that it lies on the xy-plane
    tr = trafo[-1]
    cam_pos = transformation.apply_transformation([mathutils.Vector((0.0, 0.0, 0.0))], tr, xyz_matrix)[0]
    corners = transformation.apply_transformation(coords, tr, xyz_matrix)
    # Printout for debugging
    print("Focal length:", focal)
    print("Camera rotation:", degrees(angles[0]), degrees(angles[1]), degrees(angles[2]))
    print("Camera position:", cam_pos)
    length = (coords[0] - coords[1]).length
    width = (coords[0] - coords[3]).length
    size = max(length, width)
    print("Rectangle length:", length)
    print("Rectangle width:", width)
    print("Rectangle corners:", corners)
    return (cam_pos, xyz_matrix, corners, size)

def calibrate_camera_F_PR_S(pa, pb, pc, pd, scale):
    # Calculate the focal length of the camera
    focal = solve_F_S(pa, pb, pc, pd, scale)
    # Reconstruct the rectangle using this focal length
    return (focal,) + reconstruct_rectangle(pa, pb, pc, pd, scale, focal)

def calibrate_camera_FX_PR_V(pa, pb, pc, pd, pe, pf, scale):
    # Get the focal length, the optical centre and the vertical shift of the camera
    focal, optical_centre, shift = solve_FY_V(pa, pb, pc, pd, pe, pf, scale)
    # Correct for the camera shift
    pa = pa - optical_centre
    pb = pb - optical_centre
    pc = pc - optical_centre
    pd = pd - optical_centre
    # Reconstruct the rectangle using the focal length and return the results, together with the shift value
    return (focal,) + reconstruct_rectangle(pa, pb, pc, pd, scale, focal) + (shift,)

def calibrate_camera_FXY_PR_VV(vertices, attached_vertices, dangling_vertices, scale):
    # Get the focal length, the optical centre and the vertical shift of the camera
    focal, shift_x, shift_y, optical_centre = solve_FXY_VV(vertices, attached_vertices, dangling_vertices, scale)
    # Correct for the camera shift
    for i in range(len(vertices)):
        vertices[i] -= optical_centre
    # Reconstruct the rectangle using the focal length and return the results, together with the shift values
    return (focal,) + reconstruct_rectangle(vertices[0], vertices[1], vertices[2], vertices[3], scale, focal) + (shift_x, shift_y)

def calibrate_camera_FXY_P_S(pa, pb, pc, pd, scale, focal, W, L):
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

### Operator F PR S ##############################################################

class CameraCalibration_F_PR_S_Operator(bpy.types.Operator):
    """Calibrates the focal length, position and rotation of the active camera."""
    bl_idname = "camera.camera_calibration_f_pr_s"
    bl_label = "Solve Focal"
    bl_options = {"REGISTER", "UNDO"}

    # Properties
    vertical_property : bpy.props.BoolProperty(name = "Vertical orientation", description = "Places the reconstructed rectangle in vertical orientation", default = False)
    size_property : bpy.props.FloatProperty(name="Size", description = "Size of the reconstructed rectangle", default = 1.0, min = 0.0, soft_min = 0.0, unit = "LENGTH")

    @classmethod
    def poll(cls, context):
        return context.active_object is not None and context.space_data.type == "PROPERTIES"

    def execute(self, context):
        # Get the camere of the scene
        scn = context.scene
        # Get the currently selected object
        obj = bpy.context.active_object
        # Check whether a mesh with 4 vertices in one polygon is selected
        if not obj.data.name in bpy.data.meshes or not len(obj.data.vertices) == 4 or not len(obj.data.polygons) == 1 or not len(obj.data.polygons[0].vertices) == 4:
            self.report({'ERROR'}, "Selected object must be a mesh with 4 vertices in 1 polygon.")
            return {'CANCELLED'}
        # Get the vertex coordinates and transform them to get the global coordinates, then project to 2d
        pa = transformation.vertex_apply_transformation(obj.data.vertices[obj.data.polygons[0].vertices[0]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        pb = transformation.vertex_apply_transformation(obj.data.vertices[obj.data.polygons[0].vertices[1]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        pc = transformation.vertex_apply_transformation(obj.data.vertices[obj.data.polygons[0].vertices[2]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        pd = transformation.vertex_apply_transformation(obj.data.vertices[obj.data.polygons[0].vertices[3]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        # Check whether the polygon is convex (this also checks for degnerate polygons)
        if not cameraplane.is_convex(pa, pb, pc, pd):
            self.report({'ERROR'}, "The polygon in the mesh must be convex and may not be degenerate.")
            return {'CANCELLED'}
        # Check for parallel edges
        if cameraplane.is_trapezoid(pa, pb, pc, pd):
            self.report({'ERROR'}, "Edges of the input rectangle must not be parallel.")
            return {'CANCELLED'}
        print("Vertices:", pa, pb, pc, pd)

        # Get the properties
        props = bpy.context.scene.camera_calibration_pvr_properties
        # Reference image
        image_obj = props.image
        if not image_obj:
            self.report({'ERROR'}, "Set a reference image.")
            return {'CANCELLED'}
        image, w, h, scale, offx, offy = reference.get_reference_image_data(image_obj)
        # Scale is the horizontal dimension. If in portrait mode, use the vertical dimension.
        if h > w:
            scale = scale / w * h

        # Perform the actual calibration
        cam_focal, cam_pos, cam_rot, coords, rec_size = calibrate_camera_F_PR_S(pa, pb, pc, pd, scale)

        if self.size_property > 0:
            size_factor = self.size_property / rec_size
        else:
            size_factor = 1.0 / rec_size
        cam_obj, cam = scene.get_or_create_camera(scn)
        # Set intrinsic camera parameters
        scene.set_camera_parameters(cam, lens = cam_focal)
        # Set background image
        reference.camera_apply_reference_image(cam, image)
        # Set extrinsic camera parameters and add a new rectangle
        scene.update_scene(cam_obj, cam_pos, cam_rot, self.vertical_property, scn, w, h, obj.name, coords, size_factor)

        # Switch to the active camera
        area = bpy.context.area.type
        bpy.context.area.type = "VIEW_3D"
        bpy.ops.view3d.view_camera()
        bpy.context.area.type = area

        return {'FINISHED'}

### Operator FX PR V #############################################################

class CameraCalibration_FX_PR_V_Operator(bpy.types.Operator):
    """Calibrates the focal length, vertical lens shift, position and rotation of the active camera."""
    bl_idname = "camera.camera_calibration_fx_pr_v"
    bl_label = "Solve Focal+Y"
    bl_options = {"REGISTER", "UNDO"}

    # Properties
    vertical_property : bpy.props.BoolProperty(name = "Vertical orientation", description = "Places the reconstructed rectangle in vertical orientation", default = False)
    size_property : bpy.props.FloatProperty(name="Size", description = "Size of the reconstructed rectangle", default = 1.0, min = 0.0, soft_min = 0.0, unit = "LENGTH")

    @classmethod
    def poll(cls, context):
        return context.active_object is not None and context.space_data.type == "PROPERTIES"

    def execute(self, context):
        # Get the camere of the scene
        scn = context.scene
        # Get the currently selected object
        obj = bpy.context.active_object
        # Check whether it is a mesh with 5 vertices, 4 in a polygon, 1 dangling at an edge
        if not obj.data.name in bpy.data.meshes or not len(obj.data.vertices) == 5 or not len(obj.data.polygons) == 1 or not len(obj.data.polygons[0].vertices) == 4 or not len(obj.data.edges) == 5:
            self.report({'ERROR'}, "Selected object must be a mesh with 4 vertices in 1 polygon and one dangling vertex.")
            return {'CANCELLED'}
        # Get the edge that is not part of the polygon
        dangling_edge = None
        for edge in obj.data.edges:
            if not edge.key in obj.data.polygons[0].edge_keys:
                dangling_edge = edge
                break
        print("Dangling edge:", dangling_edge.key)
        # Get the index to the attached and dangling vertex
        if dangling_edge.key[0] in obj.data.polygons[0].vertices:
            dangling_vertex = dangling_edge.key[1]
            attached_vertex = dangling_edge.key[0]
        else:
            dangling_vertex = dangling_edge.key[0]
            attached_vertex = dangling_edge.key[1]
        print("Dangling vertex:", dangling_vertex)
        print("Attached vertex:", attached_vertex)
        print("Polygon edges:", obj.data.polygons[0].edge_keys)
        # Get the vertex coordinates and apply the transformation to get global coordinates, then project to 2d
        pa = transformation.vertex_apply_transformation(obj.data.vertices[obj.data.polygons[0].vertices[0]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        pb = transformation.vertex_apply_transformation(obj.data.vertices[obj.data.polygons[0].vertices[1]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        pc = transformation.vertex_apply_transformation(obj.data.vertices[obj.data.polygons[0].vertices[2]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        pd = transformation.vertex_apply_transformation(obj.data.vertices[obj.data.polygons[0].vertices[3]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        pe = transformation.vertex_apply_transformation(obj.data.vertices[attached_vertex].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        pf = transformation.vertex_apply_transformation(obj.data.vertices[dangling_vertex].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        # Check whether the polygon is convex (this also checks for degnerate polygons)
        if not cameraplane.is_convex(pa, pb, pc, pd):
            self.report({'ERROR'}, "The polygon in the mesh must be convex and may not be degenerate.")
            return {'CANCELLED'}
        # Check for parallel edges
        if not cameraplane.is_trapezoid_but_not_rectangle(pa, pb, pc, pd):
            self.report({'ERROR'}, "Exactly one opposing edge pair of the input rectangle must be parallel.")
            return {'CANCELLED'}
        print("Vertices:", pa, pb, pc, pd, pe, pf)

        # Get the properties
        props = bpy.context.scene.camera_calibration_pvr_properties
        # Reference image
        image_obj = props.image
        if not image_obj:
            self.report({'ERROR'}, "Set a reference image.")
            return {'CANCELLED'}
        image, w, h, scale, offx, offy = reference.get_reference_image_data(image_obj)
        # Scale is the horizontal dimension. If in portrait mode, use the vertical dimension.
        if h > w:
            scale = scale / w * h
        # Perform the actual calibration
        calibration_data = calibrate_camera_FX_PR_V(pa, pb, pc, pd, pe, pf, scale)
        cam_focal, cam_pos, cam_rot, coords, rec_size, camera_shift = calibration_data
        if self.size_property > 0:
            size_factor = self.size_property / rec_size
        else:
            size_factor = 1.0 / rec_size
        cam_obj, cam = scene.get_or_create_camera(scn)
        # Set intrinsic camera parameters
        scene.set_camera_parameters(cam, lens = cam_focal, shift_y = camera_shift)
        # Set background image
        reference.camera_apply_reference_image(cam, image)
        # Set extrinsic camera parameters and add a new rectangle
        scene.update_scene(cam_obj, cam_pos, cam_rot, self.vertical_property, scn, w, h, obj.name, coords, size_factor)

        # Switch to the active camera
        area = bpy.context.area.type
        bpy.context.area.type = "VIEW_3D"
        bpy.ops.view3d.view_camera()
        bpy.context.area.type = area

        return {'FINISHED'}

### Operator FXY PR VV ############################################################

class CameraCalibration_FXY_PR_VV_Operator(bpy.types.Operator):
    """Calibrates the focal length, lens shift (horizontal and vertical), position and rotation of the active camera."""
    bl_idname = "camera.camera_calibration_fxy_pr_vv"
    bl_label = "Solve Focal+X+Y"
    bl_options = {"REGISTER", "UNDO"}

    # Properties
    vertical_property : bpy.props.BoolProperty(name = "Vertical orientation", description = "Places the reconstructed rectangle in vertical orientation", default = False)
    size_property : bpy.props.FloatProperty(name="Size", description = "Size of the reconstructed rectangle", default = 1.0, min = 0.0, soft_min = 0.0, unit = "LENGTH")

    @classmethod
    def poll(cls, context):
        return context.active_object is not None and context.space_data.type == "PROPERTIES"

    def execute(self, context):
        # Get the camere of the scene
        scn = context.scene
        # Get the currently selected object
        obj = bpy.context.active_object
        # Check whether it is a mesh with 6 vertices, 1 polygon, with 4 vertices and 2 dangling vertices
        if not obj.data.name in bpy.data.meshes or not len(obj.data.vertices) == 6 or not len(obj.data.polygons) == 1 or not len(obj.data.polygons[0].vertices) == 4 or not len(obj.data.edges) == 6:
            self.report({'ERROR'}, "Selected object must be a mesh with one polygon of 4 vertices with two dangling vertices.")
            return {'CANCELLED'}
        # Get the edges that are not part of the polygon
        print("Polygon edges:", obj.data.polygons[0].edge_keys)
        dangling_edges = []
        for edge in obj.data.edges:
            if not edge.key in obj.data.polygons[0].edge_keys:
                dangling_edges.append(edge)
        print("Dangling edges:", dangling_edges[0].key, dangling_edges[1].key)
        # Get the indices of the attached and dangling vertices
        dangling_vertices = [0, 0]
        attached_vertices = [0, 0]
        for i in range(2):
            if dangling_edges[i].key[0] in obj.data.polygons[0].vertices:
                dangling_vertices[i] = dangling_edges[i].key[1]
                attached_vertices[i] = dangling_edges[i].key[0]
            else:
                dangling_vertices[i] = dangling_edges[i].key[0]
                attached_vertices[i] = dangling_edges[i].key[1]
        print("Dangling vertices:", dangling_vertices)
        print("Attached vertices:", attached_vertices)
        # Convert indices to vertices
        for i in range(2):
            dangling_vertices[i] = transformation.vertex_apply_transformation(obj.data.vertices[dangling_vertices[i]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
            attached_vertices[i] = transformation.vertex_apply_transformation(obj.data.vertices[attached_vertices[i]].co, obj.scale, obj.rotation_euler, obj.location).to_2d()
        # Get the vertex coordinates and apply the transformation to get global coordinates, then project to 2d
        vertices = []
        for i in range(4):
            index = obj.data.polygons[0].vertices[i]
            vertices.append(transformation.vertex_apply_transformation(obj.data.vertices[index].co, obj.scale, obj.rotation_euler, obj.location).to_2d())
        # Check whether the polygon is convex (this also checks for degnerate polygons)
        if not cameraplane.is_convex(vertices[0], vertices[1], vertices[2], vertices[3]):
            self.report({'ERROR'}, "The polygon in the mesh must be convex and may not be degenerate.")
            return {'CANCELLED'}
        # Check for parallel edges
        if cameraplane.is_collinear(vertices[0] - vertices[1], vertices[3] - vertices[2]) or cameraplane.is_collinear(vertices[0] - vertices[3], vertices[1] - vertices[2]) or cameraplane.is_collinear(dangling_vertices[0] - attached_vertices[0], dangling_vertices[1] - attached_vertices[1]):
            self.report({'ERROR'}, "Edges must not be parallel.")
            return {'CANCELLED'}

        # Get the properties
        props = bpy.context.scene.camera_calibration_pvr_properties
        # Reference image
        image_obj = props.image
        if not image_obj:
            self.report({'ERROR'}, "Set a reference image.")
            return {'CANCELLED'}
        image, w, h, scale, offx, offy = reference.get_reference_image_data(image_obj)
        # Scale is the horizontal dimension. If in portrait mode, use the vertical dimension.
        if h > w:
            scale = scale / w * h

        # Perform the actual calibration
        calibration_data = calibrate_camera_FXY_PR_VV(vertices, attached_vertices, dangling_vertices, scale)
        cam_focal, cam_pos, cam_rot, coords, rec_size, camera_shift_x, camera_shift_y = calibration_data

        if self.size_property > 0:
            size_factor = self.size_property / rec_size
        else:
            size_factor = 1.0 / rec_size
        cam_obj, cam = scene.get_or_create_camera(scn)
        # Set intrinsic camera parameters
        scene.set_camera_parameters(cam, lens = cam_focal, shift_x = camera_shift_x, shift_y = camera_shift_y)
        # Set extrinsic camera parameters and add a new rectangle
        scene.update_scene(cam_obj, cam_pos, cam_rot, self.vertical_property, scn, w, h, obj.name, coords, size_factor)
        # Set background image
        reference.camera_apply_reference_image(cam, image)

        # Switch to the active camera
        area = bpy.context.area.type
        bpy.context.area.type = "VIEW_3D"
        bpy.ops.view3d.view_camera()
        bpy.context.area.type = area

        return {'FINISHED'}

### Operator FXY P S ###############################################################

class CameraCalibration_FXY_P_S_Operator(bpy.types.Operator):
    """Calibrates the focal length, lens shift (horizontal and vertical), position and rotation of the active camera."""
    bl_idname = "camera.camera_calibration_fxy_p_s"
    bl_label = "Solve 1-point"
    bl_options = {"REGISTER", "UNDO"}

    # Properties
    mode_property : bpy.props.EnumProperty(items=[("use_focal", "Focal Length Mode", "Uses a fixed focal length for the camera", 0), ("use_length", "Rectangle Length Mode", "Uses a fixed rectangle length", 1)], name="Mode")
    focal_property : bpy.props.FloatProperty(name="Focal Length", description = "Focal length of the camera (mm)", default = 35, min = 0.0, soft_min = 0.0, unit = "LENGTH")
    width_property : bpy.props.FloatProperty(name="Width", description = "Width of the reconstructed rectangle", default = 1.0, min = 0.0, soft_min = 0.0, unit = "LENGTH")
    length_property : bpy.props.FloatProperty(name="Length", description = "Length of the reconstructed rectangle", default = 1.0, min = 0.0, soft_min = 0.0, unit = "LENGTH")
    vertical_property : bpy.props.BoolProperty(name = "Vertical orientation", description = "Places the reconstructed rectangle in vertical orientation", default = False)

    @classmethod
    def poll(cls, context):
        return context.active_object is not None and context.space_data.type == "PROPERTIES"

    def execute(self, context):
        # Get the camere of the scene
        scn = context.scene
        # Get the currently selected object
        obj = bpy.context.active_object
        # Check whether it is a mesh with 4 vertices in 1 polygon
        if not obj.data.name in bpy.data.meshes or not len(obj.data.vertices) == 4 or not len(obj.data.polygons) == 1 or not len(obj.data.polygons[0].vertices) == 4:
            self.report({'ERROR'}, "Selected object must be a mesh with one polygon of 4 vertices with two horizontal edges.")
            return {'CANCELLED'}
        # Get the vertex coordinates and apply the transformation to get global coordinates, then project to 2d
        vertices = []
        for i in range(4):
            index = obj.data.polygons[0].vertices[i]
            vertices.append(transformation.vertex_apply_transformation(obj.data.vertices[index].co, obj.scale, obj.rotation_euler, obj.location).to_2d())
        # Check whether the polygon is convex (this also checks for degnerate polygons)
        if not cameraplane.is_convex(vertices[0], vertices[1], vertices[2], vertices[3]):
            self.report({'ERROR'}, "The polygon in the mesh must be convex and may not be degenerate.")
            return {'CANCELLED'}
        # Check for parallel edges
        if not cameraplane.is_trapezoid_but_not_rectangle(*vertices):
            self.report({'ERROR'}, "Exactly two opposing edges must be parallel.")
            return {'CANCELLED'}

        # Get the properties
        props = bpy.context.scene.camera_calibration_pvr_properties
        # Reference image
        image_obj = props.image
        if not image_obj:
            self.report({'ERROR'}, "Set a reference image.")
            return {'CANCELLED'}
        image, w, h, scale, offx, offy = reference.get_reference_image_data(image_obj)
        # Scale is the horizontal dimension. If in portrait mode, use the vertical dimension.
        if h > w:
            scale = scale / w * h

        # Prepare the calibration
        width = self.width_property
        length = self.length_property
        if width <= 0:
            width = 1
        if length <= 0:
            length = 1
        if self.mode_property == "use_focal":
            focal = self.focal_property
            if focal <= 0:
                focal = 35
            length = None
        else:
            focal = None
        # For testing only
        if False:
            verts = vertices.copy()
            for j in range(2):
                for i in range(4):
                    # Test the calibration with the given order of vertices
                    calibration_data = calibrate_camera_FXY_P_S(*verts, scale, focal, width, length)
                    print("Test:", *calibration_data)
                    # Roll the vertices by one
                    verts = verts[1:] + verts[:1]
                # Reverse the list of vertices
                verts.reverse()
        # Perform calibration
        calibration_data = calibrate_camera_FXY_P_S(*vertices, scale, focal, width, length)
        cam_focal, cam_pos, cam_rot, camera_shift_x, camera_shift_y, coords, length = calibration_data

        if self.mode_property == "use_focal":
            self.length_property = length
        else:
            self.focal_property = cam_focal

        cam_rot = mathutils.Euler((pi / 2, cam_rot, 0), "XYZ")

        cam_obj, cam = scene.get_or_create_camera(scn)
        # Set intrinsic camera parameters
        scene.set_camera_parameters(cam, lens = cam_focal, shift_x = camera_shift_x, shift_y = camera_shift_y)
        # Set background image
        reference.camera_apply_reference_image(cam, image)
        # Set extrinsic camera parameters and add a new rectangle
        scene.update_scene(cam_obj, cam_pos, cam_rot, self.vertical_property, scn, w, h, obj.name, coords, 1.0)

        # Switch to the active camera
        area = bpy.context.area.type
        bpy.context.area.type = "VIEW_3D"
        bpy.ops.view3d.view_camera()
        bpy.context.area.type = area

        return {'FINISHED'}

    def draw(self, context):
        layout = self.layout
        layout.props_enum(self, "mode_property")
        row = layout.row()
        row.prop(self, "focal_property")
        row.enabled = self.mode_property == "use_focal"
        layout.prop(self, "width_property")
        row = layout.row()
        row.prop(self, "length_property")
        row.enabled = self.mode_property == "use_length"
        layout.prop(self, "vertical_property")
