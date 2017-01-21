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
from math import sqrt, pi, atan2

bl_info = {
    "name": "Camera Calibration using Perspective Views of Rectangles",
    "author": "Marco Rossini",
    "version": (0, 0, 1),
    "blender": (2, 7, 0),
    "location": "3D View > Tools Panel > Misc > Camera Calibration",
    "description": "Calibrates position, rotation and focal length of a camera using a single image of a rectangle.",
    "tracker_url": "https://github.com/mrossini-ethz/camera-calibration-pvr/issues",
    "support": "COMMUNITY",
    "category": "3D View"
}

### Polynomials ##################################################################

def make_poly(coeffs):
    """Make a new polynomial"""
    return list(coeffs)

def poly_norm(poly):
    """Normalizes a given polynomial"""
    f = poly[-1]
    result = []
    for coeff in poly:
        result.append(coeff / f)
    return result

def poly_sub(a, b):
    """Subtract the two polynomials"""
    n = max(len(a), len(b))
    _a = [0] * n
    _b = [0] * n
    for i in range(len(a)):
        _a[i] = a[i]
    for i in range(len(b)):
        _b[i] = b[i]
    result = []
    for i in range(n):
        result.append(_a[i] - _b[i])
    return result

def poly_scale(poly, factor):
    """Normalizes a given polynomial"""
    f = poly[-1]
    result = []
    for coeff in poly:
        result.append(coeff * factor)
    return result

def poly_reduce(poly):
    """Removes leading coefficients that are zero"""
    result = []
    for i in range(len(poly) - 1, -1, -1):
        if poly[i] != 0 or len(result) > 0:
            result.append(poly[i])
    result.reverse()
    return result

def poly_derivative(poly):
    """Calculates the derivative of the polynomial"""
    result = []
    for i in range(1, len(poly)):
        result.append(i * poly[i])
    return result

def poly_eval(poly, x):
    """Evaluate the polynomial"""
    result = 0.0
    for i in range(len(poly)):
        result += poly[i] * x ** i
    return result

def poly_order(poly):
    """Get the order of the polynomial"""
    return len(poly) - 1

def poly_coeff(poly, idx):
    """Get the nth coefficient of the polynomial"""
    if idx > len(poly) - 1:
        return 0.0
    elif idx >= 0:
        return poly[idx]

def poly_div(a, b):
    """Calculate the polynom division of a and b"""
    na = poly_order(a)
    nb = poly_order(b)
    result = [0] * (na - nb + 1)
    for n in range(na, nb - 1, -1):
        f = a[n] / b[-1]
        result[n - nb] = f
        a = poly_sub(a, [0] * (n - nb) + poly_scale(b, f))
    return result

### Root Finder ##################################################################

def find_root(f, df, ddf, initial_guess = 0.0, limit = 0.00001, max_iterations = 1000):
    """Find the root of the function f using Halley's method"""
    xn_1 = initial_guess
    i = 0
    while i < max_iterations:
        fx = f(xn_1)
        dfx = df(xn_1)
        ddfx = ddf(xn_1)
        xn = xn_1 - 2 * fx * dfx / (2 * dfx ** 2 - fx * ddfx)
        if abs(xn - xn_1) < limit:
            return xn
        xn_1 = xn
        i += 1
    return None

def find_poly_root(poly, initial_guess = 0.0, limit = 0.00001, max_iterations = 1000):
    """Find a root of the given polynomial"""
    # Calculate the polynomial derivatives
    dpoly = poly_derivative(poly)
    ddpoly = poly_derivative(dpoly)
    # Closures !!!
    f = lambda x: poly_eval(poly, x)
    df = lambda x: poly_eval(dpoly, x)
    ddf = lambda x: poly_eval(ddpoly, x)
    # Call the generic root finder
    return find_root(f, df, ddf, initial_guess, limit, max_iterations)

def find_poly_roots(poly, initial_guess = 0.0, limit = 0.00001, max_iterations = 1000):
    """Find all roots of the given polynomial"""
    solutions = []
    # Find solutions numerically for n > 0, split them off until n = 2
    for q in range(poly_order(poly) - 2):
        x = find_poly_root(poly, initial_guess, limit, max_iterations)
        if not x:
            break
        poly = poly_div(poly, make_poly([-x, 1]))
        solutions.append(x)
    # Find the rest of the roots analytically
    if poly_order(poly) == 1:
        solutions.append(- poly_coeff(poly, 1) / poly_coeff(poly, 0))
    elif poly_order(poly) == 2:
        a = poly_coeff(poly, 2)
        b = poly_coeff(poly, 1)
        c = poly_coeff(poly, 0)
        d = b ** 2 - 4 * a * c
        if d == 0:
            solutions.append(-b / (2 * a))
        elif d > 0:
            solutions.append((- b + sqrt(d)) / (2 * a))
            solutions.append((- b - sqrt(d)) / (2 * a))
    return solutions

### Algorithm ####################################################################

def intersect_2d(pa, pb, pc, pd):
    """Find the intersection point of the lines AB and CD (2 dimensions)"""
    # FIXME: proper solution of linear system to avoid division by zero
    d1 = pb - pa
    d2 = pd - pc
    r = pc - pa
    s = (r[0] - r[1] * d1[0] / d1[1]) / (d1[0] *  d2[1] / d1[1] - d2[0])
    return pc + s * d2

def get_vanishing_points(pa, pb, pc, pd):
    """Get the two vanishing points of the rectangle defined by the corners pb pb pc pd"""
    return (intersect_2d(pa, pb, pd, pc), intersect_2d(pa, pd, pb, pc))

def get_camera_plane_vector(p, scale, focal_length = 1.0):
    """Convert a 2d point in the camera plane into a 3d vector from the camera onto the camera plane"""
    # field_of_view = 2 * atan(sensor_size / 2 * focal_length), assume sensor_size = 32
    s = (16.0 / focal_length) / (scale / 2.0)
    return mathutils.Vector((p[0] * s, p[1] * s, -1.0))

def calculate_focal_length(pa, pb, pc, pd, scale):
    """Get the vanishing points of the rectangle as defined by pa, pb, pc and pd"""
    pm, pn = get_vanishing_points(pa, pb, pc, pd)
    # Calculate the vectors from camera to the camera plane where the vanishing points are located
    vm = get_camera_plane_vector(pm, scale)
    vn = get_camera_plane_vector(pn, scale)
    # Calculate the focal length
    return sqrt(abs(vm.dot(vn)))

def get_lambda_d_poly_a(qab, qac, qad, qbc, qbd, qcd):
    """Equation A (see paper)"""
    d4 = qac * qbd ** 2 - qad * qbc * qbd
    d3 = qab * qad * qbc + qad ** 2 * qbc * qbd + qad ** 2 * qcd + qbc * qbd - 2 * qab * qac * qbd - qab * qad * qbd * qcd - qac * qad * qbd ** 2
    d2 = qab ** 2 * qac + qab ** 2 * qad * qcd + 3 * qab * qac * qad * qbd + qab * qbd * qcd - qab * qad ** 2 * qbc - qab * qbc - qac * qad ** 2 - qad * qbc * qbd - 2 * qad * qcd
    d1 = qab * qad * qbc + 2 * qac * qad + qcd - 2 * qab ** 2 * qac * qad - qab ** 2 * qcd - qab * qac * qbd
    d0 = qab ** 2 * qac - qac
    return make_poly([d0, d1, d2, d3, d4])

def get_lambda_d_poly_b(qab, qac, qad, qbc, qbd, qcd):
    """Equation B (see paper)"""
    d4 = qbd - qbd * qcd ** 2
    d3 = qab * qcd ** 2 + qac * qbd * qcd + 2 * qad * qbd * qcd ** 2 - qab - 2 * qad * qbd - qad * qbc * qcd
    d2 = 2 * qab * qad + qac * qad * qbc + qad ** 2 * qbc * qcd + qad **2 * qbd + qbc * qcd - qab * qac * qcd - qab * qad * qcd ** 2 - 3 * qac * qad * qbd * qcd - qbd * qcd ** 2
    d1 = qab * qac * qad * qcd + qac ** 2 * qad * qbd + 2 * qac * qbd * qcd - qab * qad ** 2 - qac * qad ** 2 * qbc - qac * qbc - qad * qbc * qcd
    d0 = qac * qad * qbc - qac ** 2 * qbd
    return make_poly([d0, d1, d2, d3, d4])

def get_lambda_d(pa, pb, pc, pd, scale, focal_length):
    """Calculate the vectors from camera to the camera plane where the rectangle corners are located"""
    va = get_camera_plane_vector(pa, scale, focal_length).normalized()
    vb = get_camera_plane_vector(pb, scale, focal_length).normalized()
    vc = get_camera_plane_vector(pc, scale, focal_length).normalized()
    vd = get_camera_plane_vector(pd, scale, focal_length).normalized()
    # Calculate dot products
    qab = va.dot(vb)
    qac = va.dot(vc)
    qad = va.dot(vd)
    qbc = vb.dot(vc)
    qbd = vb.dot(vd)
    qcd = vc.dot(vd)
    # Determine the equation that needs to be solved
    pa = poly_norm(get_lambda_d_poly_a(qab, qac, qad, qbc, qbd, qcd))
    pb = poly_norm(get_lambda_d_poly_b(qab, qac, qad, qbc, qbd, qcd))
    print("A:", pa)
    print("B:", pb)
    p = poly_reduce(poly_sub(pa, pb))
    print("P:", p)
    # Solve the equation
    roots = find_poly_roots(p)
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
        angles = [(rb - ra).angle(rd - ra) * 180 / pi, (ra - rb).angle(rc - rb) * 180 / pi, (rb - rc).angle(rd - rc) * 180 / pi, (rc - rd).angle(ra - rd) * 180 / pi]
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

def calibrate_camera_from_rectangle(pa, pb, pc, pd, scale):
    # Calculate the focal length of the camera
    focal = calculate_focal_length(pa, pb, pc, pd, scale)
    # Calculate the coordinates of the rectangle in 3d
    coords = get_lambda_d(pa, pb, pc, pd, scale, focal)
    # Calculate the transformation of the rectangle
    trafo = get_transformation(coords[0], coords[1], coords[2], coords[3])
    # Reconstruct the rotation angles of the transformation
    angles = get_rot_angles(trafo[0], trafo[1], trafo[2])
    xyz_matrix = mathutils.Euler((angles[0], angles[1], angles[2]), "XYZ")
    # Reconstruct the camera position
    cam_pos = -trafo[-1]
    cam_pos.rotate(xyz_matrix)
    # Calculate the corners of the rectangle in 3d such that it lies on the xy-plane
    tr = trafo[-1]
    ca = coords[0] - tr
    cb = coords[1] - tr
    cc = coords[2] - tr
    cd = coords[3] - tr
    ca.rotate(xyz_matrix)
    cb.rotate(xyz_matrix)
    cc.rotate(xyz_matrix)
    cd.rotate(xyz_matrix)
    # Printout for debugging
    print("Focal length:", focal)
    print("Camera Rx:", angles[0] * 180 / pi)
    print("Camera Ry:", angles[1] * 180 / pi)
    print("Camera Rz:", angles[2] * 180 / pi)
    print("Camera x:", cam_pos[0])
    print("Camera y:", cam_pos[1])
    print("Camera z:", cam_pos[2])
    length = (coords[0] - coords[1]).length
    width = (coords[0] - coords[3]).length
    size = max(length, width)
    print("Rectangle length:", length)
    print("Rectangle width:", width)
    print("Rectangle A:", ca)
    print("Rectangle B:", cb)
    print("Rectangle C:", cc)
    print("Rectangle D:", cd)
    return (focal, cam_pos, xyz_matrix, [ca, cb, cc, cd], size)

### Operator #####################################################################

class CameraCalibrationOperator(bpy.types.Operator):
    """Calibrates the active camera"""
    bl_idname = "camera.camera_calibration"
    bl_label = "Camera Calibration"
    bl_options = {"REGISTER", "UNDO"}

    # Properties
    vertical_property = bpy.props.BoolProperty(name = "Vertical orientation", description = "Places the reconstructed rectangle in vertical orientation", default = False)
    size_property = bpy.props.FloatProperty(name="Size", description = "Size of the reconstructed rectangle", default = 1.0, min = 0.0, soft_min = 0.0, unit = "LENGTH")

    @classmethod
    def poll(cls, context):
        return context.active_object is not None and context.space_data.type == "VIEW_3D"

    def execute(self, context):
        # Get the camere of the scene
        scene = bpy.data.scenes["Scene"]
        cam_obj = scene.camera
        if not cam_obj:
            self.report({'ERROR'}, "There is no active camera.")
            return {'CANCELLED'}
        cam = bpy.data.cameras[cam_obj.data.name]
        # Get the currently selected object
        obj = bpy.context.object
        # Check whether a mesh with 4 vertices in one polygon is selected
        if not obj.data.name in bpy.data.meshes or not len(obj.data.vertices) == 4 or not len(obj.data.polygons) == 1 or not len(obj.data.polygons[0].vertices) == 4:
            self.report({'ERROR'}, "Selected object must be a mesh with 4 vertices in 1 polygon.")
            return {'CANCELLED'}
        # Get the vertex coordinates
        pa = obj.data.vertices[obj.data.polygons[0].vertices[0]].co.copy()
        pb = obj.data.vertices[obj.data.polygons[0].vertices[1]].co.copy()
        pc = obj.data.vertices[obj.data.polygons[0].vertices[2]].co.copy()
        pd = obj.data.vertices[obj.data.polygons[0].vertices[3]].co.copy()
        # Apply the scale
        obj_scale = obj.scale
        for i in range(3):
            pa[i] *= obj_scale[i]
            pb[i] *= obj_scale[i]
            pc[i] *= obj_scale[i]
            pd[i] *= obj_scale[i]
        # Apply rotation
        obj_rotation = obj.rotation_euler
        pa.rotate(obj_rotation)
        pb.rotate(obj_rotation)
        pc.rotate(obj_rotation)
        pd.rotate(obj_rotation)
        # Apply translation and project to x-y-plane
        obj_translation = obj.location
        pa = (pa + obj_translation).to_2d()
        pb = (pb + obj_translation).to_2d()
        pc = (pc + obj_translation).to_2d()
        pd = (pd + obj_translation).to_2d()
        print("Vertices:", pa, pb, pc, pd)
        # Get the background images
        bkg_images = bpy.context.space_data.background_images
        if len(bkg_images) == 1:
            # If there is only one background image, take that one
            img = bkg_images[0]
        else:
            # Get the visible background images with view axis 'top'
            bkg_images_top = []
            for img in bkg_images:
                if (img.view_axis == "TOP" or img.view_axis == "ALL") and img.show_background_image:
                    bkg_images_top.append(img)
            # Check the number of images
            if len(bkg_images_top) != 1:
                # Check only the TOP images
                bkg_images_top = []
                for img in bkg_images:
                    if img.view_axis == "TOP" and img.show_background_image:
                        bkg_images_top.append(img)
                if len(bkg_images_top) != 1:
                    self.report({'ERROR'}, "Exactly 1 visible background image required in top view.")
                    return {'CANCELLED'}
            # Get the background image properties
            img = bkg_images_top[0]
        offx = img.offset_x
        offy = img.offset_y
        rot = img.rotation
        scale = img.size
        flipx = img.use_flip_x
        flipy = img.use_flip_y
        w, h = img.image.size
        # Scale is the horizontal dimension. If in portrait mode, use the vertical dimension.
        if h > w:
            scale = scale / w * h
        # Perform the actual calibration
        cam_focal, cam_pos, cam_rot, coords, rec_size = calibrate_camera_from_rectangle(pa, pb, pc, pd, scale)
        if self.size_property > 0:
            size_factor = self.size_property / rec_size
        else:
            size_factor = 1.0 / rec_size
        cam.lens = cam_focal
        cam_obj.location = cam_pos * size_factor
        cam_obj.rotation_euler = cam_rot
        # Perform rotation to obtain vertical orientation, if necessary.
        if self.vertical_property:
            # Get the up direction of the camera
            up_vec = mathutils.Vector((0.0, 1.0, 0.0))
            up_vec.rotate(cam_rot)
            # Decide around which axis to rotate
            vert_mode_rotate_x = abs(up_vec[0]) < abs(up_vec[1])
            # Create rotation matrix
            if vert_mode_rotate_x:
                vert_angle = pi / 2 if up_vec[1] > 0 else -pi / 2
                vert_matrix = mathutils.Matrix().Rotation(vert_angle, 3, "X")
            else:
                vert_angle = pi / 2 if up_vec[0] < 0 else -pi / 2
                vert_matrix = mathutils.Matrix().Rotation(vert_angle, 3, "Y")
            # Apply matrix
            cam_obj.location.rotate(vert_matrix)
            cam_obj.rotation_euler.rotate(vert_matrix)
            for i in range(4):
                coords[i].rotate(vert_matrix)
        # Set the render resolution
        scene.render.resolution_x = w
        scene.render.resolution_y = h
        # Add the rectangle to the scene
        bpy.ops.mesh.primitive_plane_add()
        rect = bpy.context.object
        rect.name = obj.name
        for i in range(4):
            rect.data.vertices[rect.data.polygons[0].vertices[i]].co = coords[i] * size_factor
        # Switch to the active camera
        if not bpy.context.space_data.region_3d.view_perspective == "CAMERA":
            bpy.ops.view3d.viewnumpad(type="CAMERA")
        return {'FINISHED'}

### Panel ########################################################################

class CameraCalibrationPanel(bpy.types.Panel):
    """Creates a Panel in the scene context of the properties editor"""
    bl_label = "Camera Calibration"
    bl_idname = "VIEW_3D_camera_calibration"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'

    def draw(self, context):
        layout = self.layout
        layout.operator("camera.camera_calibration")

### Register #####################################################################

def register():
    bpy.utils.register_class(CameraCalibrationPanel)
    bpy.utils.register_class(CameraCalibrationOperator)

def unregister():
    bpy.utils.unregister_class(CameraCalibrationPanel)
    bpy.utils.unregister_class(CameraCalibrationOperator)

if __name__ == "__main__":
    register()
