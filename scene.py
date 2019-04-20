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

def object_name_append(name, suffix):
    # Check whether the object name is numbered
    if len(name) > 4 and name[-4] == "." and name[-3:].isdecimal():
        return name[:-4] + suffix + name[-4:]
    return name + suffix

def get_or_create_camera(scn):
    cam_obj = scn.camera
    if not cam_obj:
        bpy.ops.object.camera_add()
        cam_obj = bpy.context.active_object
    cam = bpy.data.cameras[cam_obj.data.name]
    return (cam_obj, cam)

def set_camera_parameters(camera, lens = 35.0, shift_x = 0.0, shift_y = 0.0, sensor_size = 32.0):
    """Sets the intrinsic camera parameters, including default ones (to get reliable results)"""
    camera.lens = lens
    camera.lens_unit = "MILLIMETERS"
    camera.shift_x = shift_x
    camera.shift_y = shift_y
    camera.sensor_width = sensor_size
    camera.sensor_fit = "AUTO"
    camera.type = "PERSP"

def set_camera_transformation(camera_obj, translation, rotation):
    """Transforms the camera object according to the given parameters"""
    camera_obj.location = translation
    camera_obj.rotation_euler = rotation

def get_vertical_mode_matrix(is_vertical, camera_rotation):
    if is_vertical:
    # Get the up direction of the camera
        up_vec = mathutils.Vector((0.0, 1.0, 0.0))
        up_vec.rotate(camera_rotation)
        # Decide around which axis to rotate
        vert_mode_rotate_x = abs(up_vec[0]) < abs(up_vec[1])
        # Create rotation matrix
        if vert_mode_rotate_x:
            vert_angle = pi / 2 if up_vec[1] > 0 else -pi / 2
            return mathutils.Matrix().Rotation(vert_angle, 3, "X")
        else:
            vert_angle = pi / 2 if up_vec[0] < 0 else -pi / 2
            return mathutils.Matrix().Rotation(vert_angle, 3, "Y")
    else:
        return mathutils.Matrix().Identity(3)

def update_scene(camera, cam_pos, cam_rot, is_vertical, scn, img_width, img_height, object_name, coords, size_factor):
    """Updates the scene by moving the camera and creating a new rectangle"""
    # Get the 3D cursor location
    cursor_pos = mathutils.Vector((bpy.context.scene.cursor.location))
    # Get transformation matrix for vertical orientation
    vert_matrix = get_vertical_mode_matrix(is_vertical, cam_rot)
    # Set the camera position and rotation
    cam_rot = cam_rot.copy()
    cam_rot.rotate(vert_matrix)
    set_camera_transformation(camera, vert_matrix @ cam_pos * size_factor + cursor_pos, cam_rot)
    # Apply the transformation matrix for vertical orientation
    for i in range(4):
        coords[i].rotate(vert_matrix)
    # Set the render resolution
    scn.render.resolution_x = img_width
    scn.render.resolution_y = img_height
    # Add the rectangle to the scene (at the 3D cursor location)
    bpy.ops.mesh.primitive_plane_add()
    rect = bpy.context.active_object
    # Rename the rectangle
    rect.name = object_name_append(object_name, "_Cal")
    # Set the correct size (local coordinates)
    for i in range(4):
        rect.data.vertices[rect.data.polygons[0].vertices[i]].co = coords[i] * size_factor
