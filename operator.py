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

from . import reference

def get_or_create_camera():
    cam_obj = bpy.context.scene.camera
    if not cam_obj:
        bpy.ops.object.camera_add()
        cam_obj = bpy.context.object
    cam = bpy.data.cameras[cam_obj.data.name]
    return (cam_obj, cam)

class CameraCalibrationPVROperator(bpy.types.Operator):
    bl_idname = "camera.camera_calibration_pvr"
    bl_label = "Solve Camera Calibration"

    def execute(self, context):
        # Get the properties
        props = bpy.context.scene.camera_calibration_pvr_properties

        # Reference image
        image_obj = props.image
        if not image_obj:
            self.report({'ERROR'}, "Set a reference image.")
            return {'CANCELLED'}
        image, width, height, size, offset_x, offset_y = reference.get_reference_image_data(image_obj)

        # Get the active camera (or create one)
        camera_obj, camera = get_or_create_camera()

        # Set the background image for the active camera (if necessary)
        reference.camera_apply_reference_image(camera, image)

        # Set the render size orientation (landscape/portrait) according to the reference image
        # Do NOT change the render size
        rwidth = context.scene.render.resolution_x
        rheight = context.scene.render.resolution_y
        if width > height and rwidth < rheight or width < height and rwidth > rheight:
            # Switch orientation
            context.scene.render.resolution_x, context.scene.render.resolution_y = rheight, rwidth

        return {'FINISHED'}

### LEGACY Operators #############################################################

from math import sqrt, pi, degrees

from . import polynomial
from . import rootfinder
from . import algebra
from . import cameraplane
from . import transformation
from . import reference
from . import solverectangle
from . import scene
from . import onepoint
from . import twopoint
from . import threepoint

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
        cam_focal, cam_pos, cam_rot, coords, rec_size = threepoint.calibrate_camera(pa, pb, pc, pd, scale)

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

        # Switch to the active camera (if there is a 3D view)
        for area in bpy.context.screen.areas:
            if area.type == "VIEW_3D" and len(area.spaces) > 0:
                area.spaces[0].region_3d.view_perspective = "CAMERA"

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
        calibration_data = twopoint.calibrate_camera(pa, pb, pc, pd, pe, pf, scale)
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

        # Switch to the active camera (if there is a 3D view)
        for area in bpy.context.screen.areas:
            if area.type == "VIEW_3D" and len(area.spaces) > 0:
                area.spaces[0].region_3d.view_perspective = "CAMERA"

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
        calibration_data = threepoint.calibrate_camera_shifted(vertices, attached_vertices, dangling_vertices, scale)
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

        # Switch to the active camera (if there is a 3D view)
        for area in bpy.context.screen.areas:
            if area.type == "VIEW_3D" and len(area.spaces) > 0:
                area.spaces[0].region_3d.view_perspective = "CAMERA"

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
        calibration_data = onepoint.calibrate_camera(*vertices, scale, focal, width, length)
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

        # Switch to the active camera (if there is a 3D view)
        for area in bpy.context.screen.areas:
            if area.type == "VIEW_3D" and len(area.spaces) > 0:
                area.spaces[0].region_3d.view_perspective = "CAMERA"

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
