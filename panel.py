import bpy

from . import images

class PreviewsExamplePanel(bpy.types.Panel):
    """Creates a Panel in the Object properties window"""
    bl_label = "Camera Calibration PVR"
    bl_idname = "OBJECT_PT_camera_calibration_pvr"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "scene"

    def draw(self, context):
        obj = context.object
        props = bpy.context.scene.camera_calibration_pvr_properties
        icon_1_point = images.preview_collection["1_point"]
        icon_2_point = images.preview_collection["2_point"]
        icon_3_point = images.preview_collection["3_point"]

        layout = self.layout
        layout.prop(props, "image")

        layout.prop_tabs_enum(props, "perspective")

        box = layout.box()
        if props.perspective == "3_point":
            box.template_icon(icon_value=icon_3_point.icon_id, scale=5)
        if props.perspective == "2_point":
            box.template_icon(icon_value=icon_2_point.icon_id, scale=5)
        if props.perspective == "1_point":
            box.template_icon(icon_value=icon_1_point.icon_id, scale=5)

        split_factor = 0.8
        layout.separator()
        row = layout.split(factor = split_factor)
        sub = row.row()
        sub.enabled = not props.focal_auto
        sub.prop(props, "focal", expand = True)
        sub = row.row()
        sub.enabled = not(props.perspective == "1_point" and props.length_auto == True and props.length_auto != props.focal_auto)
        sub.prop(props, "focal_auto", expand = False, icon_only = False, toggle = True)

        row = layout.split(factor = split_factor)
        sub = row.row()
        sub.enabled = not props.x_shift_auto or props.perspective == "2_point"
        sub.prop(props, "x_shift", expand = True)
        sub = row.row()
        sub.enabled = props.perspective != "2_point"
        sub.prop(props, "x_shift_auto", expand = False, icon_only = False, toggle = True)

        row = layout.split(factor = split_factor)
        sub = row.row()
        sub.enabled = not props.y_shift_auto
        sub.prop(props, "y_shift", expand = True)
        sub = row.row()
        sub.enabled = True
        sub.prop(props, "y_shift_auto", expand = False, icon_only = False, toggle = True)

        row = layout.split(factor = split_factor)
        sub = row.row()
        sub.enabled = not props.width_auto
        sub.prop(props, "width", expand = True)
        sub = row.row()
        sub.enabled = True
        sub.prop(props, "width_auto", expand = False, icon_only = False, toggle = True)

        row = layout.split(factor = split_factor)
        sub = row.row()
        sub.enabled = not props.length_auto or not (props.perspective != "1_point" or props.focal_auto != True)
        sub.prop(props, "length", expand = True)
        sub = row.row()
        sub.enabled = props.perspective != "1_point" or props.focal_auto != True
        sub.prop(props, "length_auto", expand = False, icon_only = False, toggle = True)

        layout.separator()
        row = layout.split(factor = split_factor)
        sub = row.row()
        sub.enabled = not props.rect_orientation_auto
        sub.prop_tabs_enum(props, "rect_orientation")
        sub = row.row()
        sub.enabled = True
        sub.prop(props, "rect_orientation_auto", expand = False, icon_only = False, toggle = True)

        layout.separator()
        layout.operator("camera.camera_calibration_pvr")
