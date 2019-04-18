import bpy

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
