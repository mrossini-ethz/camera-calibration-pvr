import bpy

class CameraCalibrationPVROperator(bpy.types.Operator):
    bl_idname = "camera.camera_calibration_pvr"
    bl_label = "Solve Camera Calibration"

    def execute(self, context):
        print("Hello World")
        return {'FINISHED'}
