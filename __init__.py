bl_info = {
    "name": "Camera Calibration using Perspective Views of Rectangles",
    "author": "Marco Rossini",
    "version": (0, 5, 0),
    # "warning": "This is an unreleased development version.",
    "blender": (2, 80, 0),
    "location": "3D View > Properties > Scene",
    "description": "Calibrates position, rotation and focal length of a camera using a single image of a rectangle.",
    "wiki_url": "https://github.com/mrossini-ethz/camera-calibration-pvr",
    "tracker_url": "https://github.com/mrossini-ethz/camera-calibration-pvr/issues",
    "support": "COMMUNITY",
    "category": "Camera"
}

if "bpy" in locals():
    import importlib as imp
    imp.reload(images)
    imp.reload(reference)
    imp.reload(operator)
    imp.reload(properties)
    imp.reload(panel)
else:
    from . import images
    from . import reference
    from . import operator
    from . import properties
    from . import panel

import bpy

def register():
    images.load_images()

    bpy.utils.register_class(properties.CameraCalibrationPVRProperties)
    bpy.types.Scene.camera_calibration_pvr_properties = bpy.props.PointerProperty(type=properties.CameraCalibrationPVRProperties)
    bpy.utils.register_class(operator.CameraCalibrationPVROperator)
    bpy.utils.register_class(panel.PreviewsExamplePanel)


def unregister():
    bpy.utils.previews.remove(images.preview_collection)

    bpy.utils.unregister_class(panel.PreviewsExamplePanel)
    bpy.utils.unregister_class(operator.CameraCalibrationPVROperator)
    bpy.utils.unregister_class(properties.CameraCalibrationPVRProperties)

if __name__ == "__main__":
    register()
