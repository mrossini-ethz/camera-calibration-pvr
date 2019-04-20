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
    imp.reload(polynomial)
    imp.reload(rootfinder)
    imp.reload(algebra)
    imp.reload(cameraplane)
    imp.reload(transformation)
    imp.reload(main)
    imp.reload(panel)
else:
    from . import images
    from . import reference
    from . import operator
    from . import properties
    from . import polynomial
    from . import rootfinder
    from . import algebra
    from . import cameraplane
    from . import transformation
    from . import main
    from . import panel

import bpy

def register():
    images.load_images()

    bpy.utils.register_class(properties.CameraCalibrationPVRProperties)
    bpy.types.Scene.camera_calibration_pvr_properties = bpy.props.PointerProperty(type=properties.CameraCalibrationPVRProperties)
    bpy.utils.register_class(operator.CameraCalibrationPVROperator)
    bpy.utils.register_class(main.CameraCalibration_F_PR_S_Operator)
    bpy.utils.register_class(main.CameraCalibration_FX_PR_V_Operator)
    bpy.utils.register_class(main.CameraCalibration_FXY_PR_VV_Operator)
    bpy.utils.register_class(main.CameraCalibration_FXY_P_S_Operator)
    bpy.utils.register_class(panel.PreviewsExamplePanel)


def unregister():
    bpy.utils.previews.remove(images.preview_collection)

    bpy.utils.unregister_class(panel.PreviewsExamplePanel)
    bpy.utils.unregister_class(main.CameraCalibration_F_PR_S_Operator)
    bpy.utils.unregister_class(main.CameraCalibration_FX_PR_V_Operator)
    bpy.utils.unregister_class(main.CameraCalibration_FXY_PR_VV_Operator)
    bpy.utils.unregister_class(main.CameraCalibration_FXY_P_S_Operator)
    bpy.utils.unregister_class(operator.CameraCalibrationPVROperator)
    bpy.utils.unregister_class(properties.CameraCalibrationPVRProperties)

if __name__ == "__main__":
    register()
