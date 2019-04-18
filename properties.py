import bpy

class CameraCalibrationPVRProperties(bpy.types.PropertyGroup):
    image : bpy.props.PointerProperty(type = bpy.types.Object,
                                      name = "Image",
                                      description = "An 'Empty' object of type 'Image'. The camera will be calibrated using this reference image")
    perspective: bpy.props.EnumProperty(items = [("3_point", "3 point", "Reference image is a 3 point perspective", 0),
                                                 ("2_point", "2 point", "Reference image is a 2 point perspective", 1),
                                                 ("1_point", "1 point", "Reference image is a 1 point perspective", 2)],
                                        name="Mode",
                                        default="3_point")
    focal: bpy.props.FloatProperty(name="Focal length", default=35)
    x_shift: bpy.props.FloatProperty(name="X-shift", default=0)
    y_shift: bpy.props.FloatProperty(name="Y-shift", default=0)
    width: bpy.props.FloatProperty(name="Rectangle width", default=1)
    length: bpy.props.FloatProperty(name="Rectangle length", default=1)
    focal_auto: bpy.props.BoolProperty(name="Solve", default=True)
    x_shift_auto: bpy.props.BoolProperty(name="Solve", default=True)
    y_shift_auto: bpy.props.BoolProperty(name="Solve", default=True)
    width_auto: bpy.props.BoolProperty(name="Solve", default=True)
    length_auto: bpy.props.BoolProperty(name="Solve", default=True)
    rect_orientation: bpy.props.EnumProperty(items = [("xy", "XY Plane", "Rectangle should be parallel to the XY plane", 0),
                                                      ("xz", "XZ Plane", "Rectangle should be parallel to the XZ plane", 1),
                                                      ("yz", "YZ Plane", "Rectangle should be parallel to the YZ plane", 2)],
                                             name="Rectangle orientation",
                                             default="xy")
    rect_orientation_auto: bpy.props.BoolProperty(name="Auto", default=True)
