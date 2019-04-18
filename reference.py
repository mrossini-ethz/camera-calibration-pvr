import bpy

def get_empty_image_data(obj):
    # Return image, width, height, scale, x-offset, y-offset for the given emtpy object (of image type)
    return (obj.data, obj.data.size[0], obj.data.size[1], obj.empty_display_size, obj.empty_image_offset[0], obj.empty_image_offset[1])
