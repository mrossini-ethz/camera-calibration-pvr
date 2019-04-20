import bpy

def get_reference_image_data(obj):
    # Return image, width, height, scale, x-offset, y-offset for the given emtpy object (of image type)
    return (obj.data, obj.data.size[0], obj.data.size[1], obj.empty_display_size, obj.empty_image_offset[0], obj.empty_image_offset[1])

def camera_apply_reference_image(camera, image):
    # Get the background images from the active camera
    needs_setting = True
    images = camera.background_images
    bkg_img = None

    # Check whether it needs setting
    for img in images:
        if img.image == image:
            bkg_img = img
            break

    if not bkg_img:
        # Add a background image to the camera
        area = bpy.context.area.type
        bpy.context.area.type = "VIEW_3D"
        bpy.ops.view3d.background_image_add()
        bpy.context.area.type = area

        # Set the background to the reference image
        images = camera.background_images
        bkg_img = images[-1]
        bkg_img.image = image

    # Ensure correct settings
    if bkg_img.alpha == 0.0:
        bkg_img.alpha = 0.5
    bkg_img.display_depth = "FRONT"
    bkg_img.frame_method = "CROP"
    bkg_img.offset[0] = 0.0
    bkg_img.offset[1] = 0.0
    bkg_img.rotation = 0.0
    bkg_img.scale = 1.0
    bkg_img.use_flip_x = False
    bkg_img.use_flip_y = False

    # Show background images
    camera.show_background_images = True
