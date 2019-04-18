import bpy
import os

preview_collection = None

def load_images():
    global preview_collection

    import bpy.utils.previews
    preview_collection = bpy.utils.previews.new()

    my_icons_dir = os.path.join(os.path.dirname(__file__), "figures")

    preview_collection.load("1_point", os.path.join(my_icons_dir, "1-point.png"), 'IMAGE')
    preview_collection.load("2_point", os.path.join(my_icons_dir, "2-point.png"), 'IMAGE')
    preview_collection.load("3_point", os.path.join(my_icons_dir, "3-point.png"), 'IMAGE')
