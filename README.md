# Camera Calibration using Perspective Views of Rectangles

Camera Calibrabion PVR is a [Blender](http://www.blender.org) add-on for matching the 3D camera to the perspective seen in a given photograph.
A rectangle reference in the image is required for this calibration.

![Dragon Demo](https://raw.githubusercontent.com/wiki/mrossini-ethz/camera-calibration-pvr/images/dragon-demo.png)

Application example for this add-on.
The focal length, position and rotation of the camera are determined from a rectangle (graph paper) in the photograph.
The dragon (credit: [Stanford University](https://graphics.stanford.edu/data/3Dscanrep/)) and the array of cubes are rendered on top of the image using the calculated perspective.

### Blender Artist thread
[Link](https://blenderartists.org/forum/showthread.php?414359-Add-on-Camera-Calibration-using-Perspective-Views-of-Rectangles&p=3145913&viewfull=1#post3145913 "Link") to Blender Artist page

## Installation
1. Download the latest [release](https://github.com/mrossini-ethz/camera-calibration-pvr/releases) and save it in a directory of your convenience.
2. Open Blender.
3. In th menu go to Edit -> Preferences -> Addons.
4. At the top of the window, chose *Install*.
5. Select the file downloaded zip file and press *Install Addon*.
6. Search for *Camera: Camera Calibration using Perspective Views of Rectangles* in the Addon list.
7. Activate the checkbox for the plugin.
8. If you want to keep the addon activated when blender restarts, open the menu (bottom left menu button) and choose *Save Preferences*.

![Screenshot Installation](https://github.com/mrossini-ethz/camera-calibration-pvr/blob/master/doc/ui1.png "Screenshot Installation")

## Usage

### Addon Panel
The controls for the addon can be found in the Properties panel under Scene > Camera Calibration PVR.

### Mode overview
There are three modes that can be used to calibrate the camera:

- **Solve Focal**: Uses the image of a single rectangle to calibrate the camera. This is the mode most often used. However, it does not work if the photograph was taken with a tilt-shift lens and/or was cropped. Furthermore, it can not be used if any two edges of the rectangle appear parallel in the image.
- **Solve Focal+Y**: This is a special mode often used for architectural images where a tilt-shift lens was used. This mode assumes the lens was set up in such a way that vertical lines appear vertical in the image. Therefore, the reference rectangle is required to have exactly two parallel edges.
- **Solve Focal+X+Y**: This is the more general mode that allows for shift in both vertical and horizontal direction. This is often useful when the image was cropped. The rectangle used for reference may not have any parallel edges.
- **Solve 1-point**: This mode solves the 1-point perspective including camera shift. The reference rectangle size or the focal length of the camera must be known.

### Calibration of focal length, position and rotation (Solve Focal)
Perform the following steps:

1. Change the 3D view into *Top Orthographic* mode (Top view in Orthographic View mode: press `Numpad 7`, then `Numpad 5`).
2. In the 3D view, create a background image (`Shift-A` -> Image -> Background). Choose your reference image file.
3. In the properties panel under Scene > Camera Calibration PVR select the background image (default name: Empty).
3. Optional: Got to wireframe mode (press `Z` in the 3D view, then `4`).
4. Create a new Plane (`Shift-A` -> Mesh -> Plane).
5. Enter Edit mode (`Tab` key) and move the plane vertices to the corners of the rectangle in the reference image.
6. Leave Edit mode (`Tab` key).
7. Ensure that the plane is selected. In the Camera Calibration PVR panel press *Solve Focal*.

This sets the position, the rotation and the focal length of the camera.
A new object is created that represents the reconstructed rectangle.
By design, this rectangle is parallel to the x-y-plane, at the location of the 3D cursor.
(This can be changed in an option, see below.)
The name of the reconstructed rectangle is the one of the original plane, with "_Cal" appended.

### Calibration of focal length, vertical lens shift, position and rotation (Solve Focal+Y)
Vertical lens shift means that the optical centre is not in the centre of the image, but shifted in the vertical direction.
This can happen when using a tilt-shift camera lens or when cropping the image.
The algorithm used for the calibration in this case has two restrictions:

- The rectangle used for calibration needs to have two edges that are parallel
- The shift may only be along vertical direction in the image.

This applies very often in architectural photographs.

To perform the calibration, go through the following steps:

1. Change the 3D view into *Top Orthographic* mode (Top view in Orthographic View mode: press `Numpad 7`, then `Numpad 5`).
2. In the 3D view, create a background image (`Shift-A` -> Image -> Background). Choose your reference image file.
3. In the properties panel under Scene > Camera Calibration PVR select the background image (default name: Empty).
4. Optional: Got to wireframe mode (press `Z` in the 3D view, then `4`).
5. Create a new Plane (`Shift-A` -> Mesh -> Plane).
6. Enter Edit mode (`Tab` key) and move the plane vertices to the corners of the rectangle in the reference image. Ensure that two sides are perfectly parallel.
7. Select one of the vertices and extrude it by pressing `E` as shown in the image below. In the 3D space of the reference image, the new edge must be perpendicular to the rectangle.
8. Leave Edit mode (`Tab` key).
9. Ensure that the plane is selected. In the Camera Calibration PVR panel press *Solve Focal+Y*.

![Screenshot Usage](https://github.com/mrossini-ethz/camera-calibration-pvr/blob/master/doc/shifted-perspective.png "Screenshot: Usage for Solve Focal+Y")

### Calibration of focal length, horizontal and vertical lens shift, position and rotation (Solve Focal+X+Y)
Lens shift means that the optical centre is not in the centre of the image, but shifted off-centre.
This can happen when using a tilt-shift camera lens or when cropping the image.

To perform the calibration, go through the following steps:

1. Change the 3D view into *Top Orthographic* mode (Top view in Orthographic View mode: press `Numpad 7`, then `Numpad 5`).
2. In the 3D view, create a background image (`Shift-A` -> Image -> Background). Choose your reference image file.
3. In the properties panel under Scene > Camera Calibration PVR select the background image (default name: Empty).
4. Optional: Got to wireframe mode (press `Z` in the 3D view, then `4`).
5. Create a new Plane (`Shift-A` -> Mesh -> Plane).
6. Enter Edit mode (`Tab` key) and move the plane vertices to the corners of the rectangle in the reference image.
7. Select one of the vertices and extrude it by pressing `E`. In the 3D space of the reference image, the new edge must be perpendicular to the rectangle.
7. Repeat the last step for a second vertex.
8. Leave Edit mode (`Tab` key).
9. Ensure that the plane is selected. In the Camera Calibration PVR panel press *Solve Focal+X+Y*.

### Options
There are options to the camera calibration:

- **Size** determines the size of the reconstructed rectangle.
- **Vertical orientation** places the reconstructed rectangle in vertical orientation instead of parallel to the x-y-Plane.

### Demos
Here are screen casts that demonstrate the usage.

#### Focal
![Screen Cast Usage](https://raw.githubusercontent.com/wiki/mrossini-ethz/camera-calibration-pvr/images/focal.gif "Screen cast: Usage for Solve Focal")

#### Focal+Y
![Screen Cast Usage](https://raw.githubusercontent.com/wiki/mrossini-ethz/camera-calibration-pvr/images/focal+y.gif "Screen cast: Usage for Solve Focal+Y")

## Things you should know
You should be aware of a few things when using this add-on.

### Technical aspects of reference images
- **Solve Focal:** Images used for camera calibration should have the optical centre in the centre of the image.
  This means that cropping the image unevenly or using tilt-shift lenses will make the add-on fail to work.
  **Solve Focal+Y:** See the restrictions above.
- Scaling of the image prior to the use for calibration should be done only if the aspect ratio of the image is preserved.
- Lens distortion effects will negatively affect the result.
- Higher image resolution is beneficial to the accuracy of the calibration.

### Reference image content
- Rectangles in the image used for calibration should be distorted by perspective.
  **Solve Focal**, **Solve Focal+X+Y**: The more parallel the sides of the rectangles appear in the image, the worse the result.
  Completely parallel sides will not work. (See also **Solve Focal+Y**.)
- Right angles in the real world are often not perfect.

## Changelog

### Version 0.5
- Added 1-point perspective support.
- Fixed small technical issues that appeared with newer versions of blender.

### Version 0.4
- Improved the add-on preferences.
- Moved the tool location from the "Misc" tab to the "Tools" tab in the Tool
  shelf. The location can be configured now.
- The scene name is not hard coded anymore.

### Version 0.3
- Added an algorithm for calculating vertical and horizontal lens shift along
  with focal length, position and rotation of the camera. Additional
  information is taken from two dangling vertices in the mesh.
- Separate button for the new calibration method.
- Renaming of the calibration methods.

### Version 0.2
- Added an algorithm for calculating vertical lens shift along with focal
  length, position and rotation of the camera. Additional information is taken
  from one dangling vertex in the mesh. The rectangle in the image is required
  to have one pair of parallel edges.
- Separate button for the new calibration method.
- Adding a camera to the scene when none exist.
- Bugfixes

### Version 0.1
First official release. It includes the following features:

- Algorithm to calculate focal length, position and rotation for the camera
  used to make a picture of a rectangle.
- The image used for calibration is taken from the viewport background from Top
  View.
- A mesh of four vertices in one polygon is used to determine the coordinates
  of the rectangle corners.
- A button in the tool shelf performs the calibration when the mesh is
  selected.
- A new rectangle is created during the calibration which represents the
  undistorted rectangle in the image.
- The active camera is transformed according to the calculation results.
- The view automatically changes to camera view.
- Options for the calibration are:
    - the size of the reconstructed rectangle
    - vertical alignment of the reconstructed rectangle

## License
Camera Calibration with Perspective Views of Rectangles

Copyright (C) 2017  Marco Rossini

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
version 2 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

This Blender plugin is based on the research paper "Recovery of Intrinsic
and Extrinsic Camera Parameters Using Perspective Views of Rectangles" by
T. N. Tan, G. D. Sullivan and K. D. Baker, Department of Computer Science,
The University of Reading, Berkshire RG6 6AY, UK,
from the Proceedings of the British Machine Vision Conference, published by
the BMVA Press.
