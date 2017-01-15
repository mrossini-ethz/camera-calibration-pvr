# Camera Calibration using Perspective Views of Rectangles

Camera-Calibrabion-PVR is a [Blender](http://www.blender.org) plugin for calibrating the 3d camera using a single image of a rectangle.

## Installation
1. Clone or download the repository into a directory of your convenience.
2. Open Blender.
3. Goto File -> User Preferences -> Addons.
4. At the bottom of the window, chose Install From File.
5. Select the file camera-calibration-pvr.py from the directory into which you cloned the repository.
6. Activate the checkbox for the plugin that you will now find in the list.

![Screenshot Installation](https://github.com/mrossini-ethz/camera-calibration-pvr/blob/master/doc/ui1.png "Schreenshot Installation")

## Usage
Ensure that your scene has a camera. Then perform the following steps:
1. Change the 3D view into Top Ortho mode (Top view in Orthographic View mode: press `Numpad 7`, then `Numpad 5`).
2. Load a background image for the 3D View by the 3D view properties (on the right side of the 3D view). Ensure that the view axis 'Top' is selected.
3. Create a new Plane (`Shift-A` -> Mesh -> Plane).
4. Enter Edit mode (`Tab` key) and move the plane vertices to the corners of the rectangle in the background image.
5. Leave Edit mode (`Tab` key).
6. Enusre that the plane is selected. In the 3D view tools (on the left side of the 3D view) select the Misc tab and find the Camera Calibration menu. Press Camera Calibration.

This sets the position, the rotation and the focal length of the camera. The view changes automatically to the camera and a new object called CalRect is created that represents the reconstructed rectangle. By design, this rectangle is flat on the x-y-plane.

To reposition, rotate or scale the generated rectangle, be sure to reposition, rotate or scale the camera along with it.

![Screenshot Usage](https://github.com/mrossini-ethz/camera-calibration-pvr/blob/master/doc/ui2.png "Schreenshot Usage")

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