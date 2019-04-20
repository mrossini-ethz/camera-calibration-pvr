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

def solve_linear_system_2d(a, b, c, d, e, f):
    """Solves the system of equations a*x + b*y = c, d*x + e*y = e using Gaussian Elimination (for numerical stability)."""
    # Pivoting (to obtain stability)
    if abs(d) > abs(a):
        a, b, c, d, e, f = (d, e, f, a, b, c)
    # Check for singularity
    if a == 0:
        return None
    tmp = e - d * b / a
    if tmp == 0:
        return None
    # This is final answer of the gaussian elimination
    y = (f - d * c / a) / tmp
    x = (c - b * y) / a
    return (x, y)
