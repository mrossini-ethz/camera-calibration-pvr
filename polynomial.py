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

def make_poly(coeffs):
    """Make a new polynomial"""
    return list(coeffs)

def norm(poly):
    """Normalizes a given polynomial"""
    f = poly[-1]
    result = []
    for coeff in poly:
        result.append(coeff / f)
    return result

def sub(a, b):
    """Subtract the two polynomials"""
    n = max(len(a), len(b))
    _a = [0] * n
    _b = [0] * n
    for i in range(len(a)):
        _a[i] = a[i]
    for i in range(len(b)):
        _b[i] = b[i]
    result = []
    for i in range(n):
        result.append(_a[i] - _b[i])
    return result

def scale(poly, factor):
    """Normalizes a given polynomial"""
    f = poly[-1]
    result = []
    for coeff in poly:
        result.append(coeff * factor)
    return result

def reduce(poly):
    """Removes leading coefficients that are zero"""
    result = []
    for i in range(len(poly) - 1, -1, -1):
        if poly[i] != 0 or len(result) > 0:
            result.append(poly[i])
    result.reverse()
    return result

def derivative(poly):
    """Calculates the derivative of the polynomial"""
    result = []
    for i in range(1, len(poly)):
        result.append(i * poly[i])
    return result

def eval(poly, x):
    """Evaluate the polynomial"""
    result = 0.0
    for i in range(len(poly)):
        result += poly[i] * x ** i
    return result

def order(poly):
    """Get the order of the polynomial"""
    return len(poly) - 1

def coeff(poly, idx):
    """Get the nth coefficient of the polynomial"""
    if idx > len(poly) - 1:
        return 0.0
    elif idx >= 0:
        return poly[idx]

def div(a, b):
    """Calculate the polynom division of a and b"""
    na = order(a)
    nb = order(b)
    result = [0] * (na - nb + 1)
    for n in range(na, nb - 1, -1):
        f = a[n] / b[-1]
        result[n - nb] = f
        a = sub(a, [0] * (n - nb) + scale(b, f))
    return result
