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

from math import sqrt
from . import polynomial

def find_root(f, df, ddf, initial_guess = 0.0, limit = 0.00001, max_iterations = 1000):
    """Find the root of the function f using Halley's method"""
    xn_1 = initial_guess
    i = 0
    while i < max_iterations:
        fx = f(xn_1)
        dfx = df(xn_1)
        ddfx = ddf(xn_1)
        xn = xn_1 - 2 * fx * dfx / (2 * dfx ** 2 - fx * ddfx)
        if abs(xn - xn_1) < limit:
            return xn
        xn_1 = xn
        i += 1
    return None

def find_poly_root(poly, initial_guess = 0.0, limit = 0.00001, max_iterations = 1000):
    """Find a root of the given polynomial"""
    # Calculate the polynomial derivatives
    dpoly = polynomial.derivative(poly)
    ddpoly = polynomial.derivative(dpoly)
    # Closures !!!
    f = lambda x: polynomial.eval(poly, x)
    df = lambda x: polynomial.eval(dpoly, x)
    ddf = lambda x: polynomial.eval(ddpoly, x)
    # Call the generic root finder
    return find_root(f, df, ddf, initial_guess, limit, max_iterations)

def find_poly_roots(poly, initial_guess = 0.0, limit = 0.00001, max_iterations = 1000):
    """Find all roots of the given polynomial"""
    solutions = []
    # Find solutions numerically for n > 0, split them off until n = 2
    for q in range(polynomial.order(poly) - 2):
        x = find_poly_root(poly, initial_guess, limit, max_iterations)
        if not x:
            break
        poly = polynomial.div(poly, polynomial.make_poly([-x, 1]))
        solutions.append(x)
    # Find the rest of the roots analytically
    if polynomial.order(poly) == 1:
        solutions.append(- polynomial.coeff(poly, 1) / polynomial.coeff(poly, 0))
    elif polynomial.order(poly) == 2:
        a = polynomial.coeff(poly, 2)
        b = polynomial.coeff(poly, 1)
        c = polynomial.coeff(poly, 0)
        d = b ** 2 - 4 * a * c
        if d == 0:
            solutions.append(-b / (2 * a))
        elif d > 0:
            solutions.append((- b + sqrt(d)) / (2 * a))
            solutions.append((- b - sqrt(d)) / (2 * a))
    return solutions
