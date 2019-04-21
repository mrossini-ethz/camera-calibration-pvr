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

from math import degrees
import mathutils
from . import polynomial
from . import rootfinder
from . import cameraplane
from . import transformation

def get_lambda_d_poly_a(qab, qac, qad, qbc, qbd, qcd):
    """Equation A (see paper)"""
    d4 = qac * qbd ** 2 - qad * qbc * qbd
    d3 = qab * qad * qbc + qad ** 2 * qbc * qbd + qad ** 2 * qcd + qbc * qbd - 2 * qab * qac * qbd - qab * qad * qbd * qcd - qac * qad * qbd ** 2
    d2 = qab ** 2 * qac + qab ** 2 * qad * qcd + 3 * qab * qac * qad * qbd + qab * qbd * qcd - qab * qad ** 2 * qbc - qab * qbc - qac * qad ** 2 - qad * qbc * qbd - 2 * qad * qcd
    d1 = qab * qad * qbc + 2 * qac * qad + qcd - 2 * qab ** 2 * qac * qad - qab ** 2 * qcd - qab * qac * qbd
    d0 = qab ** 2 * qac - qac
    return polynomial.make_poly([d0, d1, d2, d3, d4])

def get_lambda_d_poly_b(qab, qac, qad, qbc, qbd, qcd):
    """Equation B (see paper)"""
    d4 = qbd - qbd * qcd ** 2
    d3 = qab * qcd ** 2 + qac * qbd * qcd + 2 * qad * qbd * qcd ** 2 - qab - 2 * qad * qbd - qad * qbc * qcd
    d2 = 2 * qab * qad + qac * qad * qbc + qad ** 2 * qbc * qcd + qad **2 * qbd + qbc * qcd - qab * qac * qcd - qab * qad * qcd ** 2 - 3 * qac * qad * qbd * qcd - qbd * qcd ** 2
    d1 = qab * qac * qad * qcd + qac ** 2 * qad * qbd + 2 * qac * qbd * qcd - qab * qad ** 2 - qac * qad ** 2 * qbc - qac * qbc - qad * qbc * qcd
    d0 = qac * qad * qbc - qac ** 2 * qbd
    return polynomial.make_poly([d0, d1, d2, d3, d4])

def get_lambda_d(pa, pb, pc, pd, scale, focal_length):
    """Calculate the vectors from camera to the camera plane where the rectangle corners are located"""
    va = cameraplane.get_camera_plane_vector(pa, scale, focal_length).normalized()
    vb = cameraplane.get_camera_plane_vector(pb, scale, focal_length).normalized()
    vc = cameraplane.get_camera_plane_vector(pc, scale, focal_length).normalized()
    vd = cameraplane.get_camera_plane_vector(pd, scale, focal_length).normalized()
    # Calculate dot products
    qab = va.dot(vb)
    qac = va.dot(vc)
    qad = va.dot(vd)
    qbc = vb.dot(vc)
    qbd = vb.dot(vd)
    qcd = vc.dot(vd)
    # Determine the equation that needs to be solved
    pa = polynomial.norm(get_lambda_d_poly_a(qab, qac, qad, qbc, qbd, qcd))
    pb = polynomial.norm(get_lambda_d_poly_b(qab, qac, qad, qbc, qbd, qcd))
    print("A:", pa)
    print("B:", pb)
    p = polynomial.reduce(polynomial.sub(pa, pb))
    print("P:", p)
    # Solve the equation
    roots = rootfinder.find_poly_roots(p)
    print("Solutions:")
    # Iterate over all roots
    solutions = []
    for ld in roots:
        # Calculate the other parameters
        #ld = 1.10201
        lb = (qad * ld - 1) / (qbd * ld - qab)
        lc = (qad * ld - ld ** 2) / (qac - qcd * ld)
        # Scale the vectors pointing to the corners from the camera plane to 3d space
        ra = va
        rb = vb * lb
        rc = vc * lc
        rd = vd * ld
        # Printout for debugging
        print("x:", ld)
        # Corner angles
        angles = [degrees((rb - ra).angle(rd - ra)), degrees((ra - rb).angle(rc - rb)), degrees((rb - rc).angle(rd - rc)), degrees((rc - rd).angle(ra - rd))]
        print("Corner angles:", angles)
        # Rectangle size
        width = (rb - ra).length
        height = (rd - ra).length
        # Flatness (normal distance of point rd to plane defined by ra, rb, rc
        n = (ra - rb).cross(rc - rb)
        d = n.dot(ra)
        dist = abs(n.dot(rd) - d) / n.length
        print("Flatness:", dist, "=", dist / max(width, height) * 100, "%")
        # Calculate badness
        badness = 0.0
        # FIXME: angle badness and flatness badness should be weighted somehow
        for ang in angles:
            badness += abs(ang - 90)
        badness += abs(dist / max(width, height) * 100)
        print("Badness:", badness)
        solutions.append((badness, [ra, rb, rc, rd]))
    # Chose solution with best score
    best_badness = solutions[0][0]
    best_index = 0
    for i in range(1, len(solutions)):
        if best_badness > solutions[i][0]:
            best_index = i
            best_badness = solutions[i][0]
    # Return the best solution
    return solutions[best_index][1]

def reconstruct_rectangle(pa, pb, pc, pd, scale, focal):
    # Calculate the coordinates of the rectangle in 3d
    coords = get_lambda_d(pa, pb, pc, pd, scale, focal)
    # Calculate the transformation of the rectangle
    trafo = transformation.get_transformation(coords[0], coords[1], coords[2], coords[3])
    # Reconstruct the rotation angles of the transformation
    angles = transformation.get_rot_angles(trafo[0], trafo[1], trafo[2])
    xyz_matrix = mathutils.Euler((angles[0], angles[1], angles[2]), "XYZ")
    # Reconstruct the camera position and the corners of the rectangle in 3d such that it lies on the xy-plane
    tr = trafo[-1]
    cam_pos = transformation.apply_transformation([mathutils.Vector((0.0, 0.0, 0.0))], tr, xyz_matrix)[0]
    corners = transformation.apply_transformation(coords, tr, xyz_matrix)
    # Printout for debugging
    print("Focal length:", focal)
    print("Camera rotation:", degrees(angles[0]), degrees(angles[1]), degrees(angles[2]))
    print("Camera position:", cam_pos)
    length = (coords[0] - coords[1]).length
    width = (coords[0] - coords[3]).length
    size = max(length, width)
    print("Rectangle length:", length)
    print("Rectangle width:", width)
    print("Rectangle corners:", corners)
    return (cam_pos, xyz_matrix, corners, size)
