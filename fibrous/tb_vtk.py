#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)

import numpy as nm
import yaml
import argparse
import sys

if sys.version_info.major == 3:
    xrange = range
import numpy as np
import vtk
from scipy.interpolate import InterpolatedUnivariateSpline
from . import tree

# new interface


class TBVTK(tree.TubeSkeletonBuilder):
    """
    This generator is called by generateTree() function as a general form.
    Other similar generator is used for generating LAR outputs.
    """

    def __init__(self,
                 # gtree,
                 cylinder_resolution=30, sphere_resolution=30,
                 polygon_radius_selection_method="inscribed",
                 cylinder_radius_compensation_factor=1.0,
                 sphere_radius_compensation_factor=1.0,
                 # tube_shape=True
                 tube_shape=False
                 ):
        """

        :param gtree:
        :param cylinder_resolution:
        :param sphere_resolution:
        :param polygon_radius_selection_method:
        :param cylinder_radius_compensation_factor:
        :param sphere_radius_compensation_factor:
        :param tube_shape: If true, the tube shape is generated.
                Otherwise the cylinder shape is used.
        """
        # self.shape = gtree.shape
        # self.data3d = np.zeros(gtree.shape, dtype=np.int)
        # self.voxelsize_mm = gtree.voxelsize_mm
        # super(tree.FiberSkeletBuilder, self).__init__()
        tree.TubeSkeletonBuilder.__init__(self)
        # make comapatible with old system
        self.polygon_radius_selection_method = polygon_radius_selection_method
        self.cylinder_radius_compensation_factor = cylinder_radius_compensation_factor
        self.sphere_radius_compensation_factor = sphere_radius_compensation_factor
        #self.tree_data = gtree.tree_data

        self.cylinder_resolution = cylinder_resolution
        self.sphere_resolution = sphere_resolution
        logger.debug("tube shape " + str(tube_shape))
        self.tube_shape = tube_shape

    def add_cylinder(self, p1m, p2m, rad, id):
        """
        Funkce na vykresleni jednoho segmentu do 3D dat
        """
        pass

    def finish(self):
        self.tube_skeleton_old = compatibility_processing(self.tube_skeleton)
        # import ipdb; ipdb.set_trace()
        self.polyData = gen_tree(
            self.tube_skeleton_old, self.cylinder_resolution, self.sphere_resolution,
            polygon_radius_selection_method=self.polygon_radius_selection_method,
            cylinder_radius_compensation_factor=self.cylinder_radius_compensation_factor,
            sphere_radius_compensation_factor=self.sphere_radius_compensation_factor,
            tube_shape=self.tube_shape
        )
        # import ipdb; ipdb.set_trace()

    def get_output(self):
        return self.polyData

    def save(self, outputfile, lc_all="C"):

        import vtk
        logger.debug("vtk version " + str(vtk.VTK_BUILD_VERSION))
        # import ipdb; ipdb.set_trace()
        if lc_all is not None:
            import locale
            locale.setlocale(locale.LC_ALL, lc_all)
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName(outputfile)
        try:
            writer.SetInputData(self.polyData)
        except:
            logger.warning("old vtk is used")
            writer.SetInput(self.polyData)
        writer.Write()

    def show(self):
        logger.info("there is no show implemented")


# old interface

def move_to_position(src, upper, direction, axis0=2, axis1=1, axis2=0):

    # along axis 1 axis1=2, axis2=1

    r1 = [0.0, 0.0, 0.0]
    r2 = [0.0, 0.0, 0.0]

    r1[axis0] = 1.0
    r2[axis1] = 1.0


    rot1 = vtk.vtkTransform()
    fi = nm.arccos(direction[axis1])

    rot1.RotateWXYZ(-nm.rad2deg(fi), r1[0], r1[1], r1[2])
    u = nm.abs(nm.sin(fi))
    rot2 = vtk.vtkTransform()
    if u > 1.0e-6:

        # sometimes d[0]/u little bit is over 1
        d0_over_u = direction[axis2] / u
        if d0_over_u > 1:
            psi = 0
        elif d0_over_u < -1:
            psi = 2 * nm.pi
        else:
            psi = nm.arccos(direction[axis2] / u)

        # logger.debug('d0 ' + str(direction[axis2]) + '  u ' + str(u) + ' psi ' + str(psi))
        if direction[axis0] < 0:
            psi = 2 * nm.pi - psi

        rot2.RotateWXYZ(-nm.rad2deg(psi), r2[0], r2[1], r2[2])

    tl = vtk.vtkTransform()
    tl.Translate(upper)

    tr1a = vtk.vtkTransformFilter()
    if "SetInputConnection" in dir(tr1a):
        tr1a.SetInputConnection(src.GetOutputPort())
    else:
        tr1a.SetInput(src.GetOutput())
    tr1a.SetTransform(rot1)

    tr1b = vtk.vtkTransformFilter()
    if "SetInputConnection" in dir(tr1b):
        tr1b.SetInputConnection(tr1a.GetOutputPort())
    else:
        tr1b.SetInput(tr1a.GetOutput())
    # tr1b.SetInput(tr1a.GetOutput())
    tr1b.SetTransform(rot2)

    tr2 = vtk.vtkTransformFilter()
    if "SetInputConnection" in dir(tr2):
        tr2.SetInputConnection(tr1b.GetOutputPort())
    else:
        tr2.SetInput(tr1b.GetOutput())
    # tr2.SetInput(tr1b.GetOutput())
    tr2.SetTransform(tl)

    tr2.Update()

    return tr2
    # return tr2.GetOutput()


def get_tube(radius=1.0, point=[0.0, 0.0, 0.0],
             direction=[0.0, 0.0, 1.0], length=1.0,
             sphere_resolution=10, cylinder_resolution=10,
             cylinder_radius_compensation_factor=1.0,
             sphere_radius_compensation_factor=1.0,
             tube_shape=True, axis=1,
             make_ladder_even=True
             ):
    point1 = [0.0, 0.0, 0.0]
    center = [0.0, 0.0, 0.0]
    point2 = [0.0, 0.0, 0.0]

    center[axis] = length / 2.0
    point2[axis] = length

    cylinder_radius = radius * cylinder_radius_compensation_factor
    sphere_radius = radius * sphere_radius_compensation_factor

    direction /= nm.linalg.norm(direction)
    lv = point + direction * length

    cylinderTri = vtk.vtkTriangleFilter()
    sphere1Tri = vtk.vtkTriangleFilter()
    sphere2Tri = vtk.vtkTriangleFilter()

    cylinder = vtk.vtkCylinderSource()
    cylinder.SetCenter(center)
    cylinder.SetHeight(length)
    cylinder.SetRadius(cylinder_radius)
    cylinder.SetResolution(cylinder_resolution)
    cylinder.Update()
    cylinderTri.SetInputData(cylinder.GetOutput())
    cylinderTri.Update()

    # make ladder even
    if make_ladder_even:
        if sphere_resolution % 2 == 0:
            phi_resolution = sphere_resolution + 1
        else:
            phi_resolution = sphere_resolution

    if not tube_shape:
        tube = move_to_position(cylinderTri, point, direction, 2, 1, 0)
        return tube.GetOutput()

    sphere1 = get_sphere(
        center=point1,
        radius=sphere_radius,
        resolution=sphere_resolution,
        start_phi=0,
        #end_phi=90,
        end_phi=180,
        axis=1,
        phi_resolution=phi_resolution
    )
    # sphere1.Update()

    sphere1Tri.SetInputData(sphere1)
    sphere1Tri.Update()

    sphere2 = get_sphere(
        center=point2,
        # radius= 1. - (cylinder_radius - sphere_radius),
        radius=sphere_radius,
        resolution=sphere_resolution,
        start_phi=0,
        end_phi=180,
        axis=1,
        phi_resolution=phi_resolution

    )
    sphere2Tri.SetInputData(sphere2)
    sphere2Tri.Update()

    boolean_operation1 = vtk.vtkBooleanOperationPolyDataFilter()
    boolean_operation2 = vtk.vtkBooleanOperationPolyDataFilter()
    boolean_operation1.SetOperationToUnion()
    boolean_operation2.SetOperationToUnion()

    # booleanOperation.SetInputData(0, cyl)
    boolean_operation1.SetInputData(0, cylinderTri.GetOutput())
    boolean_operation1.SetInputData(1, sphere1Tri.GetOutput())
    boolean_operation1.Update()
    boolean_operation2.SetInputData(0, boolean_operation1.GetOutput())
    boolean_operation2.SetInputData(1, sphere2Tri.GetOutput())
    # booleanOperation.SetInputData(2, sph2)
    boolean_operation2.Update()
    # tube_in_base_position = boolean_operation2.GetOutput()

    #tube = move_to_position(boolean_operation2, point, direction, 1, 2)
    tube = move_to_position(boolean_operation2, point, direction, 2, 1, 0)
    return tube.GetOutput()


def get_cylinder(upper, height, radius,
                 direction,
                 resolution=10):
    import vtk
    src = vtk.vtkCylinderSource()
    src.SetCenter((0, height / 2, 0))
    # src.SetHeight(height + radius/2.0)
    src.SetHeight(height)
    src.SetRadius(radius)
    src.SetResolution(resolution)
    return move_to_position(src, upper, direction).GetOutput()


def get_sphere(center, radius, resolution=10, start_phi=None, end_phi=None, axis=0, **kwargs):
    sph = get_sphere_source(center, radius, resolution, start_phi, end_phi, axis, **kwargs)
    #sph.Update()
    return sph

def get_sphere_source(center, radius, resolution=10, start_phi=None, end_phi=None, axis=0, theta_resolution=None, phi_resolution=None):
    # create source
    if theta_resolution is None:
        theta_resolution=resolution
    if phi_resolution is None:
        phi_resolution=resolution
    import vtk
    sphere = vtk.vtkSphereSource()
    sphere.SetPhiResolution(phi_resolution)
    sphere.SetThetaResolution(theta_resolution)
    # sphere.SetCenter(center[0], center[1], center[2])
    sphere.SetCenter(.0, .0, .0)
    sphere.SetRadius(radius)
    if start_phi is not None:
        sphere.SetStartPhi(start_phi)
    if end_phi is not None:
        sphere.SetEndPhi(end_phi)

    if axis == 0:
        translate = vtk.vtkTransform()
        translate.Translate(center)

        tr2 = vtk.vtkTransformFilter()
        tr2.SetInputConnection(sphere.GetOutputPort())
        tr2.SetTransform(translate)
        sphere = tr2

    if axis == 1:
        #sphere = move_to_position(sphere, center, [1., 1., 0.])

        rot1 = vtk.vtkTransform()
        rot1.RotateWXYZ(90, 1, 0, 0)
        translate = vtk.vtkTransform()
        translate.Translate(center)

        tr1 = vtk.vtkTransformFilter()
        # tr1a.SetInputConnection(src.GetOutputPort())
        tr1.SetInputConnection(sphere.GetOutputPort())
        tr1.SetTransform(rot1)
        tr1.Update()

        tr2 = vtk.vtkTransformFilter()
        tr2.SetInputConnection(tr1.GetOutputPort())
        tr2.SetTransform(translate)
        sphere = tr2
    sphere.Update()
    return sphere.GetOutput()

def polygon_radius_compensation_factos(
        polygon_radius_selection_method,
        cylinder_radius_compensation_factor,
        sphere_radius_compensation_factor,
        cylinder_resolution,
        sphere_resolution
):
    # cylinder volume + sphere error
    # x_cvse = [6, 8, 10, 12, 14, 17, 21, 23, 29, 31, 39, 46, 50, 60, 70, 80, 90, 100, 150, 200]
    # y_cvse = [0.907761087455, 0.949394472294, 0.968085949491, 0.978051451209, 0.983984947219, 0.989216674476, 0.992978216638, 0.994159974337, 0.996345087831, 0.996805451193, 0.997989063564, 0.998557649514, 0.998780374797, 0.999154602003, 0.999379702886, 0.999525551526, 0.999625411186, 0.999696768773, 0.999865478984, 1.0]
    x_cvse = [5, 6, 7, 8, 10, 12, 14, 17, 19, 21, 23, 29, 31, 34, 36, 39, 42, 44, 46, 48, 50, 53, 58, 60, 63, 67, 70, 73, 78, 80, 84, 87, 90, 95, 200, ]
    y_cvse = [0.86445499173, 0.907761069358, 0.933196212368, 0.94939448372, 0.968085943193, 0.978051440662, 0.983984934568, 0.9892166767, 0.991397728183, 0.992978212333, 0.994159975923, 0.996345087892, 0.996805450921, 0.997348549124, 0.997637128575, 0.997989061078, 0.998267848559, 0.998422687049, 0.998557649457, 0.998675999854, 0.998780372347, 0.998915211425, 0.999094999814, 0.999154603983, 0.999233534268, 0.999322673214, 0.999379704358, 0.999429824527, 0.999500820765, 0.999525551368, 0.999569806342, 0.999599052092, 0.99962541127, 0.999663914279, 1.0, ]
    # cylinder surface + sphere error
    # x_csse = [6, 8, 10, 12, 14, 17, 21, 23, 29, 31, 39, 46, 50, 60, 70, 80, 90, 100, 150, 200]
    # y_csse = [0.975228602567, 0.987423714247, 0.992397574304, 0.994910516827, 0.996355306056, 0.997593221093, 0.998459193066, 0.998726453568, 0.999213600518, 0.999314918173, 0.999572948543, 0.999695445089, 0.999743138504, 0.99982282109, 0.999870450728, 0.999901168474, 0.999922129169, 0.999937063217, 0.999972213153, 1.0]
    # x_csse = [5, 6, 7, 8, 10, 12, 14, 17, 19, 21, 23, 29, 31, 34, 36, 39, 42, 44, 46, 48, 50, 53, 58, 60, 63, 67, 70, 73, 78, 80, 84, 87, 90, 95, 200, ]
    # y_csse = [0.960899219376, 0.975228595976, 0.982858510101, 0.987423713758, 0.992397573545, 0.994910509499, 0.996355309909, 0.99759321248, 0.998098030942, 0.998459189081, 0.998726456981, 0.999213598605, 0.999314919675, 0.999433785576, 0.999496625338, 0.999572950856, 0.999633144062, 0.999666470021, 0.999695450011, 0.999720810685, 0.999743136717, 0.999771914668, 0.999810174968, 0.999822823213, 0.999839549217, 0.9998584046, 0.999870450859, 0.999881016339, 0.999895967889, 0.999901168437, 0.999910466349, 0.999916604662, 0.999922127085, 0.999930191796, 1.0, ]
    # lot of paper figure is based on thhis setting 2017-07
    # x_csse = [5, 6, 7, 8, 10,
    #           12, 14, 17, 19,
    #           21, 23, 29, 31,
    #           34, 36, 39, 42,
    #           44, 46, 48, 50, 53, 58, 60, 63, 67, 70, 73, 78, 80, 84, 87, 90, 95, 200, ]
    # y_csse = [0.898910918632, 0.931274723443, 0.950194578295, 0.962239811246, 0.976153670948,
    #           0.983584378196, 0.988014411793, 0.991924785771, 0.993556301246, 0.994739096095,
    #           0.995623800682, 0.997260361529, 0.997605276155, 0.998012238947, 0.998228507342,
    #           0.998492277871, 0.998701247359, 0.998817315734, 0.998918487886, 0.999007210034,
    #           0.999085459703, 0.999186550556, 0.999321349616, 0.999366042071, 0.999425223485, 0.999492061445, 0.999534825661, 0.999572408455, 0.999625645831, 0.999644190807, 0.999677378202, 0.999699309004, 0.99971907481, 0.999747949314, 1.0, ]
    x_csse = [5, 6, 7, 8, 10, 12, 14, 17, 19, 21, 23, 29, 30, 31, 34, 36, 39, 42, 44, 46, 48, 50, 53, 58, 60, 63, 67,
              70, 73, 78, 80, 84, 87, 90, 95, 200, ]
    y_csse = [0.900576771582, 0.936869836439, 0.950985318762, 0.964590043529, 0.977386494167, 0.984323254476,
              0.988498485354, 0.992048202584, 0.993654556051, 0.994819174883, 0.995690300838, 0.997301856778,
              0.997508981673, 0.997641534309, 0.998061908546, 0.998271760101, 0.998515046246, 0.998731160009,
              0.998844108505, 0.998942638136, 0.999029070544, 0.999105325336, 0.99919881978, 0.999335448082,
              0.99937908559, 0.99943387055, 0.999500058848, 0.999544024835, 0.999578847301, 0.999632860646,
              0.999651016253, 0.99968349269, 0.999703829344, 0.999724328432, 0.999751740122, 1.0, ]

    # x_csseje = [6, 8, 10, 12, 16, 20, 22, 23, 26, 27, 29, 33, 37, 41, 51, 57, 61, 67, 73, 77, 81, 87, 100, 200, ]
    # y_csseje = [1.01816186937, 1.0105899069, 1.00695437881, 1.00642661083, 1.00279254142, 1.00199576241, 1.0011686604, 1.00019346955, 1.00174728137, 1.00013695628, 1.00011940911, 1.00009067128, 1.0000715922, 1.00005867129, 1.00003740911, 1.00002975563, 1.00002600566, 1.00002147955, 1.00001804574, 1.00001619267, 1.00001462116, 1.00001264225, 1.0, 1.0, ]
    x_csseje = [5, 7, 9, 11, 13, 15, 19, 23, 25, 27, 29, 33, 35, 37, 41, 43, 45, 51, 55, 57, 61, 65, 67, 69, 73, 77, 81, 85, 87, 100, 200, ]
    y_csseje = [1.00280461843, 1.0008915901, 1.0003335719, 1.00019869299, 1.00011629061, 1.00007253692, 1.00003596271, 1.00003716463, 1.00002536491, 1.00002053443, 1.00001996422, 1.00002226941, 1.00002018493, 1.00001772083, 1.00001419362, 1.00001265241, 1.00001143422, 1.00000843313, 1.0000074585, 1.00000670063, 1.00000585268, 1.00006477878, 1.00006803908, 1.00007184803, 1.00006644751, 1.0000630609, 1.00005628084, 1.00005190342, 1.00007740341, 1.0, 1.0, ]

    # x_cvseje = [5, 7, 9, 11, 13, 15, 19, 23, 25, 27, 29, 33, 35, 37, 41, 43, 45, 51, 55, 57, 61, 65, 67, 69, 73, 77, 81, 85, 87, 91, 95, 200, ]
    # y_cvseje = [1.00013588959, 1.00000045991, 0.999781214117, 1.00002153432, 0.999984546002, 0.999986251724, 1.00000199896, 0.999989169617, 0.999974198814, 0.999985921625, 0.9999951955, 0.999998502971, 1.00000178949, 1.00000359342, 0.999999259851, 0.999998648343, 0.999997430649, 1.00000576793, 1.00000262582, 1.00000426568, 1.00000121655, 1.0000049034, 1.00000506061, 1.00000511941, 1.00000306134, 1.00000425606, 1.00000622994, 1.00000603282, 1.00000469495, 1.00000393567, 1.00000215059, 1.0, ]
    x_cvseje = [5, 7, 9, 11, 13, 15, 19, 23, 25, 27, 29, 33, 35, 37, 41, 43, 45, 51, 55, 57, 61, 65, 67, 69, 73, 77, 81, 85, 200, ]
    y_cvseje = [1.00016060109, 1.00000045991, 0.999781214117, 1.00002153432, 0.999984546002, 0.999986251724, 1.00000199896, 0.999989169617, 0.999974198814, 0.999985921625, 0.9999951955, 0.999998502971, 1.00000178949, 1.00000359342, 0.999999259851, 0.999998648343, 0.999997430649, 1.00000576793, 1.00000262582, 1.00000426568, 1.00000121655, 1.0000049034, 1.00000506061, 1.00000511941, 1.00000306134, 1.00000425606, 1.00000622994, 1.00000603282, 1.0, ]

    if polygon_radius_selection_method == "inscribed":
        cylinder_radius_compensation_factor = 1.0
        sphere_radius_compensation_factor = 1.0
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "circumscribed":
        # from .. import geometry3d as g3
        factor = circumscribed_polygon_radius(cylinder_resolution)
        cylinder_radius_compensation_factor = factor
        sphere_radius_compensation_factor = factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "compensation factors":
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor
        pass

    elif polygon_radius_selection_method == "cylinder surface":
        # from .. import geometry3d as g3
        radius_compensation_factor =  regular_polygon_perimeter_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "cylinder volume":
        # from .. import geometry3d as g3
        radius_compensation_factor =  regular_polygon_area_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor = radius_compensation_factor
        sphere_radius_compensation_factor = radius_compensation_factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "average":
        # from .. import geometry3d as g3
        factor = circumscribed_polygon_radius(cylinder_resolution)
        factor = (factor + 1.0) / 2.0
        cylinder_radius_compensation_factor = factor
        sphere_radius_compensation_factor = factor
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor


    elif polygon_radius_selection_method == "cylinder volume + sphere error":
        # analytically compensated cylinder + sphere compensate by measurement
        cylinder_radius_compensation_factor *=  regular_polygon_area_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor

        spl1 = InterpolatedUnivariateSpline(x_cvse, y_cvse)
        sphere_radius_compensation_factor *= 1. / spl1(cylinder_resolution)
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "cylinder surface + sphere error":
        # analytically compensated cylinder + sphere compensate by measurement
        cylinder_radius_compensation_factor =  regular_polygon_perimeter_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor

        spl1 = InterpolatedUnivariateSpline(x_csse, y_csse)
        sphere_radius_compensation_factor *= 1. / spl1(cylinder_resolution)
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor

    elif polygon_radius_selection_method == "cylinder surface + sphere error + join error":
        # sphere like objects
        # analytically compensated cylinder + sphere compensate by measurement

        # new
        # analytically compensated cylinder + sphere compensate by measurement
        cylinder_radius_compensation_factor *=  regular_polygon_area_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor

        spl1 = InterpolatedUnivariateSpline(x_csse, y_csse)
        sphere_radius_compensation_factor = 1. / spl1(cylinder_resolution)

        spl2 = InterpolatedUnivariateSpline(x_csseje, y_csseje)
        # radius_compensation_factor *= 1. / spl1(cylinder_resolution)
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor / spl2(cylinder_resolution)

    elif polygon_radius_selection_method == "cylinder volume + sphere error + join error":
        # sphere like objects
        # analytically compensated cylinder + sphere compensate by measurement
        cylinder_radius_compensation_factor =  regular_polygon_area_equivalent_radius(cylinder_resolution)
        cylinder_radius_compensation_factor_long = cylinder_radius_compensation_factor

        spl1 = InterpolatedUnivariateSpline(x_cvse, y_cvse)
        sphere_radius_compensation_factor *= 1. / spl1(cylinder_resolution)

        spl2 = InterpolatedUnivariateSpline(x_cvseje, y_cvseje)
        sphere_radius_compensation_factor_long = sphere_radius_compensation_factor / spl2(cylinder_resolution)

    else:
        logger.error("Unknown compensation method '{}'".format(polygon_radius_selection_method))
        print("Unknown compensation method '{}'".format(polygon_radius_selection_method))

    return cylinder_radius_compensation_factor, sphere_radius_compensation_factor, cylinder_radius_compensation_factor_long, sphere_radius_compensation_factor_long


def gen_tree(tree_data, cylinder_resolution=10, sphere_resolution=10,
             polygon_radius_selection_method="inscribed",
             cylinder_radius_compensation_factor=1.0,
             sphere_radius_compensation_factor=1.0,
             tube_shape=True
             ):
    """

    :param polygon_radius_selection_method:
        "inscribed":
        "compensation factors"
        "cylinder volume"
        "cylinder surface"
    :param tree_data:
    :param cylinder_resolution:
    :param sphere_resolution:
    :param cylinder_radius_compensation_factor: is used to change radius of cylinder and spheres
    :return:
    """
    import vtk
    # appendFilter = vtk.vtkAppendPolyData()
    appended_data = None
    if vtk.VTK_MAJOR_VERSION <= 5:
        logger.error("VTK 6 required")
    factors = polygon_radius_compensation_factos(
        polygon_radius_selection_method,
        cylinder_radius_compensation_factor,
        sphere_radius_compensation_factor,
        cylinder_resolution,
        sphere_resolution
    )

    cylinder_radius_compensation_factor, sphere_radius_compensation_factor,\
    cylinder_radius_compensation_factor_long, sphere_radius_compensation_factor_long = factors

    # import ipdb; ipdb.set_trace()
    print("len tree_data", len(tree_data))
    for br in tree_data[:]:
        # import ipdb;
        # ipdb.set_trace()
        something_to_add = True
        radius = br['radius']
        length = br["length"]
        direction = br["direction"]
        uv = br['upperVertex']
        print(length, radius, uv, direction, nm.linalg.norm(direction))
        if direction[1] == 0:
            print("nula", direction)
            # direction[1] = -0.1
            # direction[0] = 1
            # direction[1] = 0
            # direction[2] = 0
            # direction = direction / np.linalg.norm(direction)
            # direction[1] = -0.1

        dbg_msg = "generating edge with length: " + str(br["length"])
        logger.debug(dbg_msg)
        # tube = get_tube_old(radius, uv, direction, length,
        if length == 0:
            tube = get_sphere(uv, radius * sphere_radius_compensation_factor, sphere_resolution)
        else:
            tube = get_tube(radius, uv, direction, length,
                            sphere_resolution, cylinder_resolution,
                            cylinder_radius_compensation_factor=cylinder_radius_compensation_factor_long,
                            sphere_radius_compensation_factor=sphere_radius_compensation_factor_long,
                            tube_shape=tube_shape)
            print("get_tube", something_to_add)
        # this is simple version
        # appendFilter.AddInputData(boolean_operation2.GetOutput())
        # print "object connected, starting addind to general space " + str(br["length"])
        if something_to_add:
            if appended_data is None:
                # appended_data = boolean_operation2.GetOutput()
                appended_data = tube
                print("append_tube")
            else:
                print("booolean")
                boolean_operation3 = vtk.vtkBooleanOperationPolyDataFilter()
                boolean_operation3.SetOperationToUnion()
                boolean_operation3.SetInputData(0, appended_data)
                boolean_operation3.SetInputData(1, tube)
                print("before update")
                boolean_operation3.Update()
                appended_data = boolean_operation3.GetOutput()

    # import ipdb; ipdb.set_trace()
    print("konec")

    # del (cylinderTri)
    # del (sphere1Tri)
    # del (sphere2Tri)
    # del (boolean_operation1)
    # del (boolean_operation2)
    # del (boolean_operation3)
    logger.debug("konec gen_tree()")
    # appendFilter.Update()
    # appended_data = appendFilter.GetOutput()
    return appended_data

def get_tube_old(radius, point, direction, length,
                 sphere_resolution, cylinder_resolution,
                 cylinder_radius_compensation_factor=1.0,
                 sphere_radius_compensation_factor=1.0,
                 tube_shape=True
                 ):
    """ Create a tube with ending half-spherese in Z axis.

    :param radius:
    :param point:
    :param direction:
    :param length:
    :param sphere_resolution:
    :param cylinder_resolution:
    :param cylinder_radius_compensation_factor:
    :param sphere_radius_compensation_factor:
    :param tube_shape:
    :return:
    """
    cylinderTri = vtk.vtkTriangleFilter()
    sphere1Tri = vtk.vtkTriangleFilter()
    sphere2Tri = vtk.vtkTriangleFilter()
    boolean_operation1 = vtk.vtkBooleanOperationPolyDataFilter()
    boolean_operation2 = vtk.vtkBooleanOperationPolyDataFilter()
    boolean_operation1.SetOperationToUnion()
    boolean_operation2.SetOperationToUnion()

    # print(dbg_msg)
    cylinder_radius = radius * cylinder_radius_compensation_factor
    sphere_radius = radius * sphere_radius_compensation_factor
    retval = None

    if tube_shape:
        sphere1 = get_sphere(point, sphere_radius, resolution=sphere_resolution)
        sphere1Tri.SetInputData(sphere1)
        sphere1Tri.Update()
    if length > 0:
        cylinder = get_cylinder(point,
                                length,
                                cylinder_radius,
                                direction,
                                resolution=cylinder_resolution)

        cylinderTri.SetInputData(cylinder)
        cylinderTri.Update()
        direction /= nm.linalg.norm(direction)

        lv = point + direction * length
        if tube_shape:
            sphere2 = get_sphere(lv, sphere_radius, resolution=sphere_resolution)
            sphere2Tri.SetInputData(sphere2)
            sphere2Tri.Update()

            # booleanOperation.SetInputData(0, cyl)
            boolean_operation1.SetInputData(0, cylinderTri.GetOutput())
            boolean_operation1.SetInputData(1, sphere1Tri.GetOutput())
            boolean_operation1.Update()
            boolean_operation2.SetInputData(0, boolean_operation1.GetOutput())
            boolean_operation2.SetInputData(1, sphere2Tri.GetOutput())
            # booleanOperation.SetInputData(2, sph2)
            boolean_operation2.Update()
        else:
            boolean_operation2 = cylinderTri

        retval = boolean_operation2.GetOutput()
    else:
        if tube_shape:
            boolean_operation2 = sphere1Tri
            retval = boolean_operation2.GetOutput()
        else:
            # length == 0 but no spheres
            # so we are generating just flat shape
            # boolean_operation2 = cylinderTri

            # return empty space
            retval = vtk.vtkTriangleFilter().GetOutput()
            # something_to_add = False
    return retval


def gen_tree_old(tree_data):
    import vtk
    points = vtk.vtkPoints()
    polyData = vtk.vtkPolyData()
    polyData.Allocate(1000, 1)
    polyData.SetPoints(points)
    poffset = 0

    for br in tree_data:
        logger.debug("generating edge " + str(br["length"]))
        cyl = get_cylinder(br['upperVertex'],
                           br['length'],
                           br['radius'],
                           br['direction'],
                           resolution=16)

        for ii in xrange(cyl.GetNumberOfPoints()):
            points.InsertPoint(poffset + ii, cyl.GetPoint(ii))

        for ii in xrange(cyl.GetNumberOfCells()):
            cell = cyl.GetCell(ii)
            cellIds = cell.GetPointIds()
            for jj in xrange(cellIds.GetNumberOfIds()):
                oldId = cellIds.GetId(jj)
                cellIds.SetId(jj, oldId + poffset)

            polyData.InsertNextCell(cell.GetCellType(),
                                    cell.GetPointIds())
        # spheres part
        # cyl = get_sphere(br['upperVertex'],
        #                    br['radius']
        #                    )
        # for ii in xrange(cyl.GetNumberOfPoints()):
        #     points.InsertPoint(poffset + ii, cyl.GetPoint(ii))
        #
        # for ii in xrange(cyl.GetNumberOfCells()):
        #     cell = cyl.GetCell(ii)
        #     cellIds = cell.GetPointIds()
        #     for jj in xrange(cellIds.GetNumberOfIds()):
        #         oldId = cellIds.GetId(jj)
        #         cellIds.SetId(jj, oldId + poffset)
        #
        #     polyData.InsertNextCell(cell.GetCellType(),
        #                             cell.GetPointIds())

        poffset += cyl.GetNumberOfPoints()

    return polyData


def compatibility_processing(indata):
    scale = 1e-3
    scale = 1

    outdata = []
    for key in indata:
        ii = indata[key]
        # logger.debug(ii)
        br = {}

        lengthEstimation = None
        try:
            # old version of yaml tree
            vA = ii['upperVertexXYZmm']
            vB = ii['lowerVertexXYZmm']
            radi = ii['radius']
            lengthEstimation = ii['length']
        except:
            # new version of yaml
            try:
                vA = ii['nodeA_ZYX_mm']
                vB = ii['nodeB_ZYX_mm']
                radi = ii['radius_mm']
                if "lengthEstimation" in ii.keys():
                    lengthEstimation = ii['lengthEstimation']
            except:
                import traceback
                logger.debug(traceback.format_exc())
                continue

        br['upperVertex'] = nm.array(vA) * scale
        br['radius'] = radi * scale
        if lengthEstimation is None:

            br['real_length'] = None
        else:
            br['real_length'] = lengthEstimation * scale

        vv = nm.array(vB) * scale - br['upperVertex']
        br['direction'] = vv / nm.linalg.norm(vv)
        br['length'] = nm.linalg.norm(vv)
        outdata.append(br)

    return outdata


def fix_tree_structure(tree_raw_data):
    if 'graph' in tree_raw_data:
        trees = tree_raw_data['graph']
    else:
        trees = tree_raw_data['Graph']
    return trees


def vt_file_2_vtk_file(infile, outfile, text_label=None):
    """
    From vessel_tree.yaml to output.vtk

    :param vessel_tree:  vt structure
    :param outfile: filename with .vtk extension
    :param text_label: text label like 'porta' or 'hepatic_veins'
    :return:

    """
    yaml_file = open(infile, 'r')
    tree_raw_data = yaml.load(yaml_file)
    vt2vtk_file(tree_raw_data, outfile, text_label)


def vt2vtk_file(vessel_tree, outfile, text_label=None):
    """
    vessel_tree structure
    :param vessel_tree:  vt structure
    :param outfile: filename with .vtk extension
    :param text_label: text label like 'porta' or 'hepatic_veins'
    :return:
    """
    import vtk
    trees = fix_tree_structure(vessel_tree)

    tkeys = trees.keys()
    if text_label is None:
        text_label = tkeys[0]

    tree_data = compatibility_processing(trees[text_label])
    polyData = gen_tree(tree_data)

    # import ipdb;
    # ipdb.set_trace()
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(outfile)
    try:
        writer.SetInputData(polyData)
    except:
        logger.warning("old vtk is used")
        writer.SetInput(polyData)
    writer.Write()


# From  geometry3.py

def circumscribed_polygon_radius(n, radius=1.0):
    """ Get circumscribed polygon radius.

    :param n: number of polygon elements
    :param radius: radius of inscribed circle
    :return: radius (distance from center to the corner) of polygon circumscribed to the
     circle
    """

    theta = 2 * np.pi / n
    radius_out = radius / np.cos(theta/2)

    return radius_out

def inscribed_polygon_radius(radius, n):
    """

    :param radius:
    :return:
    """
    pass

def regular_polygon_area_equivalent_radius(n, radius=1.0):
    """ Compute equivalent radius to obtain same surface area as circle.

    \theta = \frac{2 \pi}{n}

    r_{eqS} = \sqrt{\frac{\theta r^2}{\sin{\theta}}}

    simplier form:

    r_{eqS} = \sqrt{\frac{2 \pi r^2}{n \sin{\frac{2 \pi}{n}}}}

    :param radius: circle radius
    :param n:  number of regular polygon segments
    :return:  equivalent regular polynom surface
    """

    theta = 2 * np.pi / n

    r = np.sqrt((theta * radius**2) / np.sin(theta))
    return r

def regular_polygon_perimeter_equivalent_radius(n, radius=1.0):
    """ Compute equivalent radius to obtain same perimeter as circle.

    \theta = \frac{2 \pi}{n}

    r_{eqP} = \frac{\theta r}{2 \sin{\frac{\theta}{2}}}

    siplier form:
    r_{eqP} = \frac{\pi r}{n \sin{\frac{\pi}{n}}}

    :param radius: circle radius
    :param n:  number of regular polygon segments
    :return:  equivalent regular polynom surface
    """

    theta = 2 * np.pi / n

    r = (theta * radius) / (2 * np.sin(theta/2.0))
    return r

def main():
    logger = logging.getLogger()

    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    logger.addHandler(ch)

    # create file handler which logs even debug messages
    # fh = logging.FileHandler('log.txt')
    # fh.setLevel(logging.DEBUG)
    # formatter = logging.Formatter(
    #     '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # fh.setFormatter(formatter)
    # logger.addHandler(fh)
    # logger.debug('start')

    # input parser
    parser = argparse.ArgumentParser(
        description=__doc__
    )
    parser.add_argument(
        'inputfile',
        default=None,
        help='input file'
    )
    parser.add_argument(
        'outputfile',
        default='output.vtk',
        nargs='?',
        help='output file'
    )
    parser.add_argument(
        '-l', '--label',
        default=None,
        help='text label of vessel tree. f.e. "porta" or "hepatic_veins". \
        First label is used if it is set to None'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')
    args = parser.parse_args()

    if args.debug:
        ch.setLevel(logging.DEBUG)

    vt_file_2_vtk_file(args.inputfile, args.outputfile, args.label)


if __name__ == "__main__":
    main()
