#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
Generator of histology report

"""
import logging

logger = logging.getLogger(__name__)

import sys
import os.path

path_to_script = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(path_to_script, "../extern/dicom2fem/src"))

import argparse
import numpy as np
import scipy.ndimage

if sys.version_info.major == 3:
    xrange = range

from . import tree
# import datareader

# import tb_vtk

# import vtk
# from vtk.util import numpy_support

# from datetime import datetime
# import collections
UNDERDEBUG = 8


class TBVolume(tree.TubeSkeletonBuilder):
    """
    This generator is called by generateTree() function as a general form.
    Other similar generator is used for generating LAR outputs.
    """

    def __init__(self,
                 **kwargs
                 ):

        # super(tree.FiberSkeletBuilder, self).__init__()
        tree.TubeSkeletonBuilder.__init__(self)
        self.init(**kwargs)

    def init_data3d(self, shape, background_intensity=20, dtype=np.int):
        self.shape = np.asarray(shape, dtype=np.int)
        self.data3d = (np.ones(shape, dtype=dtype) * background_intensity).astype(dtype=dtype)

    def init(self, tube_skeleton=None, shape=None, voxelsize_mm=None,
             background_intensity=20, dtype=np.int, intensity_profile=None):

        self.tube_skeleton = tube_skeleton
        if shape is None:
            shape = [100, 100, 100]
        if voxelsize_mm is None:
            voxelsize_mm = [1., 1., 1.]

        self.init_data3d(shape, background_intensity, dtype=dtype)
        self.voxelsize_mm = voxelsize_mm
        if intensity_profile is not None:
            self.intensity_profile = intensity_profile
        else:
            # self.intensity_profile = {1:200, 0.6: 100}
            self.intensity_profile = {1: 200}

        # self.intensity_profile = incollections.OrderedDict(sorted(intensity_profile, reverse=True))
        self._cylinders_params = []
        self._temp_intensity = 10
        self.finish_progress_callback = None
        # self.output_intensity = 200

    def add_cylinder(self, p1m, p2m, rad, id):
        """
        Funkce na vykresleni jednoho segmentu do 3D dat
        """
        self._cylinders_params.append([p1m, p2m, rad, id])

    def _add_cylinder(self, p1m, p2m, rad, id):

        cyl_data3d = np.ones(self.shape, dtype=np.bool)
        # prvni a koncovy bod, ve pixelech
        p1 = [p1m[0] / self.voxelsize_mm[0], p1m[1] /
              self.voxelsize_mm[1], p1m[2] / self.voxelsize_mm[2]]
        p2 = [p2m[0] / self.voxelsize_mm[0], p2m[1] /
              self.voxelsize_mm[1], p2m[2] / self.voxelsize_mm[2]]
        logger.log(UNDERDEBUG,
                   "p1_px: " + str(p1[0]) + " " + str(p1[1]) + " " + str(p1[2]))
        logger.log(UNDERDEBUG,
                   "p2_px: " + str(p2[0]) + " " + str(p2[1]) + " " + str(p2[2]))
        logger.log(UNDERDEBUG, "radius_mm:" + str(rad))

        # vzdalenosti mezi prvnim a koncovim bodem (pro jednotlive osy)
        pdiff = [abs(p1[0] - p2[0]), abs(p1[1] - p2[1]), abs(p1[2] - p2[2])]

        # generovani hodnot pro osu segmentu
        num_points = max(pdiff) * \
                     2  # na jeden "pixel nejdelsi osy" je 2 bodu primky (shannon)
        zvalues = np.linspace(p1[0], p2[0], num_points)
        yvalues = np.linspace(p1[1], p2[1], num_points)
        xvalues = np.linspace(p1[2], p2[2], num_points)

        # drawing a line
        no_index_error_occured = True
        for i in range(0, len(xvalues)):
            # TODO make something with indexes out of requested area
            try:
                cyl_data3d[int(zvalues[i])][int(yvalues[i])][int(xvalues[i])] = 0
            except IndexError:
                if no_index_error_occured:
                    import traceback
                    traceback.print_exc()
                    logger.warning("Cylinder drawing out of bounds. Other same type warnings are suppressed.")
                    no_index_error_occured = False
            except:
                import traceback
                traceback.print_exc()
                print("except in drawing line")
                logger.warning("Cylinder drawing problem")
                # import ipdb; ipdb.set_trace() #  noqa BREAKPOINT

        # cuting size of 3d space needed for calculating distances (smaller ==
        # a lot faster)
        cut_up = max(
            0, round(min(p1[0], p2[0]) - (rad / min(self.voxelsize_mm)) - 2))
        # ta 2 je kuli tomu abyh omylem nurizl
        cut_down = min(self.shape[0], round(
            max(p1[0], p2[0]) + (rad / min(self.voxelsize_mm)) + 2))
        cut_yu = max(
            0, round(min(p1[1], p2[1]) - (rad / min(self.voxelsize_mm)) - 2))
        cut_yd = min(self.shape[1], round(
            max(p1[1], p2[1]) + (rad / min(self.voxelsize_mm)) + 2))
        cut_xl = max(
            0, round(min(p1[2], p2[2]) - (rad / min(self.voxelsize_mm)) - 2))
        cut_xr = min(self.shape[2], round(
            max(p1[2], p2[2]) + (rad / min(self.voxelsize_mm)) + 2))
        logger.log(UNDERDEBUG, "cutter_px: z_up-" + str(cut_up) + " z_down-" + str(cut_down) + " y_up-" + str(
            cut_yu) + " y_down-" + str(cut_yd) + " x_left-" + str(cut_xl) + " x_right-" + str(cut_xr))
        cyl_data3d_cut = cyl_data3d[
                         int(cut_up):int(cut_down),
                         int(cut_yu):int(cut_yd),
                         int(cut_xl):int(cut_xr)]

        # calculating distances
        # spotrebovava naprostou vetsinu casu (pro 200^3  je to kolem 1.2
        # sekundy, proto jsou data osekana)
        lineDst = scipy.ndimage.distance_transform_edt(
            cyl_data3d_cut, self.voxelsize_mm)

        # zkopirovani vyrezu zpet do celeho rozsahu dat
        for z in xrange(0, len(cyl_data3d_cut)):
            for y in xrange(0, len(cyl_data3d_cut[z])):
                for x in xrange(0, len(cyl_data3d_cut[z][y])):
                    if lineDst[z][y][x] <= rad:
                        iX = int(z + cut_up)
                        iY = int(y + cut_yu)
                        iZ = int(x + cut_xl)
                        self.data3d[iX][iY][iZ] = self._temp_intensity

    def get_output(self):
        return self.data3d

    def finish(self):
        """
        :param self.finish_progress_callback(self, progress): function with iprogress parameter from 0.0 to 1.0
        :return:
        """
        progress_step = 1.0 / (len(self.intensity_profile) * len(self._cylinders_params))
        progress = 0.0

        for radk in sorted(self.intensity_profile, reverse=True):
            radk_intensity = self.intensity_profile[radk]

            for cyl in self._cylinders_params:
                self._add_cylinder(cyl[0], cyl[1], cyl[2] * radk, cyl[3])

                if self.finish_progress_callback is not None:
                    self.finish_progress_callback(self, progress)
                    progress += progress_step

            self.data3d[self.data3d == self._temp_intensity] = radk_intensity

        if self.finish_progress_callback is not None:
            self.finish_progress_callback(self, 1.0)

    def save(self, outputfile, filetype='pklz'):
        import io3d
        import io3d.misc
        import numpy as np
        data = {
            'data3d': self.data3d.astype(np.uint8),  # * self.output_intensity,
            'voxelsize_mm': self.voxelsize_mm,
            # 'segmentation': np.zeros_like(self.data3d, dtype=np.int8)
        }

        # data3d = np.zeros([10,10,10])
        # segmentation = np.zeros([10,10,10])
        #
        # data3d [2:7,:3:5, :6] = 100
        # datap = {
        #     "data3d": data3d,
        #     "segmentation": segmentation,
        #     "voxelsize_mm": [1,1,1]
        # }
        # io3d.write(datap, "file1.pklz")
        # import ipdb; ipdb.set_trace()

        io3d.write(data, outputfile)
        # io3d.misc.obj_to_file(data, outputfile, filetype=filetype)
        # dw = datawriter.DataWriter()
        # dw.Write3DData(self.data3d, outputfile, filetype)

    def show(self):
        import sed3 as se
        pyed = se.sed3(self.data3d)
        pyed.show()
