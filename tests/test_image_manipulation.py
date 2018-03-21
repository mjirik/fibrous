#! /usr/bin/python
# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)
import unittest
import numpy as np
from fibrous import image_manipulation as imma



class ImageManipulationTest(unittest.TestCase):
    def assertListsEqual(self, lst1, lst2):

        lst1 = np.asarray(lst1)
        lst2 = np.asarray(lst2)

        err = np.sum(np.abs(lst1 - lst2))
        self.assertEqual(err, 0)

    def test_area_sampling_no_voxelsize(self):
        voxelsize_mm, areasize_mm, shape = imma.get_all_area_sampling_parameters(
            areasize_mm=[10, 15, 12]
        )
        self.assertListsEqual(voxelsize_mm, [1, 1, 1])
        self.assertListsEqual(shape, [10, 15, 12])

    def test_area_sampling_from_area_mm(self):
        voxelsize_mm, areasize_mm, shape = imma.get_all_area_sampling_parameters(
            areasize_mm=[10, 15, 12], voxelsize_mm=[2, 3, 1]
        )
        self.assertListsEqual(shape, [5, 5, 12])

    def test_area_sampling_from_area_px(self):
        voxelsize_mm, areasize_mm, shape = imma.get_all_area_sampling_parameters(
            areasize_px=[5, 5, 12], voxelsize_mm=[2, 3, 1]
        )
        self.assertListsEqual(shape, [5, 5, 12])
        self.assertListsEqual(areasize_mm, [10, 15, 12])

    def test_area_sampling_from_data3d(self):
        voxelsize_mm, areasize_mm, shape = imma.get_all_area_sampling_parameters(
            data3d=np.zeros([5, 5, 12]), voxelsize_mm=[2, 3, 1]
        )
        self.assertListsEqual(shape, [5, 5, 12])
        self.assertListsEqual(areasize_mm, [10, 15, 12])

    def test_area_sampling_calculate_voxelsize(self):
        voxelsize_mm, areasize_mm, shape = imma.get_all_area_sampling_parameters(
            areasize_mm=[10, 15, 12], areasize_px=[5, 5, 12]
        )
        self.assertListsEqual(voxelsize_mm, [2, 3, 1])

