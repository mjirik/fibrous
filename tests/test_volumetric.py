#! /usr/bin/python
# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger(__name__)
# import funkcí z jiného adresáře
import os
import os.path

import pytest

path_to_script = os.path.dirname(os.path.abspath(__file__))
import unittest
import numpy as np
import sys
import io3d
from matplotlib import pyplot as plt

def test_synthetic_volumetric_data_generation():
    """
    Generovani umeleho stromu do 3D a jeho evaluace.
    V testu dochazi ke kontrole predpokladaneho objemu a deky cev


    """
    # from skelet3d.tree import TreeGenerator
    # from fibrous.tree import TubeSkeletonBuilder as TreeGenerator
    # from fibrous.tree import TreeBuilder as TreeGenerator
    from fibrous.tb_volume import TBVolume as TreeGenerator

    # print("zacatek podezreleho testu")
    # import segmentation
    # import misc

    # generate 3d data from yaml for testing
    tvg = TreeGenerator()
    yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
    tvg.importFromYaml(yaml_path)
    tvg.set_area_sampling(voxelsize_mm=[1, 1, 1], shape=[100, 101, 102])
    # tvg.voxelsize_mm = [1, 1, 1]
    # tvg.shape = [100, 100, 100]
    # data3d = tvg.generateTree()
    data3d = tvg.buildTree()
    assert data3d.shape[0] == 100
    assert data3d.shape[1] == 101
    assert data3d.shape[2] == 102

    # plt.imshow(data3d[:,:,50])
    # plt.colorbar()
    # plt.show()
    io3d.write(data3d.astype(np.uint8), path="./slice{:04d}.jpg")
    # self.assertEqual(data3d.shape[0], 100)
    # self.assertEqual(data3d.shape[1], 101)
    # self.assertEqual(data3d.shape[2], 102)

    tvg = TreeGenerator()
    yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
    tvg.importFromYaml(yaml_path)
    tvg.set_area_sampling(voxelsize_mm=[0.5, 0.5, 0.5], shape=[100, 101, 102])
    io3d.write(data3d.astype(np.uint8), path="./ssslice{:04d}.jpg")
