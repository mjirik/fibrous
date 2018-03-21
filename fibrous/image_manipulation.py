#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def get_all_area_sampling_parameters(
        voxelsize_mm=None,
        areasize_mm=None,
        areasize_px=None,
        data3d=None,
):
    """
    Calculate all area sampling dimensions
    :param voxelsize_mm:
    :param areasize_mm:
    :param areasize_px:
    :param data3d:
    :return: voxelsize_mm, areasize_mm, areasize_px
    """

    if data3d is not None:
        areasize_px = data3d.shape

    if areasize_mm is None and voxelsize_mm is None:
        voxelsize_mm = [1, 1, 1]
    if areasize_px is None and voxelsize_mm is None:
        voxelsize_mm = [1, 1, 1]
    if voxelsize_mm is None:
        voxelsize_mm = (np.asarray(areasize_mm) / areasize_px).astype(np.float)

    voxelsize_mm = np.asarray(voxelsize_mm, dtype=np.float)

    if areasize_mm is None:
        areasize_mm = np.asarray(areasize_px) * voxelsize_mm

    if areasize_px is None:
        areasize_px = (np.asarray(areasize_mm) / voxelsize_mm).astype(np.int)

    return voxelsize_mm, areasize_mm, areasize_px
