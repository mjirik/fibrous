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

try:
    import skelet3d

    data3d = np.ones([3, 7, 9])
    data3d[:, 3, 3:6] = 0
    skelet3d.skelet3d(data3d)
    skelet3d_installed = True
    # skelet3d
except:
    skelet3d_installed = False
    logger.warning("skelet3d is not working")

try:
    import larcc

    larcc_installed = True
except:
    larcc_installed = False
    logger.warning("larcc is not working")

from fibrous.tree import TreeBuilder

from fibrous.tb_vtk import TBVTK
from fibrous.tb_volume import TBVolume
# There is some problem with VTK. Code seams to be fine but it fails
#  Generic Warning: In /tmp/vtk20150408-2435-1y7p97u/VTK-6.2.0/Common/Core/vtkObjectBase.cxx, line 93
#  Trying to delete object with non-zero reference count.
#  ERROR: In /tmp/vtk20150408-2435-1y7p97u/VTK-6.2.0/Common/Core/vtkObject.cxx, line 156
#  vtkObject (0x11a26e760): Trying to delete object with non-zero reference count.


VTK_MALLOC_PROBLEM = True


#

class TubeTreeTest(unittest.TestCase):
    def setUp(self):
        self.interactivetTest = False

    # interactivetTest = True

    @pytest.mark.LAR
    @unittest.skipIf(not ("larcc" in sys.modules), "larcc is not installed")
    def test_vessel_tree_lar(self):
        import fibrous.tb_lar
        tvg = TreeBuilder(fibrous.tb_lar.TBLar)
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa
        if self.interactiveTests:
            tvg.show()

    def test_import_all(self):
        import fibrous.image_manipulation
        import fibrous.tb_lar
        import fibrous.tb_lar_kunes
        import fibrous.tb_lar_smooth
        import fibrous.tb_volume
        import fibrous.tb_vtk
        import fibrous.tree

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_nothing(self):
        print("skelet3d_installed", skelet3d_installed)
        # import ipdb; ipdb.set_trace()
        self.assertTrue(False)

    def sample_tube_skeleton(self):
        sample_tube = {
            1:{
                "nodeA_ZYX_mm": [0, 30, 20],
                "nodeB_ZYX_mm": [0, 35, 25],
                "radius_mm": 3,
            },
            2:{
                "nodeA_ZYX_mm": [0, 30, 20],
                "nodeB_ZYX_mm": [15, 25, 20],
                "radius_mm": 3,
            },
            3:{
                "nodeA_ZYX_mm": [10, 10, 20],
                "nodeB_ZYX_mm": [5, 10, 20],
                "radius_mm": 3,
            },
            # TODO these are failing
            # 4:{
            #     "nodeA_ZYX_mm": [10, 10, 20],
            #     "nodeB_ZYX_mm": [5, 10, 21],
            #     "radius_mm": 3,
            # },
            # 5:{
            #     "nodeA_ZYX_mm": [10, 10, 20],
            #     "nodeB_ZYX_mm": [15, 14, 19],
            #     "radius_mm": 3,
            # },
            # 6:{
            #     "nodeA_ZYX_mm": [10, 10, 20],
            #     "nodeB_ZYX_mm": [-15, 10, -19],
            #     "radius_mm": 1,
            # },
        }
        return sample_tube

    def test_vessel_tree_vtk_with_new_subclass_on_artifical_sample_data(self):
        output_file = "test_output.vtk"
        if os.path.exists(output_file):
            os.remove(output_file)

        tube_skeleton = self.sample_tube_skeleton()
        tvg = TBVTK()
        tvg.set_model1d(model1d=tube_skeleton)
        # yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        # tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa
        tvg.saveToFile(output_file)

        self.assertTrue(os.path.exists(output_file))

    def test_vessel_tree_volume_with_new_subclass_on_artifical_sample_data(self):
        tube_skeleton = self.sample_tube_skeleton()
        tvg = TBVolume()
        tvg.set_model1d(model1d=tube_skeleton)
        # yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        # tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa

        self.assertTrue(type(output) == np.ndarray)

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_vessel_tree_vtk_with_new_subclass(self):
        from fibrous.tb_vtk import TBVTK
        tvg = TBVTK()
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa

    def test_vessel_tree_volume_with_new_subclass(self):
        from fibrous.tb_volume import TBVolume
        tvg = TBVolume()
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.init_data3d(shape=[100, 100, 100])
        output = tvg.buildTree()  # noqa
        # import sed3
        # ed = sed3.sed3(output)
        # ed.show()

        self.assertEqual(output[50, 20, 20], 200)
        self.assertEqual(output[5, 5, 5], 20)

    def test_get_tree_simple(self):
        from fibrous.tb_volume import TBVolume
        tvg = TBVolume()
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)

        # tvg.init_data3d(shape=[100, 100, 100])
        # output = tvg.buildTree()  # noqa
        from fibrous.tb_vtk import gen_tree_simple
        # from fibrous.tree import single_tree_compatibility_to_old
        vtk_tree = gen_tree_simple(tvg.tube_skeleton)
        self.assertTrue(vtk_tree is not None)

    def test_get_vt_from_file_and_save_it_to_vtk(self):
        import fibrous.tb_vtk
        fn="output.vtk"

        if os.path.exists(fn):
            os.remove(fn)

        # tvg = TBVolume()
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        fibrous.tb_vtk.vt_file_2_vtk_file(yaml_path, outfile=fn, tube_shape=False, use_simple_cylinder_method=True)

        self.assertTrue(os.path.exists(fn))

    @unittest.skip("This test is always failing from unknownd reason")
    def test_get_vt_from_file_and_save_it_to_vtk_with_not_simple_method(self):
        # TODO findout what is wrong with VTK export with boolean operations
        import fibrous.tb_vtk
        fn="output.vtk"

        if os.path.exists(fn):
            os.remove(fn)

        tvg = TBVolume()
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        fibrous.tb_vtk.vt_file_2_vtk_file(yaml_path, outfile=fn, tube_shape=False, use_simple_cylinder_method=False)


        self.assertTrue(os.path.exists(fn))

    def test_synthetic_volumetric_data_generation(self):
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
        self.assertEqual(data3d.shape[0], 100)
        self.assertEqual(data3d.shape[1], 101)
        self.assertEqual(data3d.shape[2], 102)

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_tube_vessel_tree_vtk_with_new_subclass_on_artifical_sample_data(self):
        tube_skeleton = self.sample_tube_skeleton()
        tvg = TBVTK(tube_shape=True)
        tvg.set_model1d(model1d=tube_skeleton)
        # yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        # tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa
        output_file = "test_output.vtk"
        tvg.saveToFile(output_file)

        self.assertTrue(os.path.exists(output_file))



    @unittest.skip("test debug")
    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_vessel_tree_vtk(self):
        tvg = TreeBuilder('vtk')
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa
        # tvg.show()
        # tvg.saveToFile("tree_output.vtk")

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    @unittest.skipIf(not ("skelet3d" in sys.modules), "skelet3d is not installed")
    @unittest.skipIf(not skelet3d_installed, "skelet3d is not installed")
    def test_vessel_tree_vtk_from_skeleton_failing(self):
        print("skelet3d_installed", skelet3d_installed)

        import skelet3d
        import skelet3d.skeleton_analyser
        import shutil

        fn_out = 'tree.vtk'
        if os.path.exists(fn_out):
            os.remove(fn_out)

        volume_data = np.zeros([3, 7, 9], dtype=np.int)
        volume_data[:, :, 1:3] = 1
        volume_data[:, 5, 2:9] = 1
        volume_data[:, 0:7, 5] = 1
        skelet = skelet3d.skelet3d(volume_data)

        skan = skelet3d.skeleton_analyser.SkeletonAnalyser(skelet, volume_data=volume_data, voxelsize_mm=[1, 1, 1])
        stats = skan.skeleton_analysis()

        # tvg = TreeBuilder('vtk')
        tvg = TBVTK()
        tvg.set_model1d(stats)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        tvg.tree_data = stats
        output = tvg.buildTree()  # noqa
        tvg.saveToFile(fn_out)
        os.path.exists(fn_out)

    # @unittest.skipIf(not skelet3d_installed, "skelet3d is not installed")
    def test_vessel_tree_vtk_from_skeleton_of_one_tube(self):
        print("skelet3d_installed", skelet3d_installed)

        import skelet3d
        import skelet3d.skeleton_analyser
        import shutil

        fn_out = 'tree_one_tube.vtk'
        if os.path.exists(fn_out):
            os.remove(fn_out)

        volume_data = np.zeros([7, 8, 9], dtype=np.int)
        volume_data[4:8, 4:6, 1:3] = 1
        volume_data[:, 5, 2:9] = 1
        volume_data[:, 0:7, 5] = 1
        skelet = skelet3d.skelet3d(volume_data)

        skan = skelet3d.skeleton_analyser.SkeletonAnalyser(skelet, volume_data=volume_data, voxelsize_mm=[1, 1, 1])
        stats = skan.skeleton_analysis()

        self.assertEqual(len(stats), 1, "There should be just one cylinder based on data with different diameter")
        # tvg = TreeBuilder('vtk')
        tvg = TBVTK()
        tvg.set_model1d(stats)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        tvg.tree_data = stats
        output = tvg.buildTree()  # noqa
        tvg.saveToFile(fn_out)
        self.assertTrue(os.path.exists(fn_out))

    @unittest.skipIf(not skelet3d_installed, "skelet3d is not installed")
    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_vessel_tree_vtk_from_skeleton_with_more_tubes(self):
        """Test dont fail but there is no right output in vtk file"""
        print("skelet3d_installed", skelet3d_installed)

        import skelet3d
        import skelet3d.skeleton_analyser
        import shutil

        fn_out = 'tree_more_tubes.vtk'
        if os.path.exists(fn_out):
            os.remove(fn_out)

        volume_data = np.zeros([20, 21, 22], dtype=np.int8)
        # croess
        volume_data[8:11, 14:17, 4:14] = 1
        volume_data[9:15, 11:15, 9:19] = 1
        volume_data[8:12, 5:9, 4:14] = 1
        volume_data[9:15, 2:7, 9:19] = 1
        # volume_data[10:12, 2:15, 13] = 1
        # volume_data[11:13, 5:12, 14] = 1
        # volume_data[12:14, 16, 13] = 1
        skelet = skelet3d.skelet3d(volume_data)

        skan = skelet3d.skeleton_analyser.SkeletonAnalyser(
            skelet, volume_data=volume_data, voxelsize_mm=[1, 1, 1],
            cut_wrong_skeleton=False
        )
        stats = skan.skeleton_analysis()

        # self.assertEqual(len(stats), 1, "There should be just one cylinder based on data with different diameter")
        # tvg = TreeBuilder('vtk')
        tvg = TBVTK()
        tvg.set_model1d(stats)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        tvg.tree_data = stats
        output = tvg.buildTree()  # noqa
        tvg.saveToFile(fn_out)
        self.assertTrue(os.path.exists(fn_out))



    # TODO finish this test
    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_vessel_tree_vol(self):
        import fibrous.tb_volume
        tvg = TreeBuilder(fibrous.tb_volume.TBVolume)
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        output = tvg.buildTree()  # noqa
        # tvg.show()
        # if self.interactiveTests:
        #     tvg.show()

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_import_new_vt_format(self):

        tvg = TreeBuilder()
        yaml_path = os.path.join(path_to_script, "vt_biodur_simple.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [150, 150, 150]
        data3d = tvg.buildTree()

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_cylinders_generator(self):
        from fibrous.generators.cylinders import CylinderGenerator

        cg = CylinderGenerator()
        cg.run()

    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_vtk_tree(self):
        import numpy as np
        tree_data = {

        }
        element_number = 1
        area_size = 100
        radius = 5
        np.random.seed(0)
        pts1 = np.random.random([element_number, 3]) * (area_size - 4 * radius) + 2 * radius
        pts2 = np.random.random([element_number, 3]) * (area_size - 4 * radius) + 2 * radius
        for i in range(element_number):
            edge = {
                # "nodeA_ZYX_mm": vor3.vertices[simplex],
                # "nodeB_ZYX_mm": vor3.vertices[simplex],
                "nodeA_ZYX_mm": pts1[i],
                "nodeB_ZYX_mm": pts2[i],
                "radius_mm": radius
            }
            tree_data[i] = edge

        tvg = TreeBuilder('vtk')
        # yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        # tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [area_size, area_size, area_size]
        tvg.tree_data = tree_data
        output = tvg.buildTree()  # noqa
        # tvg.show()
        tvg.saveToFile("test_tree_output.vtk")

    # @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    @unittest.skipIf(VTK_MALLOC_PROBLEM, "VTK malloc problem")
    def test_tree_generator(self):
        import numpy as np
        tree_data = {

        }
        element_number = 6
        np.random.seed(0)
        pts = np.random.random([element_number, 3]) * 100

        # construct voronoi
        import scipy.spatial
        import itertools
        vor3 = scipy.spatial.Voronoi(pts)

        # for i, two_points in enumerate(vor3.ridge_points):
        for i, simplex in enumerate(vor3.ridge_vertices):
            simplex = np.asarray(simplex)
            # fallowing line removes all ridges with oulayers
            simplex = simplex[simplex > 0]
            if np.all(simplex >= 0):

                x = vor3.vertices[simplex, 0]
                y = vor3.vertices[simplex, 1]
                z = vor3.vertices[simplex, 2]
                for two_points in itertools.combinations(simplex, 2):
                    edge = {
                        # "nodeA_ZYX_mm": vor3.vertices[simplex],
                        # "nodeB_ZYX_mm": vor3.vertices[simplex],
                        "nodeA_ZYX_mm": vor3.vertices[two_points[0]],
                        "nodeB_ZYX_mm": vor3.vertices[two_points[1]],
                        "radius_mm": 2
                    }
                    tree_data[i] = edge
            else:
                pass

        show_input_points = False
        if show_input_points:
            length = len(tree_data)
            for i in range(element_number):
                edge = {
                    #         #"nodeA_ZYX_mm": np.random.random(3) * 100,
                    "nodeA_ZYX_mm": pts[i - 1],
                    "nodeB_ZYX_mm": pts[i],
                    #         "nodeB_ZYX_mm": np.random.random(3) * 100,
                    "radius_mm": 1
                }
                tree_data[i + length] = edge

        tvg = TreeBuilder('vtk')
        yaml_path = os.path.join(path_to_script, "./hist_stats_test.yaml")
        # tvg.importFromYaml(yaml_path)
        tvg.voxelsize_mm = [1, 1, 1]
        tvg.shape = [100, 100, 100]
        tvg.tree_data = tree_data
        output = tvg.buildTree()  # noqa
        # tvg.show()
        tvg.saveToFile("test_tree_output.vtk")

        tvgvol = TreeBuilder('vol')
        tvgvol.voxelsize_mm = [1, 1, 1]
        tvgvol.shape = [100, 100, 100]
        tvgvol.tree_data = tree_data
        outputvol = tvgvol.buildTree()
        tvgvol.saveToFile("tree_volume.pklz")
        # self.assertTrue(False)

    def test_import_new_vt_format(self):
        from fibrous.tree import TreeBuilder

        # tvg = TreeBuilder()
        tvg = TBVolume()

        yaml_path = os.path.join(path_to_script, "vt_biodur_simple.yaml")
        tvg.importFromYaml(yaml_path)
        tvg.set_area_sampling(voxelsize_mm=[1,1,1], shape=[150, 150, 150])
        # tvg.voxelsize_mm = [1, 1, 1]
        # tvg.shape = [150, 150, 150]
        data3d = tvg.buildTree()

    @unittest.skipUnless(os.path.exists(r"e:/vessel_tree_data/ep_hcc2_porta1d.yaml"), "data are not available")
    def test_read_portal_vein_data(self):
        import fibrous.tree
        from fibrous.tree import TreeBuilder

        # tvg = TreeBuilder()
        tvg = TBVolume()

        # yaml_path = os.path.join(path_to_script, "vt_biodur_simple.yaml")
        yaml_path = r"e:/vessel_tree_data/ep_hcc2_porta1d.yaml"
        tree = fibrous.tree.read_tube_skeleton_from_yaml(yaml_path)
        # tvg.importFromYaml(yaml_path)
        # tvg.set_area_sampling(voxelsize_mm=[1,1,1], shape=[150, 150, 150])
        # tvg.voxelsize_mm = [1, 1, 1]
        # tvg.shape = [150, 150, 150]
        # data3d = tvg.buildTree()

if __name__ == "__main__":
    unittest.main()
