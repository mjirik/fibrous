#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR% %USER% <%MAIL%>
#
# Distributed under terms of the %LICENSE% license.

"""
%HERE%
"""

import logging

logger = logging.getLogger(__name__)
import argparse
import datetime
import numpy as np


class TubeSkeletonBuilder:
    def __init__(self):
        """
        This function can be used as vessel_tree iterator. Just implement generator_class

        :param generator_class: class with function add_cylinder(p1pix, p2pix, rad_mm) and get_output()
        :param generator_params:
        """
        self.rawdata = None
        self.tube_skeleton = None
        self.data3d = None
        self.voxelsize_mm = [1, 1, 1]
        self.shape = None
        self.use_lar = False
        self.tree_label = None
        self.segments_progress_callback = self.finish_progress_callback
        self.stop_processing = False

    def set_model1d(self, model1d, label=None):
        """
        Set the 1D model and make compatibility fixtures.

        :param model1d:
        :param label:
        :return:
        """
        self.tube_skeleton, self.rawdata = pick_model1d(rawdata=model1d, label=label)

    def importFromYaml(self, filename):
        tube_skeleton, rawdata = read_tube_skeleton_from_yaml(
            filename=filename,
            tree_label=self.tree_label,
            return_rawdata=True
        )
        self.rawdata = rawdata
        self.tube_skeleton = tube_skeleton

    def add_segment_to_tree(self, pointA, pointB, radius, id=None):
        """
        Before generation this can be used to add new segment
        :return:
        """
        if self.tube_skeleton is None:
            self.tube_skeleton = {}

        if id is None:
            id = len(self.tube_skeleton)

        self.tube_skeleton[id] = {
            'nodeA_ZYX_mm': pointA,
            'nodeB_ZYX_mm': pointB,
            'radius_mm': radius
        }

    def buildTree(self):
        """
        | Funkce na vygenerování objemu stromu ze zadaných dat.
        | Generates output by defined generator. If VolumeTreeGenerator is used, output is data3d.
        """
        # LAR init
        if self.use_lar:
            import lar_vessels
            self.lv = lar_vessels.LarVessels()

        # use generator init
        # if self.generator_params is None:
        #     self.generator = self.generator_class(self)
        # else:
        #     self.generator = self.generator_class(self, **self.generator_params)

        tree_output = self._build_tree_per_segments()

        logger.debug("before visualization - generateTree()")
        if self.use_lar:
            self.lv.show()

        return tree_output

    def _build_tree_per_segments(self):
        ln = len(self.tube_skeleton)
        if ln == 0:
            ln = 1
        progress_step = 1.0 / ln
        progress = 0.0

        for cyl_id in self.tube_skeleton:
            if self.stop_processing:
                break
            logger.debug("CylinderId: " + str(cyl_id))
            cyl_data = self.tube_skeleton[cyl_id]

            # try:
            #     cyl_data = self.data['graph']['porta'][cyl_id]
            # except:
            #     cyl_data = self.data['Graph'][cyl_id]

            # prvni a koncovy bod, v mm + radius v mm
            try:
                p1m = cyl_data['nodeA_ZYX_mm']  # souradnice ulozeny [Z,Y,X]
                p2m = cyl_data['nodeB_ZYX_mm']
                rad = cyl_data['radius_mm']
                self.add_cylinder(p1m, p2m, rad, cyl_id)
            except Exception as e:
                # import ipdb; ipdb.set_trace() #  noqa BREAKPOINT

                logger.error(
                    "Segment id " + str(cyl_id) + ": error reading data from yaml!: " + str(e))
                # return

                # if self.use_lar:
                #     self.generator.add_cylinder(p1m, p2m, rad, in)
            if self.segments_progress_callback is not None:
                self.segments_progress_callback(progress)
            progress += progress_step
        logger.debug("cylinders generated")

        # import ipdb; ipdb.set_trace()
        if "finish" in dir(self):
            # generator could have finish() function
            self.finish_progress_callback = self.finish_progress_callback
            self.finish()
            logger.debug("joints generated")
        else:
            import traceback
            # logger.debug(traceback.format_exc())
            logger.debug("no finish() function in tree constructor")
        # import ipdb; ipdb.set_trace()

        output = self.get_output()

        return output

    def finish_progress_callback(self, progress, *args, **kwargs):
        # if self.segments_progress_callback is not None:
        #     self.segments_progress_callback(progress)
        print("progess: {}\r".format(str(progress)))
        logger.debug(str(progress))

    def saveToFile(self, *args, **kwargs):
        self.save(*args, **kwargs)


    def stop(self):
        self.stop_processing = True

class TreeBuilder:
    def __init__(self, generator_class='volume', generator_params=None):
        """
        This function can be used as vessel_tree iterator. Just implement generator_class

        :param generator_class: class with function add_cylinder(p1pix, p2pix, rad_mm) and get_output()
        :param generator_params:
        """
        self.rawdata = None
        self.tree_data = None
        self.data3d = None
        self.voxelsize_mm = [1, 1, 1]
        self.shape = None
        self.use_lar = False
        self.tree_label = None
        self.segments_progress_callback = self.finish_progress_callback
        self.stop_processing = False
        logger.warning("TreeBuilder is deprecated. Use TubeSkeletonBuilder instead.")

        if generator_class in ['vol', 'volume']:
            from . import tb_volume
            generator_class = tb_volume.TBVolume
        elif generator_class in ['lar']:
            from . import tb_lar
            generator_class = tb_lar.TBLar
        elif generator_class in ['vtk']:
            from . import tb_vtk
            generator_class = tb_vtk.TBVTK
        elif generator_class in ['kunes']:
            from . import tb_lar_kunes
            generator_class = tb_lar_kunes.TBLar
        elif generator_class in ['larsm']:
            from . import tb_lar_smooth
            generator_class = tb_lar_smooth.TBLarSmooth
        elif generator_class in ['lar_nojoints']:
            from . import tb_lar
            generator_class = tb_lar.TBLar
            generator_params = {
                'endDistMultiplicator': 0,
                'use_joints': False
            }
        self.generator_class = generator_class
        self.generator_params = generator_params

    def fix_tree_structure(self, tree_raw_data):
        return backward_compatibility_tree_structure(tree_raw_data)

    def importFromYaml(self, filename):
        import yaml
        f = open(filename, 'rb')
        rawdata = yaml.load(f)
        f.close()
        self.rawdata = self.fix_tree_structure(rawdata)

        tkeys = list(self.rawdata['Graph'])
        if (self.tree_label is None) or (self.tree_label not in tkeys):
            self.tree_label = tkeys[0]
        tree_data = self.rawdata['Graph'][self.tree_label]
        self.tree_data = tree_data

    def add_segment_to_tree(self, pointA, pointB, radius, id=None):
        """
        Before generation this can be used to add new segment
        :return:
        """
        if self.tree_data is None:
            self.tree_data = {}

        if id is None:
            id = len(self.tree_data)

        self.tree_data[id] = {
            'nodeA_ZYX_mm': pointA,
            'nodeB_ZYX_mm': pointB,
            'radius_mm': radius
        }

    def buildTree(self):
        """
        | Funkce na vygenerování objemu stromu ze zadaných dat.
        | Generates output by defined generator. If VolumeTreeGenerator is used, output is data3d.
        """
        # LAR init
        if self.use_lar:
            import lar_vessels
            self.lv = lar_vessels.LarVessels()

        # use generator init
        if self.generator_params is None:
            self.generator = self.generator_class(self)
        else:
            self.generator = self.generator_class(self, **self.generator_params)

        tree_output = self._build_tree_per_segments()

        logger.debug("before visualization - generateTree()")
        if self.use_lar:
            self.lv.show()

        return tree_output

    def add_cylinder(p1m, p2m, rad, cyl_id):
        # It is mandatory implement this function
        pass

    def _build_tree_per_segments(self):
        ln = len(self.tree_data)
        if ln == 0:
            ln = 1
        progress_step = 1.0 / ln
        progress = 0.0

        for cyl_id in self.tree_data:
            if self.stop_processing:
                break
            logger.debug("CylinderId: " + str(cyl_id))
            cyl_data = self.tree_data[cyl_id]

            # try:
            #     cyl_data = self.data['graph']['porta'][cyl_id]
            # except:
            #     cyl_data = self.data['Graph'][cyl_id]

            # prvni a koncovy bod, v mm + radius v mm
            try:
                p1m = cyl_data['nodeA_ZYX_mm']  # souradnice ulozeny [Z,Y,X]
                p2m = cyl_data['nodeB_ZYX_mm']
                rad = cyl_data['radius_mm']
                self.generator.add_cylinder(p1m, p2m, rad, cyl_id)
            except Exception as e:
                # import ipdb; ipdb.set_trace() #  noqa BREAKPOINT

                logger.error(
                    "Segment id " + str(cyl_id) + ": error reading data from yaml!: " + str(e))
                # return

                # if self.use_lar:
                #     self.generator.add_cylinder(p1m, p2m, rad, in)
            if self.segments_progress_callback is not None:
                self.segments_progress_callback(progress)
            progress += progress_step
        logger.debug("cylinders generated")

        # import ipdb; ipdb.set_trace()
        if "finish" in dir(self.generator):
            # generator could have finish() function
            self.generator.finish_progress_callback = self.finish_progress_callback
            self.generator.finish()
            logger.debug("joints generated")
        else:
            import traceback
            # logger.debug(traceback.format_exc())
            logger.debug("no finish() function in tree constructor")
        # import ipdb; ipdb.set_trace()

        output = self.generator.get_output()

        return output

    def finish_progress_callback(self, progress, *args, **kwargs):
        # if self.segments_progress_callback is not None:
        #     self.segments_progress_callback(progress)
        print("progress: {}\r".format(str(progress)))
        logger.debug(str(progress))

    def saveToFile(self, *args, **kwargs):
        self.generator.save(*args, **kwargs)

    def show(self):
        self.generator.show()

    def stop(self):
        self.stop_processing = True

def parse_area_properties(rawdata):
    def find_in_general_key(general):
        area = {}
        if "voxelsize_mm" in general:
            area["voxelsize_mm"] = general["voxelsize_mm"]
        if "shape_px" in general:
            area["areasize_px"] = general["shape_px"]
        if "areasize_px" in general:
            area["areasize_px"] = general["areasize_px"]
        return area

    area = {}
    if "general" in rawdata.keys():
        general = rawdata["general"]
        area = find_in_general_key(general)
    if "areasampling" in rawdata.keys():
        general = rawdata["areasampling"]
        area = find_in_general_key(general)

    # compatibility with older paper data
    if "voxelsize_mm" in rawdata.keys():
        area["voxelsize_mm"] = rawdata["voxelsize_mm"]
    if "voxelsize_px" in rawdata.keys():
        area["areasize_px"] = rawdata["voxelsize_px"]

    area["areasize_mm"] = np.asarray(area["areasize_px"]) * np.asarray(area["voxelsize_mm"])
    return area

def read_tube_skeleton_from_yaml(filename, tree_label=None, return_rawdata=False):
    """ Get tube skeleton and raw data from yaml file.

    :param self:
    :param filename: yaml filename
    :param tree_label: label of tree. The first tree is used if None.
    :return:
    """
    import yaml
    f = open(filename, 'rb')
    rawdata = yaml.load(f)
    f.close()
    tube_skeleton, rawdataf = pick_model1d(rawdata, label=tree_label)

    if return_rawdata:
        return tube_skeleton, rawdataf
    else:
        return tube_skeleton

def pick_model1d(rawdata, label=None):
    """
    Fix all model problems and pick one of tree based on label.

    :param rawdata: structure with 1D model of fibrous material
    :param tree_label: select this tree or use the first one.
    :return: one submodel(based on label) of fibrous 1d model
    """


    rawdataf = backward_compatibility_tree_structure(rawdata)

    tkeys = list(rawdataf['Graph'])
    if (label is None) or (label not in tkeys):
        label = tkeys[0]
    tube_skeleton = rawdataf['Graph'][label]
    return tube_skeleton, rawdataf

def backward_compatibility_tree_structure(tree_raw_data):
    """
    Fix backward compatibility
    :param tree_raw_data:
    :return: fixed tree_raw_data
    """
    if 'graph' in tree_raw_data:
        gr = tree_raw_data.pop('graph')
        tree_raw_data['Graph'] = gr  # {'tree1':gr}

    if "Graph" not in tree_raw_data:
        tree_raw_data = {"Graph": tree_raw_data}

    # if all keys in Graph a
    if all([type(k) != str for k in tree_raw_data['Graph'].keys()]):
        gr = tree_raw_data.pop('Graph')
        tree_raw_data['Graph'] = {'tree1': gr}

    # else:
    #     tree_raw_data = tree_raw_data['Graph']
    return tree_raw_data

def main():
    logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel(logging.WARNING)

    # input parser
    parser = argparse.ArgumentParser(
        description='Histology analyser reporter. Try: \
python src/tb_volume.py -i ./tests/hist_stats_test.yaml'
    )
    parser.add_argument(
        '-i', '--inputfile',
        default=None,
        required=True,
        help='input file, yaml file'
    )
    parser.add_argument(
        '-o', '--outputfile',
        default=None,
        help='output file, .raw, .dcm, .tiff, given by extension '
    )
    parser.add_argument(
        '-ot', '--outputfiletype',
        default='pkl',
        help='output file type.  raw, dcm, tiff, or pkl,   default is pkl, '
    )
    parser.add_argument(
        '-vs', '--voxelsize',
        default=[1.0, 1.0, 1.0],
        type=float,
        metavar='N',
        nargs='+',
        help='size of voxel (ZYX)'
    )
    parser.add_argument(
        '-ds', '--datashape',
        default=[200, 200, 200],
        type=int,
        metavar='N',
        nargs='+',
        help='size of output data in pixels for each axis (ZYX)'
    )
    parser.add_argument(
        '-g', '--generator',
        default='vol',
        type=str,
        help='Volume or surface model can be generated by use this option. \
                Use "vol", "volume" for volumetric model. For LAR surface model\
                use "lar". For VTK file use "vtk".'
    )
    parser.add_argument(
        '-d', '--debug', action='store_true',
        help='Debug mode')
    parser.add_argument(
        '-l', '--useLar', action='store_true',
        help='Use LAR')
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    startTime = datetime.now()

    generator_params = None
    generator_class = args.generator

    # if args.generator == "vtk":
    #     import gen_vtk_tree
    #     gen_vtk_tree.vt2vtk_file(args.inputfile, args.outputfile)
    #     return

    tg = TreeBuilder(generator_class, generator_params)
    tg.importFromYaml(args.inputfile)
    tg.voxelsize_mm = args.voxelsize
    tg.shape = args.datashape
    tg.use_lar = args.useLar
    data3d = tg.buildTree()

    logger.info("TimeUsed:" + str(datetime.now() - startTime))
    # volume_px = sum(sum(sum(data3d)))
    # volume_mm3 = volume_px * \
    #     (tg.voxelsize_mm[0] * tg.voxelsize_mm[1] * tg.voxelsize_mm[2])
    # logger.info("Volume px:" + str(volume_px))
    # logger.info("Volume mm3:" + str(volume_mm3))

    # vizualizace
    logger.debug("before visualization")
    tg.show()
    logger.debug("after visualization")

    # ukládání do souboru
    if args.outputfile is not None:
        tg.saveToFile(args.outputfile, args.outputfiletype)


# class TreeGenerator(TreeConstructor):
#     """
#     back compatibility
#     """
#     pass

if __name__ == "__main__":
    main()
