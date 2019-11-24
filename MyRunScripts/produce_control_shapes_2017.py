#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True  # disable ROOT internal argument parser
ROOT.gErrorIgnoreLevel = ROOT.kError

from itertools import product

import argparse
import yaml

import logging
logger = logging.getLogger("")

from shape_producer import InputManager

from ROOT import RDataFrame

def setup_logging(output_file, level=logging.DEBUG):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(output_file, "w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Produce shapes for 2016 Standard Model analysis.")

    parser.add_argument(
        "--directory",
        required=True,
        type=str,
        help="Directory with Artus outputs.")
    parser.add_argument(
        "--mm-friend-directory",
        type=str,
        default=[],
        nargs='+',
        help=
        "Directories arranged as Artus output and containing a friend tree for mm."
    )
    parser.add_argument(
        "--et-friend-directory",
        type=str,
        default=[],
        nargs='+',
        help=
        "Directories arranged as Artus output and containing a friend tree for et."
    )
    parser.add_argument(
        "--mt-friend-directory",
        type=str,
        default=[],
        nargs='+',
        help=
        "Directories arranged as Artus output and containing a friend tree for mt."
    )
    parser.add_argument(
        "--tt-friend-directory",
        type=str,
        default=[],
        nargs='+',
        help=
        "Directories arranged as Artus output and containing a friend tree for tt."
    )
    parser.add_argument(
        "--em-friend-directory",
        type=str,
        default=[],
        nargs='+',
        help=
        "Directories arranged as Artus output and containing a friend tree for em."
    )
    parser.add_argument(
        "--fake-factor-friend-directory",
        default=None,
        type=str,
        help=
        "Directory arranged as Artus output and containing friend trees to data files with fake factors."
    )
    parser.add_argument(
        "--datasets", required=True, type=str, help="Kappa datsets database.")
    parser.add_argument(
        "--binning", required=True, type=str, help="Binning configuration.")
    parser.add_argument(
        "--channels",
        default=[],
        nargs='+',
        type=str,
        help="Channels to be considered.")
    parser.add_argument(
        "--QCD-extrap-fit",
        default=False,
        action='store_true',
        help="Create shapes for QCD extrapolation factor determination.")
    parser.add_argument(
        "--HIG16043",
        action="store_true",
        default=False,
        help="Create shapes of HIG16043 reference analysis.")
    parser.add_argument(
        "--num-threads",
        default=20,
        type=int,
        help="Number of threads to be used.")
    parser.add_argument(
        "--backend",
        default="classic",
        choices=["classic", "tdf"],
        type=str,
        help="Backend. Use classic or tdf.")
    return parser.parse_args()


def main(args):

    # Input files
    directory = args.directory
    et_friend_directory = args.et_friend_directory
    mt_friend_directory = args.mt_friend_directory
    tt_friend_directory = args.tt_friend_directory
    em_friend_directory = args.em_friend_directory
    mm_friend_directory = args.mm_friend_directory

    #### No sources
    #regex = directory + '*/*/'
    #rdf = RDataFrame('ntuple', regex)
    #col_dict = rdf.AsNumpy()
    #print(len(col_dict.items()))

    #### Python
    #chain = InputManager.CreateTChainFromPath('ntuple', directory)
    # debug
    #files_list = chain.GetListOfFiles()
    #for name in files_list:
        #print(name)
    #rdf = RDataFrame(chain)
    #col_dict = rdf.AsNumpy()
    #print(len(col_dict.items()))
    #print(col_dict)

    ##### C++
    ROOT.gInterpreter.Declare('''#include "ShapeProducer/inc/InputManager.h"''')
    ROOT.gInterpreter.Declare('''#include "ShapeProducer/inc/Cutstring.h"''')
    ROOT.gInterpreter.Declare('''#include "ShapeProducer/inc/CutScheduler.h"''')
    Cut = ROOT.Cut
    CutScheduler = ROOT.CutScheduler

    chain = ROOT.shape_producer.InputManager.CreateTChainFromPath(str, str, str, str)('ntuple', directory, *et_friend_directory)
    rdf = RDataFrame(chain)

    cut_vec = ROOT.std.vector['Cut'](
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15 && iso_2<0.15", "muon_iso"), Cut(
                "q_1*q_2<0", "os"),
            Cut("m_vis > 50","m_vis_cut"),
            Cut("(trg_singlemuon_27==1 || trg_singlemuon_24==1)", "trg_selection")
            )



if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("produce_shapes.log", logging.INFO)
    main(args)
