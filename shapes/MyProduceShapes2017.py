#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True  # disable ROOT internal argument parser
ROOT.gErrorIgnoreLevel = ROOT.kError

#from ShapeProducer.cutstring import Cut, Cuts, Weight
#from ShapeProducer.systematics import Systematics, Systematic
#from ShapeProducer.categories import Category
#from ShapeProducer.binning import ConstantBinning, VariableBinning
#from ShapeProducer.variable import Variable
#from ShapeProducer.systematic_variations import Nominal, DifferentPipeline, SquareAndRemoveWeight, create_systematic_variations, AddWeight, ReplaceWeight, Relabel
#from ShapeProducer.process import Process
#from ShapeProducer.estimation_methods import AddHistogramEstimationMethod
#from ShapeProducer.channel import ETSM2017, MTSM2017, TTSM2017, EMSM2017

from itertools import product

import argparse
import yaml

import logging
logger = logging.getLogger()


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
        description="Produce shapes for 2017 Standard Model analysis.")

    parser.add_argument(
        "--directory",
        required=True,
        type=str,
        help="Directory with Artus outputs.")
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
        "--fake-factor-friend-directory",
        default=None,
        type=str,
        help=
        "Directory arranged as Artus output and containing friend trees to data files with fake factors."
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
        "--QCD-extrap-fit",
        default=False,
        action='store_true',
        help="Create shapes for QCD extrapolation factor determination.")
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
    parser.add_argument("--era", type=str, help="Experiment era.")
    parser.add_argument(
        "--gof-channel",
        default=None,
        type=str,
        help="Channel for goodness of fit shapes.")
    parser.add_argument(
        "--gof-variable",
        type=str,
        help="Variable for goodness of fit shapes.")
    parser.add_argument(
        "--num-threads",
        default=32,
        type=int,
        help="Number of threads to be used.")
    parser.add_argument(
        "--backend",
        default="classic",
        choices=["classic", "tdf"],
        type=str,
        help="Backend. Use classic or tdf.")
    parser.add_argument(
        "--tag", default="ERA_CHANNEL", type=str, help="Tag of output files.")
    parser.add_argument(
        "--skip-systematic-variations",
        #default=False,
        default=True,
        type=str,
        help="Do not produce the systematic variations.")
    return parser.parse_args()


def main(args):
    # Container for all distributions to be drawn
    logger.info("Set up shape variations.")
    '''
    systematics = Systematics(
        "{}_shapes.root".format(args.tag),
        num_threads=args.num_threads,
        skip_systematic_variations=args.skip_systematic_variations)
    '''

    # Channels and processes
    # yapf: disable
    directory = args.directory
    et_friend_directory = args.et_friend_directory
    mt_friend_directory = args.mt_friend_directory
    tt_friend_directory = args.tt_friend_directory
    em_friend_directory = args.em_friend_directory
    ff_friend_directory = args.fake_factor_friend_directory

'''
    # Era selection
    if "2017" in args.era:
        from ShapeProducer.estimation_methods_2017 import DataEstimation, ZTTEstimation, ZTTEmbeddedEstimation, ZLEstimation, ZJEstimation, TTLEstimation, TTJEstimation, TTTEstimation, VVLEstimation, VVTEstimation, VVJEstimation, WEstimation, ggHEstimation, qqHEstimation, VHEstimation, WHEstimation, ZHEstimation, ttHEstimation, QCDEstimation_ABCD_TT_ISO2, QCDEstimation_SStoOS_MTETEM, NewFakeEstimationLT, NewFakeEstimationTT, HWWEstimation, ggHWWEstimation, qqHWWEstimation

        from ShapeProducer.era import Run2017
        era = Run2017(args.datasets)

    else:
        logger.critical("Era {} is not implemented.".format(args.era))
        raise Exception

    # Channels and processes
    # yapf: disable
    directory = args.directory
    et_friend_directory = args.et_friend_directory
    mt_friend_directory = args.mt_friend_directory
    tt_friend_directory = args.tt_friend_directory
    em_friend_directory = args.em_friend_directory
    ff_friend_directory = args.fake_factor_friend_directory
    mt = MTSM2017()
    if args.QCD_extrap_fit:
        mt.cuts.remove("muon_iso")
        mt.cuts.add(Cut("(iso_1<0.5)*(iso_1>=0.15)", "muon_iso_loose"))
    mt_processes = {
        "data"  : Process("data_obs", DataEstimation      (era, directory, mt, friend_directory=mt_friend_directory)),
        "ZTT"   : Process("ZTT",      ZTTEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),
        "EMB"   : Process("EMB",      ZTTEmbeddedEstimation  (era, directory, mt, friend_directory=mt_friend_directory)),
        "ZJ"    : Process("ZJ",       ZJEstimation        (era, directory, mt, friend_directory=mt_friend_directory)),
        "ZL"    : Process("ZL",       ZLEstimation        (era, directory, mt, friend_directory=mt_friend_directory)),
        "TTT"   : Process("TTT",      TTTEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),
        "TTJ"   : Process("TTJ",      TTJEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),
        "TTL"   : Process("TTL",      TTLEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),
        "VVT"   : Process("VVT",      VVTEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),
        "VVJ"   : Process("VVJ",      VVJEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),
        "VVL"   : Process("VVL",      VVLEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),
        "W"     : Process("W",        WEstimation         (era, directory, mt, friend_directory=mt_friend_directory)),

        "VH125"    : Process("VH125",    VHEstimation        (era, directory, mt, friend_directory=mt_friend_directory)),
        "WH125"    : Process("WH125",    WHEstimation        (era, directory, mt, friend_directory=mt_friend_directory)),
        "ZH125"    : Process("ZH125",    ZHEstimation        (era, directory, mt, friend_directory=mt_friend_directory)),
        "ttH125"   : Process("ttH125",   ttHEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),

        "ggHWW125" : Process("ggHWW125", ggHWWEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),
        "qqHWW125" : Process("qqHWW125", qqHWWEstimation       (era, directory, mt, friend_directory=mt_friend_directory)),
        }

    # Stage 0 and 1.1 signals for ggH & qqH
    for ggH_htxs in ggHEstimation.htxs_dict:
        mt_processes[ggH_htxs] = Process(ggH_htxs, ggHEstimation(ggH_htxs, era, directory, mt, friend_directory=mt_friend_directory))
    for qqH_htxs in qqHEstimation.htxs_dict:
        mt_processes[qqH_htxs] = Process(qqH_htxs, qqHEstimation(qqH_htxs, era, directory, mt, friend_directory=mt_friend_directory))

    mt_processes["FAKES"] = Process("jetFakes", NewFakeEstimationLT(era, directory, mt, [mt_processes[process] for process in ["EMB", "ZL", "TTL", "VVL"]], mt_processes["data"], friend_directory=mt_friend_directory+[ff_friend_directory]))
    mt_processes["QCD"] = Process("QCD", QCDEstimation_SStoOS_MTETEM(era, directory, mt,
            [mt_processes[process] for process in ["ZTT", "ZL", "ZJ", "W", "TTT", "TTJ", "TTL", "VVT", "VVJ", "VVL"]],
            mt_processes["data"], friend_directory=mt_friend_directory, extrapolation_factor=1.00))


    et = ETSM2017()
    et_processes = {
        "data"  : Process("data_obs", DataEstimation      (era, directory, et, friend_directory=et_friend_directory)),
        "ZTT"   : Process("ZTT",      ZTTEstimation       (era, directory, et, friend_directory=et_friend_directory)),
        "EMB"   : Process("EMB",      ZTTEmbeddedEstimation  (era, directory, et, friend_directory=et_friend_directory)),
        "ZJ"    : Process("ZJ",       ZJEstimation        (era, directory, et, friend_directory=et_friend_directory)),
        "ZL"    : Process("ZL",       ZLEstimation        (era, directory, et, friend_directory=et_friend_directory)),
        "TTT"   : Process("TTT",      TTTEstimation       (era, directory, et, friend_directory=et_friend_directory)),
        "TTJ"   : Process("TTJ",      TTJEstimation       (era, directory, et, friend_directory=et_friend_directory)),
        "TTL"   : Process("TTL",      TTLEstimation       (era, directory, et, friend_directory=et_friend_directory)),
        "VVT"   : Process("VVT",      VVTEstimation       (era, directory, et, friend_directory=et_friend_directory)),
        "VVJ"   : Process("VVJ",      VVJEstimation       (era, directory, et, friend_directory=et_friend_directory)),
        "VVL"   : Process("VVL",      VVLEstimation       (era, directory, et, friend_directory=et_friend_directory)),
        "W"     : Process("W",        WEstimation         (era, directory, et, friend_directory=et_friend_directory)),

        "VH125"    : Process("VH125",    VHEstimation        (era, directory, et, friend_directory=et_friend_directory)),
        "WH125"    : Process("WH125",    WHEstimation        (era, directory, et, friend_directory=et_friend_directory)),
        "ZH125"    : Process("ZH125",    ZHEstimation        (era, directory, et, friend_directory=et_friend_directory)),
        "ttH125"   : Process("ttH125",   ttHEstimation       (era, directory, et, friend_directory=et_friend_directory)),

        "ggHWW125" : Process("ggHWW125", ggHWWEstimation       (era, directory, et, friend_directory=et_friend_directory)),
        "qqHWW125" : Process("qqHWW125", qqHWWEstimation       (era, directory, et, friend_directory=et_friend_directory)),
        }

    # Stage 0 and 1.1 signals for ggH & qqH
    for ggH_htxs in ggHEstimation.htxs_dict:
        et_processes[ggH_htxs] = Process(ggH_htxs, ggHEstimation(ggH_htxs, era, directory, et, friend_directory=et_friend_directory))
    for qqH_htxs in qqHEstimation.htxs_dict:
        et_processes[qqH_htxs] = Process(qqH_htxs, qqHEstimation(qqH_htxs, era, directory, et, friend_directory=et_friend_directory))

    et_processes["FAKES"] = Process("jetFakes", NewFakeEstimationLT(era, directory, et, [et_processes[process] for process in ["EMB", "ZL", "TTL", "VVL"]], et_processes["data"], friend_directory=et_friend_directory+[ff_friend_directory]))
    et_processes["QCD"] = Process("QCD", QCDEstimation_SStoOS_MTETEM(era, directory, et,
            [et_processes[process] for process in ["ZTT", "ZL", "ZJ", "W", "TTT", "TTJ", "TTL", "VVT", "VVJ", "VVL"]],
            et_processes["data"], friend_directory=et_friend_directory, extrapolation_factor=1.00))


    tt = TTSM2017()
    tt_processes = {
        "data"  : Process("data_obs", DataEstimation      (era, directory, tt, friend_directory=tt_friend_directory)),
        "ZTT"   : Process("ZTT",      ZTTEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),
        "EMB"   : Process("EMB",      ZTTEmbeddedEstimation  (era, directory, tt, friend_directory=tt_friend_directory)),
        "ZJ"    : Process("ZJ",       ZJEstimation        (era, directory, tt, friend_directory=tt_friend_directory)),
        "ZL"    : Process("ZL",       ZLEstimation        (era, directory, tt, friend_directory=tt_friend_directory)),
        "TTT"   : Process("TTT",      TTTEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),
        "TTJ"   : Process("TTJ",      TTJEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),
        "TTL"   : Process("TTL",      TTLEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),
        "VVT"   : Process("VVT",      VVTEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),
        "VVJ"   : Process("VVJ",      VVJEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),
        "VVL"   : Process("VVL",      VVLEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),
        "W"     : Process("W",        WEstimation         (era, directory, tt, friend_directory=tt_friend_directory)),

        "VH125"    : Process("VH125",    VHEstimation        (era, directory, tt, friend_directory=tt_friend_directory)),
        "WH125"    : Process("WH125",    WHEstimation        (era, directory, tt, friend_directory=tt_friend_directory)),
        "ZH125"    : Process("ZH125",    ZHEstimation        (era, directory, tt, friend_directory=tt_friend_directory)),
        "ttH125"   : Process("ttH125",   ttHEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),

        "ggHWW125" : Process("ggHWW125", ggHWWEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),
        "qqHWW125" : Process("qqHWW125", qqHWWEstimation       (era, directory, tt, friend_directory=tt_friend_directory)),
        }

    # Stage 0 and 1.1 signals for ggH & qqH
    for ggH_htxs in ggHEstimation.htxs_dict:
        tt_processes[ggH_htxs] = Process(ggH_htxs, ggHEstimation(ggH_htxs, era, directory, tt, friend_directory=tt_friend_directory))
    for qqH_htxs in qqHEstimation.htxs_dict:
        tt_processes[qqH_htxs] = Process(qqH_htxs, qqHEstimation(qqH_htxs, era, directory, tt, friend_directory=tt_friend_directory))

    tt_processes["FAKES"] = Process("jetFakes", NewFakeEstimationTT(era, directory, tt, [tt_processes[process] for process in ["EMB", "ZL", "TTL", "VVL"]], tt_processes["data"], friend_directory=tt_friend_directory+[ff_friend_directory]))
    tt_processes["QCD"] = Process("QCD", QCDEstimation_ABCD_TT_ISO2(era, directory, tt,
            [tt_processes[process] for process in ["ZTT", "ZL", "ZJ", "W", "TTT", "TTJ", "TTL", "VVT", "VVJ", "VVL"]],
            tt_processes["data"], friend_directory=tt_friend_directory))

    em = EMSM2017()
    em_processes = {
        "data"  : Process("data_obs", DataEstimation      (era, directory, em, friend_directory=em_friend_directory)),
        "ZTT"   : Process("ZTT",      ZTTEstimation       (era, directory, em, friend_directory=em_friend_directory)),
        "EMB"   : Process("EMB",      ZTTEmbeddedEstimation  (era, directory, em, friend_directory=em_friend_directory)),
        "ZL"    : Process("ZL",       ZLEstimation        (era, directory, em, friend_directory=em_friend_directory)),
        "TTT"   : Process("TTT",      TTTEstimation       (era, directory, em, friend_directory=em_friend_directory)),
        "TTL"   : Process("TTL",      TTLEstimation       (era, directory, em, friend_directory=em_friend_directory)),
        "VVT"   : Process("VVT",      VVTEstimation       (era, directory, em, friend_directory=em_friend_directory)),
        "VVL"   : Process("VVL",      VVLEstimation       (era, directory, em, friend_directory=em_friend_directory)),
        "W"     : Process("W",        WEstimation         (era, directory, em, friend_directory=em_friend_directory)),

        "VH125"    : Process("VH125",    VHEstimation        (era, directory, em, friend_directory=em_friend_directory)),
        "WH125"    : Process("WH125",    WHEstimation        (era, directory, em, friend_directory=em_friend_directory)),
        "ZH125"    : Process("ZH125",    ZHEstimation        (era, directory, em, friend_directory=em_friend_directory)),
        "ttH125"   : Process("ttH125",   ttHEstimation       (era, directory, em, friend_directory=em_friend_directory)),

        "ggHWW125" : Process("ggHWW125", ggHWWEstimation       (era, directory, em, friend_directory=em_friend_directory)),
        "qqHWW125" : Process("qqHWW125", qqHWWEstimation       (era, directory, em, friend_directory=em_friend_directory)),
        }

    # Stage 0 and 1.1 signals for ggH & qqH
    for ggH_htxs in ggHEstimation.htxs_dict:
        em_processes[ggH_htxs] = Process(ggH_htxs, ggHEstimation(ggH_htxs, era, directory, em, friend_directory=em_friend_directory))
    for qqH_htxs in qqHEstimation.htxs_dict:
        em_processes[qqH_htxs] = Process(qqH_htxs, qqHEstimation(qqH_htxs, era, directory, em, friend_directory=em_friend_directory))

    em_processes["QCD"] = Process("QCD", QCDEstimation_SStoOS_MTETEM(era, directory, em, [em_processes[process] for process in ["EMB", "ZL", "W", "TTL", "VVL"]], em_processes["data"], extrapolation_factor=1.0, qcd_weight = Weight("em_qcd_extrap_up_Weight","qcd_weight")))


    # Variables and categories
    binning = yaml.load(open(args.binning))

    et_categories = []
    # Analysis shapes
    if "et" in args.channels:
        classes_et = ["ggh", "qqh", "ztt", "zll", "w", "tt", "ss", "misc"]
        for i, label in enumerate(classes_et):
            score = Variable(
                "et_max_score",
                 VariableBinning(binning["analysis"]["et"][label]))
            et_categories.append(
                Category(
                    label,
                    et,
                    Cuts(
                        Cut("et_max_index=={index}".format(index=i), "exclusive_score")),
                    variable=score))
            if label in ["ggh", "qqh"]:
                stxs = 100 if label == "ggh" else 200
                for i_e, e in enumerate(binning["stxs_stage1p1"][label]):
                    score = Variable(
                        "et_max_score",
                         VariableBinning(binning["analysis"]["et"][label]))
                    et_categories.append(
                        Category(
                            "{}_{}".format(label,str(stxs+i_e)),
                            et,
                            Cuts(Cut("et_max_index=={index}".format(index=i), "exclusive_score"),
                                 Cut(e, "stxs_stage1p1_cut")),
                            variable=score))
    # Goodness of fit shapes
    elif args.gof_channel == "et":
        score = Variable(
                args.gof_variable,
                VariableBinning(binning["gof"]["et"][args.gof_variable]["bins"]),
                expression=binning["gof"]["et"][args.gof_variable]["expression"])
        if "cut" in binning["gof"]["et"][args.gof_variable].keys():
            cuts=Cuts(Cut(binning["gof"]["et"][args.gof_variable]["cut"], "binning"))
        else:
            cuts=Cuts()
        et_categories.append(
            Category(
                args.gof_variable,
                et,
                cuts,
                variable=score))

    mt_categories = []
    # Analysis shapes
    if "mt" in args.channels:
        classes_mt = ["ggh", "qqh", "ztt", "zll", "w", "tt", "ss", "misc"]
        for i, label in enumerate(classes_mt):
            score = Variable(
                "mt_max_score",
                 VariableBinning(binning["analysis"]["mt"][label]))
            mt_categories.append(
                Category(
                    label,
                    mt,
                    Cuts(
                        Cut("mt_max_index=={index}".format(index=i), "exclusive_score")),
                    variable=score))
            if label in ["ggh", "qqh"]:
                stxs = 100 if label == "ggh" else 200
                for i_e, e in enumerate(binning["stxs_stage1p1"][label]):
                    score = Variable(
                        "mt_max_score",
                         VariableBinning(binning["analysis"]["mt"][label]))
                    mt_categories.append(
                        Category(
                            "{}_{}".format(label,str(stxs+i_e)),
                            mt,
                            Cuts(Cut("mt_max_index=={index}".format(index=i), "exclusive_score"),
                                 Cut(e, "stxs_stage1p1_cut")),
                            variable=score))
    # Goodness of fit shapes
    elif args.gof_channel == "mt":
        score = Variable(
                args.gof_variable,
                VariableBinning(binning["gof"]["mt"][args.gof_variable]["bins"]),
                expression=binning["gof"]["mt"][args.gof_variable]["expression"])
        if "cut" in binning["gof"]["mt"][args.gof_variable].keys():
            cuts=Cuts(Cut(binning["gof"]["mt"][args.gof_variable]["cut"], "binning"))
        else:
            cuts=Cuts()
        mt_categories.append(
            Category(
                args.gof_variable,
                mt,
                cuts,
                variable=score))

    tt_categories = []
    # Analysis shapes
    if "tt" in args.channels:
        classes_tt = ["ggh", "qqh", "ztt", "noniso", "misc"]
        for i, label in enumerate(classes_tt):
            score = Variable(
                "tt_max_score",
                 VariableBinning(binning["analysis"]["tt"][label]))
            tt_categories.append(
                Category(
                    label,
                    tt,
                    Cuts(
                        Cut("tt_max_index=={index}".format(index=i), "exclusive_score")),
                    variable=score))
            if label in ["ggh", "qqh"]:
                stxs = 100 if label == "ggh" else 200
                for i_e, e in enumerate(binning["stxs_stage1p1"][label]):
                    score = Variable(
                        "tt_max_score",
                         VariableBinning(binning["analysis"]["tt"][label]))
                    tt_categories.append(
                        Category(
                            "{}_{}".format(label,str(stxs+i_e)),
                            tt,
                            Cuts(Cut("tt_max_index=={index}".format(index=i), "exclusive_score"),
                                 Cut(e, "stxs_stage1p1_cut")),
                            variable=score))
    # Goodness of fit shapes
    elif args.gof_channel == "tt":
        score = Variable(
                args.gof_variable,
                VariableBinning(binning["gof"]["tt"][args.gof_variable]["bins"]),
                expression=binning["gof"]["tt"][args.gof_variable]["expression"])
        if "cut" in binning["gof"]["tt"][args.gof_variable].keys():
            cuts=Cuts(Cut(binning["gof"]["tt"][args.gof_variable]["cut"], "binning"))
        else:
            cuts=Cuts()
        tt_categories.append(
            Category(
                args.gof_variable,
                tt,
                cuts,
                variable=score))
    em_categories = []
    # Analysis shapes
    if "em" in args.channels:
        classes_em = ["ggh", "qqh", "ztt", "tt", "ss", "misc", "db"]
        for i, label in enumerate(classes_em):
            score = Variable(
                "em_max_score",
                 VariableBinning(binning["analysis"]["em"][label]))
            em_categories.append(
                Category(
                    label,
                    em,
                    Cuts(
                        Cut("em_max_index=={index}".format(index=i), "exclusive_score")),
                    variable=score))
            if label in ["ggh", "qqh"]:
                stxs = 100 if label == "ggh" else 200
                for i_e, e in enumerate(binning["stxs_stage1p1"][label]):
                    score = Variable(
                        "em_max_score",
                         VariableBinning(binning["analysis"]["em"][label]))
                    em_categories.append(
                        Category(
                            "{}_{}".format(label,str(stxs+i_e)),
                            em,
                            Cuts(Cut("em_max_index=={index}".format(index=i), "exclusive_score"),
                                 Cut(e, "stxs_stage1p1_cut")),
                            variable=score))
    # Goodness of fit shapes
    elif args.gof_channel == "em":
        score = Variable(
                args.gof_variable,
                VariableBinning(binning["gof"]["em"][args.gof_variable]["bins"]),
                expression=binning["gof"]["em"][args.gof_variable]["expression"])
        if "cut" in binning["gof"]["em"][args.gof_variable].keys():
            cuts=Cuts(Cut(binning["gof"]["em"][args.gof_variable]["cut"], "binning"))
        else:
            cuts=Cuts()
        em_categories.append(
            Category(
                args.gof_variable,
                em,
                cuts,
                variable=score))



    # Nominal histograms
    signal_nicks = ["WH125", "ZH125", "VH125", "ttH125"]
    ww_nicks = ["ggHWW125", "qqHWW125"]
    if args.gof_channel == None:
        signal_nicks += [ggH_htxs for ggH_htxs in ggHEstimation.htxs_dict] + [qqH_htxs for qqH_htxs in qqHEstimation.htxs_dict] + ww_nicks
    else:
        signal_nicks +=  ["ggH125", "qqH125"] + ww_nicks

    # yapf: enable
    if "et" in [args.gof_channel] + args.channels:
        for process, category in product(et_processes.values(), et_categories):
            systematics.add(
                Systematic(
                    category=category,
                    process=process,
                    analysis="smhtt",
                    era=era,
                    variation=Nominal(),
                    mass="125"))
    if "mt" in [args.gof_channel] + args.channels:
        for process, category in product(mt_processes.values(), mt_categories):
            systematics.add(
                Systematic(
                    category=category,
                    process=process,
                    analysis="smhtt",
                    era=era,
                    variation=Nominal(),
                    mass="125"))
    if "tt" in [args.gof_channel] + args.channels:
        for process, category in product(tt_processes.values(), tt_categories):
            systematics.add(
                Systematic(
                    category=category,
                    process=process,
                    analysis="smhtt",
                    era=era,
                    variation=Nominal(),
                    mass="125"))
    if "em" in [args.gof_channel] + args.channels:
        for process, category in product(em_processes.values(), em_categories):
            systematics.add(
                Systematic(
                    category=category,
                    process=process,
                    analysis="smhtt",
                    era=era,
                    variation=Nominal(),
                    mass="125"))

    # Produce histograms
    logger.info("Start producing shapes.")
    systematics.produce()
    logger.info("Done producing shapes.")
'''

if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("{}_produce_shapes.log".format(args.tag), logging.INFO)
    main(args)
