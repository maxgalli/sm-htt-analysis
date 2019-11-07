# -*- coding: utf-8 -*-

import copy
import os

from estimation_methods import EstimationMethod, SStoOSEstimationMethod, ABCDEstimationMethod, SumUpEstimationMethod, NewFakeEstimationMethodLT, NewFakeEstimationMethodTT
from histogram import *
from cutstring import *
from process import *
from systematics import *
from systematic_variations import *
from era import log_query

import logging
logger = logging.getLogger(__name__)

ggH_htxs = {
    "ggH125": "(htxs_stage1p1cat>=100)&&(htxs_stage1p1cat<=113)",
    "ggH_GG2H_FWDH125": "htxs_stage1p1cat == 100",
    "ggH_GG2H_PTH_GT200125": "htxs_stage1p1cat == 101",
    "ggH_GG2H_0J_PTH_0_10125": "htxs_stage1p1cat == 102",
    "ggH_GG2H_0J_PTH_GT10125": "htxs_stage1p1cat == 103",
    "ggH_GG2H_1J_PTH_0_60125": "htxs_stage1p1cat == 104",
    "ggH_GG2H_1J_PTH_60_120125": "htxs_stage1p1cat == 105",
    "ggH_GG2H_1J_PTH_120_200125": "htxs_stage1p1cat == 106",
    "ggH_GG2H_GE2J_MJJ_0_350_PTH_0_60125": "htxs_stage1p1cat == 107",
    "ggH_GG2H_GE2J_MJJ_0_350_PTH_60_120125": "htxs_stage1p1cat == 108",
    "ggH_GG2H_GE2J_MJJ_0_350_PTH_120_200125": "htxs_stage1p1cat == 109",
    "ggH_GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25125":
    "htxs_stage1p1cat == 110",
    "ggH_GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25125":
    "htxs_stage1p1cat == 111",
    "ggH_GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25125":
    "htxs_stage1p1cat == 112",
    "ggH_GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25125":
    "htxs_stage1p1cat == 113",
}

qqH_htxs = {
    "qqH125": "(htxs_stage1p1cat>=200)&&(htxs_stage1p1cat<=210)",
    "qqH_QQ2HQQ_FWDH125": "htxs_stage1p1cat == 200",
    "qqH_QQ2HQQ_0J125": "htxs_stage1p1cat == 201",
    "qqH_QQ2HQQ_1J125": "htxs_stage1p1cat == 202",
    "qqH_QQ2HQQ_GE2J_MJJ_0_60125": "htxs_stage1p1cat == 203",
    "qqH_QQ2HQQ_GE2J_MJJ_60_120125": "htxs_stage1p1cat == 204",
    "qqH_QQ2HQQ_GE2J_MJJ_120_350125": "htxs_stage1p1cat == 205",
    "qqH_QQ2HQQ_GE2J_MJJ_GT350_PTH_GT200125": "htxs_stage1p1cat == 206",
    "qqH_QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25125":
    "htxs_stage1p1cat == 207",
    "qqH_QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25125":
    "htxs_stage1p1cat == 208",
    "qqH_QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25125":
    "htxs_stage1p1cat == 209",
    "qqH_QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25125":
    "htxs_stage1p1cat == 210",
}

# TODO fix all the weights
# Definition of global weights


def get_triggerweight_for_channel(channel):
    weight = Weight("1.0", "triggerweight")

    singleMC = "singleTriggerMCEfficiencyWeightKIT_1"
    crossMC = "crossTriggerMCEfficiencyWeight_1"
    MCTau_1 = "((byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5 && byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)*crossTriggerMCEfficiencyWeight_vloose_MVA_1 + (byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)*crossTriggerMCEfficiencyWeight_tight_MVA_1)"
    MCTau_2 = MCTau_1.replace("_1", "_2")

    if "mt" in channel:
        trig_sL = "(trg_singlemuon)"
        trig_X = "(pt_1 > 20 && pt_1 < 21 && trg_mutaucross)"
        MCTau = "dummy"

        MuTauMC = "*".join([trig_sL, singleMC]) + "+" + "*".join(
            [trig_X, crossMC])  # TODO add TauCrossTrigger Weight back in
        MuTauData = MuTauMC.replace("MC", "Data")
        MuTau = "(" + MuTauData + ")/(" + MuTauMC + ")"
        weight = Weight(MuTau, "triggerWeight")

    elif "et" in channel:
        trig_sL = "(trg_singleelectron)"
        trig_X = "blank"  # TODO add ElTau Crosstrigger

        ElTauMC = "*".join([trig_sL, singleMC
                            ])  # + "+" + "*".join([trig_X, crossMCL, MCTau_2])
        ElTauData = ElTauMC.replace("MC", "Data")
        ElTau = "(" + ElTauData + ")/(" + ElTauMC + ")"
        weight = Weight(ElTau, "triggerweight")

    elif "tt" in channel:  # TODO add TauTrigger SF
        #DiTauMC = "*".join([MCTau_1,MCTau_2])
        #DiTauData = DiTauMC.replace("MC","Data")
        #DiTau = "("+DiTauData+")/("+DiTauMC+")"
        #weight = Weight(DiTau,"triggerweight")
        weight = Weight("triggerWeight_1 * triggerWeight_2", "triggerweight")

    # elif "em" in channel:
    #     weight = Weight(
    #         "(trigger_23_data_Weight_2*trigger_12_data_Weight_1*(trg_muonelectron_mu23ele12==1)+trigger_23_data_Weight_1*trigger_8_data_Weight_2*(trg_muonelectron_mu8ele23==1) - trigger_23_data_Weight_2*trigger_23_data_Weight_1*(trg_muonelectron_mu8ele23==1 && trg_muonelectron_mu23ele12==1))/(trigger_23_mc_Weight_2*trigger_12_mc_Weight_1*(trg_muonelectron_mu23ele12==1)+trigger_23_mc_Weight_1*trigger_8_mc_Weight_2*(trg_muonelectron_mu8ele23==1) - trigger_23_mc_Weight_2*trigger_23_mc_Weight_1*(trg_muonelectron_mu8ele23==1 && trg_muonelectron_mu23ele12==1))",
    #         "trigger_lepton_sf")

    return weight


def get_singlelepton_triggerweight_for_channel(channel):
    weight = Weight("1.0", "triggerweight_sl")

    MCTau_1 = "((byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5 && byMediumIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)*crossTriggerMCEfficiencyWeight_medium_MVA_1 + (byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)*crossTriggerMCEfficiencyWeight_tight_MVA_1)"
    MCTau_2 = MCTau_1.replace("_1", "_2")

    # if "mt" in channel or "et" in channel:
    #     weight = Weight("singleTriggerDataEfficiencyWeightKIT_1/singleTriggerMCEfficiencyWeightKIT_1","triggerweight")
    # elif "tt" in channel:
    #     DiTauMC = "*".join([MCTau_1,MCTau_2])
    #     DiTauData = DiTauMC.replace("MC","Data")
    #     DiTau = "("+DiTauData+")/("+DiTauMC+")"
    #     weight = Weight(DiTau,"triggerweight")

    return weight


def get_tauByIsoIdWeight_for_channel(channel):
    # WPs: Tight 0.87, urrently used: SR mt,et Tight; SR tt Tight
    # Source: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#Tau_ID_efficiency
    weight = Weight("1.0", "taubyIsoIdWeight")
    if "ID" in channel.__class__.__name__:  # this is used for the TauID measurements
        return weight
    elif "mt" in channel.name or "et" in channel.name:
        weight = Weight("((gen_match_2 == 5)*0.87 + (gen_match_2 != 5))",
                        "taubyIsoIdWeight")
    elif "tt" in channel.name:
        weight = Weight(
            "((gen_match_1 == 5)*0.87 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.87 + (gen_match_2 != 5))",
            "taubyIsoIdWeight")
    return weight


def get_eleHLTZvtxWeight_for_channel(channel):  #TODO is this needed ?
    weight = Weight("1.0", "eleHLTZvtxWeight")
    # if "et" in channel:
    #     weight = Weight("(trg_singleelectron_35 || trg_singleelectron_32 || trg_singleelectron_27 || trg_crossele_ele24tau30)*0.991 + (!(trg_singleelectron_35 || trg_singleelectron_32 || trg_singleelectron_27 || trg_crossele_ele24tau30))*1.0", "eleHLTZvtxWeight")
    return weight


def get_eleRecoWeight_for_channel(channel):
    weight = Weight("1.0", "eleRecoWeight")
    if "et" in channel:
        weight = Weight("(eleRecoWeight_1)", "eleRecoWeight")
    if "em" in channel:
        weight = Weight("(eleRecoWeight_1)",
                        "eleRecoWeight")  # TODO all this in next ntuples
    return weight


class DataEstimation(EstimationMethod):
    def __init__(self,
                 era,
                 directory,
                 channel,
                 friend_directory=None,
                 folder="nominal"):
        super(DataEstimation, self).__init__(name="data_obs",
                                             folder=folder,
                                             era=era,
                                             directory=directory,
                                             friend_directory=friend_directory,
                                             channel=channel,
                                             mc_campaign=None)
        self._channel = channel

    def get_files(self):
        return self.artus_file_names(self.era.data_files(self._channel))

    def get_cuts(self):
        return Cuts()


class FakeEstimationLT(DataEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(DataEstimation, self).__init__(
            name="fakes",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def get_weights(self):
        return Weights(Weight("ff2_nom", "fake_factor"))

    def create_root_objects(self, systematic):
        aiso_systematic = copy.deepcopy(systematic)
        aiso_systematic.category.cuts.remove("tau_iso")
        aiso_systematic.category.cuts.add(
            Cut(
                "byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5&&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                "tau_aiso"))
        return super(FakeEstimationLT,
                     self).create_root_objects(aiso_systematic)


class NewFakeEstimationLT(NewFakeEstimationMethodLT):
    def __init__(
            self,
            era,
            directory,
            channel,
            nofake_processes,
            data_process,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(NewFakeEstimationLT, self).__init__(
            name="jetFakes",
            folder="nominal",
            era=era,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            nofake_processes=nofake_processes,
            data_process=data_process,
            aisoCut=Cut(
                "byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5&&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                "tau_aiso"),
            fakeWeightstring="ff2_nom")


class NewFakeEstimationTT(NewFakeEstimationMethodTT):
    def __init__(
            self,
            era,
            directory,
            channel,
            nofake_processes,
            data_process,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(NewFakeEstimationTT, self).__init__(
            name="jetFakes",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            nofake_processes=nofake_processes,
            data_process=data_process,
            aisoCut=Cut(
                "(byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5&&byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5&&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)||(byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5&&byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5&&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2>0.5)",
                "tau_aiso"),
            fakeWeightstring=
            "(0.5*ff1_nom*(byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5)+0.5*ff2_nom*(byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5))"
        )


class FakeEstimationTT(DataEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(DataEstimation, self).__init__(
            name="fakes",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def get_weights(self):
        return Weights(
            Weight(
                "(0.5*ff1_nom*(byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5)+0.5*ff2_nom*(byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5))",
                "fake_factor"))

    def create_root_objects(self, systematic):
        aiso_systematic = copy.deepcopy(systematic)
        aiso_systematic.category.cuts.remove("tau_1_iso")
        aiso_systematic.category.cuts.remove("tau_2_iso")
        aiso_systematic.category.cuts.add(
            Cut(
                "(byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5&&byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5&&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)||(byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5&&byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5&&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2>0.5)",
                "tau_aiso"))
        return super(FakeEstimationTT,
                     self).create_root_objects(aiso_systematic)


class HTTEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="HTT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        # if self.channel.name=="em":
        #     return Weights(
        #         Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
        #         Weight("numberGeneratedEventsWeight", "numberGeneratedEventsWeight"),
        #         Weight("triggerWeight_1*triggerWeight_2*identificationWeight_1*identificationWeight_2*trackWeight_1*trackWeight_2", "eventWeight"), #TODO update
        #         Weight("puweight", "puweight"),
        #         Weight("prefiringweight", "prefireWeight"),
        #         self.era.lumi_weight)
        # else:
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("generatorWeight", "generatorWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(^GluGluHToTauTau.*125.*|^VBFHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ggHEstimation(HTTEstimation):
    htxs_dict = ggH_htxs

    def __init__(
            self,
            name,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
    ):
        super(HTTEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        # if self.channel.name=="em":
        #     return Weights(
        #         Weight("ggh_NNLO_weight", "gghNNLO"),
        #         Weight("numberGeneratedEventsWeight", "numberGeneratedEventsWeight"),
        #         Weight("8.8384e-8/numberGeneratedEventsWeight", "ggh_stitching_weight"),
        #         Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
        #         Weight("1.01", "bbh_inclusion_weight"),
        #         Weight("triggerWeight_1*triggerWeight_2*identificationWeight_1*identificationWeight_2*trackWeight_1*trackWeight_2", "eventWeight"),
        #         Weight("puweight", "puweight"), self.era.lumi_weight,
        #         Weight("prefiringweight", "prefireWeight"))
        # else:
        weights = super(ggHEstimation, self).get_weights()
        weights.remove("numberGeneratedEventsWeight")
        weights.add(Weight("8.8384e-8", "numberGeneratedEventsWeight"))
        weights.add(Weight("ggh_NNLO_weight", "gghNNLO"))
        weights.add(Weight("1.01", "bbh_inclusion_weight"))
        return weights

    def get_cuts(self):
        return Cuts(
            Cut(self.htxs_dict.get(self.name, self.htxs_dict["ggH125"]),
                "htxs_match"))

    def get_files(self):
        query = {
            "process": "^GluGluHToTauTau.*125.*",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class qqHEstimation(HTTEstimation):
    htxs_dict = qqH_htxs

    def __init__(
            self,
            name,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
    ):
        super(HTTEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(
            Cut(self.htxs_dict.get(self.name, self.htxs_dict["qqH125"]),
                "htxs_match"))

    def get_files(self):
        query = {
            "process":
            "(^VBFHToTauTau.*125.*|^W(minus|plus)HToTauTau.*125.*|^ZHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class VHEstimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="VH",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(
            Cut("(htxs_stage1p1cat>=300)&&(htxs_stage1p1cat<=404)",
                "htxs_match"))

    def get_files(self):
        query = {
            "process": "(^W(minus|plus)HToTauTau.*125.*|^ZHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class WHEstimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="WH",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(
            Cut("(htxs_stage1p1cat>=300)&&(htxs_stage1p1cat<=304)",
                "htxs_match"))

    def get_files(self):
        query = {
            "process": "(^W(minus|plus)HToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ZHEstimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="ZH",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(
            Cut("(htxs_stage1p1cat>=400)&&(htxs_stage1p1cat<=404)",
                "htxs_match"))

    def get_files(self):
        query = {
            "process": "(^ZHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ttHEstimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="ttH",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_files(self):
        query = {
            "process": "(^ttHJetToTT.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class ggHWWEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
    ):
        super(ggHWWEstimation, self).__init__(
            name="ggHWW125",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("1.1019558", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2", "idweight"),
            Weight("isoWeight_1*isoWeight_2", "isoweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            # MC weights
            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "GluGluHToWWTo2L2Nu_M125", 
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg-JHUgenv628-pythia8", 
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class qqHWWEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
    ):
        super(qqHWWEstimation, self).__init__(
            name="qqHWW125",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("0.0857883", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2", "idweight"),
            Weight("isoWeight_1*isoWeight_2", "isoweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "VBFHToWWTo2L2Nu_M125", 
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg-JHUgenv628-pythia8", 
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class bbH120Estimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="bbH120",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_files(self):
        query = {
            "process": "(^SUSYGluGluToBBHToTauTau.*120$)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class bbH130Estimation(HTTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(HTTEstimation, self).__init__(
            name="bbH130",
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_files(self):
        query = {
            "process": "(^SUSYGluGluToBBHToTauTau.*130$)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class DYJetsToLLEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        self.atNLO = atNLO
        name = "DYJetsToLLNLO" if self.atNLO else "DYJetsToLL"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        # if self.channel.name=="em":
        #     return Weights(
        #          Weight("((genbosonmass >= 50.0) * 4.255812e-05*((npartons == 0 || npartons >= 5)*1.0+(npartons == 1)*0.32123574062076404+(npartons == 2)*0.3314444833963529+(npartons == 3)*0.3389929050626262+(npartons == 4)*0.2785338687268455) + ((genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight))","z_stitching_weight"),
        #         Weight("trackWeight_1*trackWeight_2", "eventWeight"), # TODO will be replaced by get_eleRecoWeight_for_channel in new ntuple
        #         Weight("puweight", "puweight"),
        #         Weight("isoWeight_1*isoWeight_2","isoWeight"),
        #         Weight("idWeight_1*idWeight_2","idWeight"),
        #         Weight("zPtReweightWeight", "zPtReweightWeight"),
        #         Weight("prefiringweight", "prefireWeight"),
        #         self.get_triggerweight_for_channel(self.channel._name),
        #         Weight("(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))","hadronic_tau_sf"),
        #         self.era.lumi_weight)
        # else:
        z_stitching_weight = Weight("(1.0)", "z_stitching_weight")
        if self.atNLO:
            z_stitching_weight = Weight(
                "((genbosonmass >= 50.0) * 5.1551e-05 + (genbosonmass < 50.0)*((abs(crossSectionPerEventWeight - 3.987) < 0.01)*4.6936e-06 + (abs(crossSectionPerEventWeight - 10.01) < 0.01)*3.7568e-06))",
                "z_stitching_weight"
            )  # xsec_NNLO [pb] = 2075.14*3, N_inclusive_NLO = 120762939, xsec_NNLO/N_inclusive_NLO = 5.1551e-05; fraction of negative events in 'generatorWeight'
        else:
            z_stitching_weight = Weight(
                "((genbosonmass >= 50.0) * 4.255812e-05*((npartons == 0 || npartons >= 5)*1.0+(npartons == 1)*0.32123574062076404+(npartons == 2)*0.3314444833963529+(npartons == 3)*0.3389929050626262+(npartons == 4)*0.2785338687268455) + (genbosonmass < 50.0)*((abs(crossSectionPerEventWeight - 3.987) < 0.01)*4.6936e-06 + (abs(crossSectionPerEventWeight - 10.01) < 0.01)*3.7568e-06))",
                "z_stitching_weight")
        return Weights(
            # TODO add triggerweights etc
            Weight("generatorWeight", "generatorWeight"),
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_triggerweight_for_channel(self.channel._name),
            # self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            # self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("zPtReweightWeight", "zPtReweightWeight"),
            z_stitching_weight,
            # lumi weight
            self.era.lumi_weight)

    def get_files(self):
        queryM10 = {
            "process": "DYJetsToLL_M10to50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "version": "v2"
        }
        queryM50_inclusive_2_3jet = {
            "process": "DY(|2|3)JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "version": "v2"
        }
        queryM50_1jet_v1 = {
            "process": "DY1JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "^$",
            "version": "v1"
        }
        queryM50_inc = {
            "process": "DYJetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "^$",
            "version": "v2"
        }
        queryM50_4jet = {
            "process": "DY4JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "version": "v2"
        }
        queryEWKZ = {
            "process": "^EWKZ",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
        }
        queryM50NLO_inc = {
            "process": "DYJetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo\-pythia8",
        }
        if self.atNLO:
            files = self.era.datasets_helper.get_nicks_with_query(queryM50NLO_inc) + \
                    self.era.datasets_helper.get_nicks_with_query(queryEWKZ)
            log_query(self.name, queryM50NLO_inc, files)
        else:
            files = self.era.datasets_helper.get_nicks_with_query(queryM50_inc) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM50_inclusive_2_3jet) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM50_1jet_v1) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM10) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM50_4jet) + \
                    self.era.datasets_helper.get_nicks_with_query(queryEWKZ)
            log_query(self.name, queryM50_inc, files)
        return self.artus_file_names(files)


class EWKZEstimation(DYJetsToLLEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        self.atNLO = atNLO
        super(DYJetsToLLEstimation, self).__init__(
            name="EWKZ",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_files(self):
        query_ewkz = {
            "process": "^EWKZ",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query_ewkz)
        log_query(self.name, query_ewkz, files)
        return self.artus_file_names(files)


class ZTTEstimation(DYJetsToLLEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        self.atNLO = atNLO
        name = "ZTTNLO" if self.atNLO else "ZTT"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            ztt_cut = "gen_match_1==4 && gen_match_2==5"
        elif "et" in self.channel.name:
            ztt_cut = "gen_match_1==3 && gen_match_2==5"
        elif "tt" in self.channel.name:
            ztt_cut = "gen_match_1==5 && gen_match_2==5"
        elif "em" in self.channel.name:
            ztt_cut = "gen_match_1==3 && gen_match_2==4"
        elif "mm" in self.channel.name:
            ztt_cut = "gen_match_1==4 && gen_match_2==4"
        return Cuts(Cut(ztt_cut, "ztt_cut"))


class ZLEstimation(DYJetsToLLEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        self.atNLO = atNLO
        name = "ZLNLO" if self.atNLO else "ZL"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    '''def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "gen_match_2<5"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1<6&&gen_match_2<6&&!(gen_match_1==5&&gen_match_2==5))"
        elif "em" in self.channel.name:
            ct = "0 == 1"
        return Cuts(Cut(ct, "zl_genmatch"))'''
    def get_cuts(self):
        if "mt" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "et" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            emb_veto = "!(gen_match_1==5 && gen_match_2==5)"
            ff_veto = "!(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==4)"
            ff_veto = "(1.0)"
        elif "mm" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==4)"
            ff_veto = "(1.0)"
        return Cuts(Cut("%s && %s" % (emb_veto, ff_veto),
                        "dy_emb_and_ff_veto"))


class ZJEstimation(DYJetsToLLEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        self.atNLO = atNLO
        name = "ZJNLO" if self.atNLO else "ZJ"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "gen_match_2 == 6"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            ct = "0 == 1"
        return Cuts(Cut(ct, "dy_fakes"))


class ZTTEmbeddedEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(ZTTEmbeddedEstimation, self).__init__(
            name="EMB",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            friend_directory=friend_directory,
            directory=directory,
            channel=channel,
            mc_campaign=None)

    def embedding_stitchingweight(self):
        if self.channel.name == 'mt':
            comp_eff_B = "(1.0/0.899)"
            comp_eff_C = "(1.0/0.881)"
            comp_eff_D = "(1.0/0.877)"
            comp_eff_E = "(1.0/0.939)"
            comp_eff_F = "(1.0/0.936)"
            comp_eff_G = "(1.0/0.908)"
            comp_eff_H = "(1.0/0.962)"
            runB = "((run >= 272007) && (run < 275657))*" + comp_eff_B
            runC = "+((run >= 275657) && (run < 276315))*" + comp_eff_C
            runD = "+((run >= 276315) && (run < 276831))*" + comp_eff_D
            runE = "+((run >= 276831) && (run < 277772))*" + comp_eff_E
            runF = "+((run >= 277772) && (run < 278820))*" + comp_eff_F
            runG = "+((run >= 278820) && (run < 280919))*" + comp_eff_G
            runH = "+((run >= 280919) && (run < 284045))*" + comp_eff_H
            return "(" + runB + runC + runD + runE + runF + runG + runH + ")"
        elif self.channel.name == 'et':
            comp_eff_B = "(1.0/0.902)"
            comp_eff_C = "(1.0/0.910)"
            comp_eff_D = "(1.0/0.945)"
            comp_eff_E = "(1.0/0.945)"
            comp_eff_F = "(1.0/0.915)"
            comp_eff_G = "(1.0/0.903)"
            comp_eff_H = "(1.0/0.933)"
            runB = "((run >= 272007) && (run < 275657))*" + comp_eff_B
            runC = "+((run >= 275657) && (run < 276315))*" + comp_eff_C
            runD = "+((run >= 276315) && (run < 276831))*" + comp_eff_D
            runE = "+((run >= 276831) && (run < 277772))*" + comp_eff_E
            runF = "+((run >= 277772) && (run < 278820))*" + comp_eff_F
            runG = "+((run >= 278820) && (run < 280919))*" + comp_eff_G
            runH = "+((run >= 280919) && (run < 284045))*" + comp_eff_H
            return "(" + runB + runC + runD + runE + runF + runG + runH + ")"
        elif self.channel.name == 'tt':
            comp_eff_B = "(1.0/0.897)"
            comp_eff_C = "(1.0/0.908)"
            comp_eff_D = "(1.0/0.950)"
            comp_eff_E = "(1.0/0.861)"
            comp_eff_F = "(1.0/0.941)"
            comp_eff_G = "(1.0/0.908)"
            comp_eff_H = "(1.0/0.949)"
            runB = "((run >= 272007) && (run < 275657))*" + comp_eff_B
            runC = "+((run >= 275657) && (run < 276315))*" + comp_eff_C
            runD = "+((run >= 276315) && (run < 276831))*" + comp_eff_D
            runE = "+((run >= 276831) && (run < 277772))*" + comp_eff_E
            runF = "+((run >= 277772) && (run < 278820))*" + comp_eff_F
            runG = "+((run >= 278820) && (run < 280919))*" + comp_eff_G
            runH = "+((run >= 280919) && (run < 284045))*" + comp_eff_H
            return "(" + runB + runC + runD + runE + runF + runG + runH + ")"
        elif self.channel.name == 'em':
            comp_eff_B = "(1.0/0.891)"
            comp_eff_C = "(1.0/0.910)"
            comp_eff_D = "(1.0/0.953)"
            comp_eff_E = "(1.0/0.947)"
            comp_eff_F = "(1.0/0.942)"
            comp_eff_G = "(1.0/0.906)"
            comp_eff_H = "(1.0/0.950)"
            runB = "((run >= 272007) && (run < 275657))*" + comp_eff_B
            runC = "+((run >= 275657) && (run < 276315))*" + comp_eff_C
            runD = "+((run >= 276315) && (run < 276831))*" + comp_eff_D
            runE = "+((run >= 276831) && (run < 277772))*" + comp_eff_E
            runF = "+((run >= 277772) && (run < 278820))*" + comp_eff_F
            runG = "+((run >= 278820) && (run < 280919))*" + comp_eff_G
            runH = "+((run >= 280919) && (run < 284045))*" + comp_eff_H
            return "(" + runB + runC + runD + runE + runF + runG + runH + ")"
        else:
            log.error(
                "Embedded currently not implemented for channel \"%s\"!" %
                self.channel.name)

    def embedding_tauid(self):
        """         if self.channel.name == "et" or self.channel.name == "mt":
            return Weight(
                "(gen_match_2==5)*0.87+(gen_match_2!=5)",
                "emb_tau_id")  # TODO measure for now MC correction factors
        elif self.channel.name == "tt":
            return Weight(
                "((gen_match_1==5)*0.87+(gen_match_1!=5))*((gen_match_2==5)*0.870+(gen_match_2!=5))",
                "emb_tau_id"),
        else: """
        return Weight("1.0", "emb_tau_id"),

    def get_weights(self):
        emb_weights = Weights(
            self.embedding_tauid(),
            Weight("35.883/5.711", "embscaling"),
            Weight("generatorWeight*(generatorWeight<=1.0)", "simulation_sf"),
            Weight("muonEffTrgWeight*muonEffIDWeight_1*muonEffIDWeight_2",
                   "scale_factor"),
            Weight(self.embedding_stitchingweight(), "2016 stitching weight"),
            Weight("embeddedDecayModeWeight", "decayMode_SF"))
        if self.channel.name == "mt":
            emb_weights.add(
                Weight(
                    "idWeight_1*(triggerWeight_1*(pt_1>23)+((MuTau_TauLeg_DataEfficiencyWeight_2/MuTau_TauLeg_EmbeddedEfficiencyWeight_2)*(pt_1<=23)))*isoWeight_1",
                    "lepton_sf"))
            #emb_weights.add(Weight("1.0", "mutau_crosstriggerweight"))
            emb_weights.add(
                Weight("gen_match_1==4 && gen_match_2==5", "emb_veto"))

        elif self.channel.name == "et":
            #emb_weights.add(Weight("idWeight_1*triggerWeight_1*isoWeight_1","lepton_sf"))
            emb_weights.add(
                Weight("gen_match_1==3 && gen_match_2==5", "emb_veto"))

        elif self.channel.name == "tt":
            #emb_weights.add(Weight("(TriggerDataEfficiencyWeight_1/TriggerEmbeddedEfficiencyWeight_1)*(TriggerDataEfficiencyWeight_2/TriggerEmbeddedEfficiencyWeight_2)""trg_sf"))
            emb_weights.add(
                Weight("gen_match_1==5 && gen_match_2==5", "emb_veto"))

        elif self.channel.name == "em":
            emb_weights.add(
                Weight("(gen_match_1==3 && gen_match_2==4)", "emb_veto")
            )  #TODO add trigger sf as soon as they are included
            emb_weights.add(
                Weight("idWeight_1*isoWeight_1*idWeight_2*isoWeight_2",
                       "leptopn_sf"))
            emb_weights.remove(
                "decayMode_SF"
            )  # embeddedDecayModeWeight is only for tau decay modes
        return emb_weights

    def get_files(self):
        query = {"process": "Embedding2016(B)", "embedded": True}
        if self.channel.name == "mt":
            query["campaign"] = "MuTauFinalState"
            query["scenario"] = "inputDoubleMu94XlegacyminiAOD"
        elif self.channel.name == "et":
            query["campaign"] = "ElTauFinalState"
            query["scenario"] = "inputDoubleMu94XlegacyminiAOD"
        elif self.channel.name == "tt":
            query["campaign"] = "TauTauFinalState"
            query["scenario"] = ".*(v2|v3)"
        elif self.channel.name == "em":
            query["campaign"] = "ElMuFinalState"
            query["scenario"] = ".*(v2|v4)"
        files = self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)

    def get_cuts(self):
        return Cuts(
            Cut(
                "((gen_match_1>2 && gen_match_1<6) &&  (gen_match_2>2 && gen_match_2<6))",
                "dy_genuine_tau"))


class WEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
    ):
        super(WEstimation, self).__init__(
            name="W",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("zPtReweightWeight", "zPtReweightWeight"),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC related weights
            #Weight("numberGeneratedEventsWeight","numberGeneratedEventsWeight"), # to be used only for one inclusive sample
            #Weight("crossSectionPerEventWeight","crossSectionPerEventWeight"), # to be used only for one inclusive sample
            # # xsec_NNLO [pb] = 61526.7, N_inclusive = 86916455, xsec_NNLO/N_inclusive = 0.00070788321 [pb] weights: [1.0, 0.2691615837248596, 0.1532341436287767, 0.03960756033932645, 0.03969970742404736]
            Weight("generatorWeight", "generatorWeight"),
            Weight(
                "((0.00070788321*((npartons <= 0 || npartons >= 5)*1.0 + (npartons == 1)*0.2691615837248596 + (npartons == 2)*0.1532341436287767 + (npartons == 3)*0.03960756033932645 + (npartons == 4)*0.03969970742404736)) * (genbosonmass>=0.0) + numberGeneratedEventsWeight * crossSectionPerEventWeight * (genbosonmass<0.0))",
                "wj_stitching_weight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "W.?JetsToLNu",
            #"process": "WJetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        query = {
            "process": "^EWKW",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class EWKWpEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(EWKWpEstimation, self).__init__(
            name="EWKWp",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        # if self.channel.name=="em":
        #     return Weights(
        #         Weight(
        #             "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
        #             "hadronic_tau_sf"), Weight("numberGeneratedEventsWeight", "numberGeneratedEventsWeight"),
        #         Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
        #         Weight("triggerWeight_1*triggerWeight_2*identificationWeight_1*identificationWeight_2*trackWeight_1*trackWeight_2", "eventWeight"),
        #         Weight("puweight", "puweight"),
        #         Weight(
        #             "(5.190747826298e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
        #             "EWKWp_stitching_weight"),
        #         Weight("prefiringweight", "prefireWeight"),
        #         self.era.lumi_weight)
        # else:
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("zPtReweightWeight", "zPtReweightWeight"),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            Weight("generatorWeight", "generatorWeight"),
            Weight(
                "(5.190747826298e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
                "EWKWp_stitching_weight"), self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^EWKWPlus",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)


class EWKWmEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(EWKWmEstimation, self).__init__(
            name="EWKWm",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        # if self.channel.name=="em":
        #     return Weights(
        #         Weight(
        #             "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
        #             "hadronic_tau_sf"), Weight("numberGeneratedEventsWeight", "numberGeneratedEventsWeight"),
        #         Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
        #         Weight("triggerWeight_1*triggerWeight_2*identificationWeight_1*identificationWeight_2*trackWeight_1*trackWeight_2", "eventWeight"),
        #         Weight("puweight", "puweight"),
        #         Weight(
        #             "(5.190747826298e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
        #             "EWKWp_stitching_weight"),
        #         Weight("prefiringweight", "prefireWeight"),
        #         self.era.lumi_weight)
        # else:
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("generatorWeight", "generatorWeight"),
            Weight(
                "(4.200367267668e-6)/(numberGeneratedEventsWeight*crossSectionPerEventWeight)",
                "EWKW_stitching_weight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^EWKWMinus",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, query, files)
        return self.artus_file_names(files)


class WEstimationRaw(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(WEstimationRaw, self).__init__(
            name="W",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        # if self.channel.name=="em":
        #         return Weights(
        #             Weight(
        #                 "(((npartons == 0 || npartons >= 5)*7.09390278348407e-4) + ((npartons == 1)*1.90063898596475e-4) + ((npartons == 2)*5.8529964471165e-5) + ((npartons == 3)*1.9206444928444e-5) + ((npartons == 4)*1.923548021385e-5))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight)",
        #                 "wj_stitching_weight"),
        #             Weight(
        #                 "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
        #                 "hadronic_tau_sf"), Weight("numberGeneratedEventsWeight", "numberGeneratedEventsWeight"),
        #             Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
        #             Weight("triggerWeight_1*triggerWeight_2*identificationWeight_1*identificationWeight_2*trackWeight_1*trackWeight_2", "eventWeight"),
        #             Weight("puweight", "puweight"),
        #             Weight("prefiringweight", "prefireWeight"),
        #             self.era.lumi_weight)
        # else:
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("generatorWeight", "generatorWeight"),
            Weight(
                "(((npartons == 0 || npartons >= 5)*7.09390278348407e-4) + ((npartons == 1)*1.90063898596475e-4) + ((npartons == 2)*5.8529964471165e-5) + ((npartons == 3)*1.9206444928444e-5) + ((npartons == 4)*1.923548021385e-5))/(numberGeneratedEventsWeight*crossSectionPerEventWeight*sampleStitchingWeight)",
                "wj_stitching_weight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "W.*JetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


# class WEstimation(SumUpEstimationMethod):
#     def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
#             get_triggerweight_for_channel=get_triggerweight_for_channel,
#             get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
#             get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
#             get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
#         super(WEstimation, self).__init__(
#             name="W",
#             folder=folder,
#             get_triggerweight_for_channel=get_triggerweight_for_channel,
#             get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
#             get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
#             get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
#             era=era,
#             directory=directory,
#             friend_directory=friend_directory,
#             channel=channel,
#             processes=[
#                 Process(
#                     "W",
#                     WEstimationRaw(
#                         era,
#                         directory,
#                         channel,
#                         friend_directory=friend_directory)),
#                 Process(
#                     "EWKWp",
#                     EWKWpEstimation(
#                         era,
#                         directory,
#                         channel,
#                         friend_directory=friend_directory)),
#                 Process(
#                     "EWKWm",
#                     EWKWmEstimation(
#                         era,
#                         directory,
#                         channel,
#                         friend_directory=friend_directory))
#             ])


class WTEstimation(WEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(WEstimation, self).__init__(
            name="W",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(Cut("gen_match_1==3||gen_match_1==4", "wt_genmatch"))


class WLEstimation(WEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(WEstimation, self).__init__(
            name="W",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        return Cuts(Cut("!(gen_match_1==3||gen_match_1==4)", "wl_genmatch"))


class TTEstimation(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(TTEstimation, self).__init__(
            name="TT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        # if self.channel.name=="em":
        #     return Weights(
        #         Weight("0.989*topPtReweightWeightRun1", "topPtReweightWeight"), #TODO topPTRun1 reweight ?
        #         Weight("numberGeneratedEventsWeight", "numberGeneratedEventsWeight"),
        #         Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
        #         Weight("triggerWeight_1*triggerWeight_2*identificationWeight_1*identificationWeight_2*trackWeight_1*trackWeight_2", "eventWeight"),
        #         Weight("puweight", "puweight"),
        #         Weight(
        #             "(((gen_match_1 == 5)*0.95 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.95 + (gen_match_2 != 5)))",
        #             "hadronic_tau_sf"),
        #         Weight("prefiringweight", "prefireWeight"),
        #         self.era.lumi_weight)
        # else:
        return Weights(
            Weight("topPtReweightWeight", "topPtReweightWeight"),
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
            Weight("generatorWeight", "generatorWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^TT$",
            "data": False,
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)


class TTLEstimation(TTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "et" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            emb_veto = "!(gen_match_1==5 && gen_match_2==5)"
            ff_veto = "!(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==4)"
            ff_veto = "(1.0)"
        elif "mm" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==4)"
            ff_veto = "(1.0)"
        return Cuts(Cut("%s && %s" % (emb_veto, ff_veto),
                        "tt_emb_and_ff_veto"))


class TTTEstimation(TTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            tt_cut = "gen_match_1==4 && gen_match_2==5"
        elif "et" in self.channel.name:
            tt_cut = "gen_match_1==3 && gen_match_2==5"
        elif "tt" in self.channel.name:
            tt_cut = "gen_match_1==5 && gen_match_2==5"
        elif "em" in self.channel.name:
            tt_cut = "gen_match_1==3 && gen_match_2==4"
        elif "mm" in self.channel.name:
            tt_cut = "gen_match_1==4 && gen_match_2==4"
        return Cuts(Cut(tt_cut, "ttt_cut"))


class TTJEstimation(TTEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(TTEstimation, self).__init__(
            name="TTJ",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "(gen_match_2 == 6 && gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            ct = "0.0 == 1.0"
        return Cuts(Cut(ct, "tt_fakes"))


# class VVEstimation(EstimationMethod):
#     def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
#             get_triggerweight_for_channel=get_triggerweight_for_channel,
#             get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
#             get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
#             get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
#         super(VVEstimation, self).__init__(
#             name="VV",
#             folder=folder,
#             get_triggerweight_for_channel=get_triggerweight_for_channel,
#             get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
#             get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
#             get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
#             era=era,
#             directory=directory,
#             friend_directory=friend_directory,
#             channel=channel,
#             mc_campaign="RunIISummer16MiniAODv3")

#     def get_weights(self):
#         # if self.channel.name=="em":
#         #     return Weights(
#         #         Weight("numberGeneratedEventsWeight", "numberGeneratedEventsWeight"),
#         #         Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
#         #         Weight("triggerWeight_1*identificationWeight_1*identificationWeight_2*trackWeight_1*trackWeight_2", "eventWeight"),
#         #         Weight("puweight", "puweight"),
#         #         self.get_tauByIsoIdWeight_for_channel(self.channel),
#         #         Weight("prefiringweight", "prefireWeight"),
#         #         self.era.lumi_weight)
#         # else:
#         return Weights(
#             Weight("isoWeight_1*isoWeight_2","isoWeight"),
#             Weight("idWeight_1*idWeight_2","idWeight"),
#             self.get_tauByIsoIdWeight_for_channel(self.channel),
#             Weight("puweight", "puweight"),
#             Weight("trackWeight_1*trackWeight_2","trackweight"),
#             self.get_triggerweight_for_channel(self.channel._name),
#             self.get_singlelepton_triggerweight_for_channel(self.channel.name),
#             Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
#             self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
#             get_eleRecoWeight_for_channel(self.channel.name),
#             Weight("prefiringweight", "prefireWeight"),
#             # MC weights
#             #TODO doing stitching with cross-section as reference for WW, WZ, ZZ, so WATCH OUT after changing the cross-sections!!!
#             Weight("1.252790591041545e-07*(abs(crossSectionPerEventWeight - 63.21) < 0.01) + 5.029933132068942e-07*(abs(crossSectionPerEventWeight - 10.32) < 0.01) + 2.501519047441559e-07*(abs(crossSectionPerEventWeight - 22.82) < 0.01) + numberGeneratedEventsWeight*(abs(crossSectionPerEventWeight - 63.21) > 0.01 && abs(crossSectionPerEventWeight - 10.32) > 0.01 && abs(crossSectionPerEventWeight - 22.82) > 0.01)","numberGeneratedEventsWeight"),
#             #TODO correct to proper cross-section values. WILL BE DEPRECATED after fixing cross-sections in datasets.json & producing ntuples
#             Weight("118.7*(abs(crossSectionPerEventWeight - 63.21) < 0.01) + crossSectionPerEventWeight*(abs(crossSectionPerEventWeight - 63.21) > 0.01)", "crossSectionPerEventWeight"),
#             Weight("generatorWeight", "generatorWeight"),
#             self.era.lumi_weight)

#     def get_files(self):
#         query = {
#             "process": "(WW|ZZ|WZ)$",  # Query for Di-Boson samples
#             "data": False,
#             "campaign": self._mc_campaign,
#             "generator": "pythia8"
#         }
#         files = self.era.datasets_helper.get_nicks_with_query(query)

#         query = {
#             "process":
#             "(STt-channelantitop4finclusiveDecays|STt-channeltop4finclusiveDecays|STtWantitop5finclusiveDecays|STtWtop5finclusiveDecays)",
#             "data":
#             False,
#             "campaign":
#             self._mc_campaign
#         }
#         files += self.era.datasets_helper.get_nicks_with_query(query)

#         log_query(self.name, "<optimzed out>", files)
#         return self.artus_file_names(files)


class VVEstimation(
        EstimationMethod
):  # TODO swap back to newer Estimation with Di-Boson Samples as soon as samples are available
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(VVEstimation, self).__init__(
            name="VV",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_weights(self):
        # if self.channel.name=="em":
        #     return Weights(
        #         Weight("numberGeneratedEventsWeight", "numberGeneratedEventsWeight"),
        #         Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),
        #         Weight("triggerWeight_1*identificationWeight_1*identificationWeight_2*trackWeight_1*trackWeight_2", "eventWeight"),
        #         Weight("puweight", "puweight"),
        #         self.get_tauByIsoIdWeight_for_channel(self.channel),
        #         Weight("prefiringweight", "prefireWeight"),
        #         self.era.lumi_weight)
        # else:
        return Weights(
            Weight("isoWeight_1*isoWeight_2", "isoWeight"),
            Weight("idWeight_1*idWeight_2", "idWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel),
            Weight("puweight", "puweight"),
            Weight("trackWeight_1*trackWeight_2", "trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight",
                   "leptonTauFakeRateWeight"),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            get_eleRecoWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),
            # MC weights
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight(
                "(1.0+0.56*(abs(crossSectionPerEventWeight-75.769996)<0.00001))",
                "VV_NNLO_reweight"),
            Weight("generatorWeight", "generatorWeight"),
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process":
            "(WWTo1L1Nu2Q|" + "WZTo1L1Nu2Q|" + "WZTo1L3Nu|" + "WZTo2L2Q|" +
            "ZZTo2L2Q" + ")",
            "data":
            False,
            "campaign":
            self._mc_campaign,
            "generator":
            "amcatnlo-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process": "(VVTo2L2Nu|ZZTo4L)",
            "extension": "ext1",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo-pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process": "WZJToLLLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process":
            "(STt-channelantitop4finclusiveDecays|STt-channeltop4finclusiveDecays|STtWantitop5finclusiveDecays|STtWtop5finclusiveDecays)",
            "data": False,
            "campaign": self._mc_campaign
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        log_query(self.name, "<optimzed out>", files)
        return self.artus_file_names(files)


class VVLEstimation(VVEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(VVEstimation, self).__init__(
            name="VVL",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "et" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==5)"
            ff_veto = "!(gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            emb_veto = "!(gen_match_1==5 && gen_match_2==5)"
            ff_veto = "!(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            emb_veto = "!(gen_match_1==3 && gen_match_2==4)"
            ff_veto = "(1.0)"
        elif "mm" in self.channel.name:
            emb_veto = "!(gen_match_1==4 && gen_match_2==4)"
            ff_veto = "(1.0)"
        return Cuts(Cut("%s && %s" % (emb_veto, ff_veto),
                        "vv_emb_and_ff_veto"))


class VVTEstimation(VVEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(VVEstimation, self).__init__(
            name="VVT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        if "mt" in self.channel.name:
            tt_cut = "gen_match_1==4 && gen_match_2==5"
        elif "et" in self.channel.name:
            tt_cut = "gen_match_1==3 && gen_match_2==5"
        elif "tt" in self.channel.name:
            tt_cut = "gen_match_1==5 && gen_match_2==5"
        elif "em" in self.channel.name:
            tt_cut = "gen_match_1==3 && gen_match_2==4"
        elif "mm" in self.channel.name:
            tt_cut = "gen_match_1==4 && gen_match_2==4"
        return Cuts(Cut(tt_cut, "vvt_cut"))


class VVJEstimation(VVEstimation):
    def __init__(
            self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(VVEstimation, self).__init__(
            name="VVJ",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIISummer16MiniAODv3")

    def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "(gen_match_2 == 6 && gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            ct = "0.0 == 1.0"
        return Cuts(Cut(ct, "vv_fakes"))


class QCDEstimation_SStoOS_MTETEM(SStoOSEstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            friend_directory=None,
            extrapolation_factor=1.0,
            folder="nominal",
            qcd_weight=Weight("1.0", "qcd_Weight"),
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(QCDEstimation_SStoOS_MTETEM, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            friend_directory=friend_directory,
            data_process=data_process,
            extrapolation_factor=extrapolation_factor,
            qcd_weight=qcd_weight)


class QCDEstimationTT(ABCDEstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(QCDEstimationTT, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            bg_processes=bg_processes,
            data_process=data_process,
            AC_cut_names=
            [  # cuts to be removed to include region for shape derivation
                "tau_2_iso"
            ],
            BD_cuts=
            [  # cuts to be applied to restrict to region for shape derivation
                Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5",
                    "tau_2_iso"),
                Cut("byLooseIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                    "tau_2_iso_loose")
            ],
            AB_cut_names=
            [  # cuts to be removed to include region for the determination of the extrapolation derivation
                "os"
            ],
            CD_cuts=
            [  # cuts to be applied to restrict to region for the determination of the extrapolation derivation
                Cut("q_1*q_2>0", "ss")
            ])


class WEstimationWithQCD(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            w_process,
            qcd_ss_to_os_extrapolation_factor,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(WEstimationWithQCD, self).__init__(
            name="WJets",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign=None)
        self._bg_processes = bg_processes
        self._data_process = data_process
        self._w_process = w_process
        self._qcd_ss_to_os_extrapolation_factor = qcd_ss_to_os_extrapolation_factor

    def create_root_objects(self, systematic):

        # create category for MC WJets shape estimation in the signal region
        signal_region = copy.deepcopy(systematic.category)
        signal_region.name = (signal_region.name +
                              "_for_wjets_mc").lstrip(self.channel.name + "_")

        # create control regions for W yield estimation
        high_mt_ss_control_region = copy.deepcopy(systematic.category)
        high_mt_ss_control_region.name = "wjets_high_mt_ss_cr"
        high_mt_ss_control_region._variable = None

        high_mt_ss_control_region.cuts.remove("m_t")
        high_mt_ss_control_region.cuts.remove("os")

        high_mt_ss_control_region.cuts.add(Cut("mt_1>70", "m_t"))
        high_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W high mt to low mt extrapolation factor
        high_mt_os_control_region = copy.deepcopy(
            systematic.category)  # this one also used for W yield estimation
        high_mt_os_control_region.name = "wjets_high_mt_os_cr"
        high_mt_os_control_region._variable = None

        high_mt_os_control_region.cuts.remove("m_t")

        high_mt_os_control_region.cuts.add(Cut("mt_1>70", "m_t"))

        low_mt_os_control_region = copy.deepcopy(systematic.category)
        low_mt_os_control_region.name = "wjets_low_mt_os_cr"
        low_mt_os_control_region._variable = None

        low_mt_ss_control_region = copy.deepcopy(systematic.category)
        low_mt_ss_control_region.name = "wjets_low_mt_ss_cr"
        low_mt_ss_control_region._variable = None

        low_mt_ss_control_region.cuts.remove("os")

        low_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W ss to os extrapolation factor
        inclusive_os_control_region = copy.deepcopy(systematic.category)
        inclusive_os_control_region.name = "wjets_os_cr"
        inclusive_os_control_region._variable = None

        inclusive_os_control_region.cuts.remove("m_t")

        inclusive_ss_control_region = copy.deepcopy(systematic.category)
        inclusive_ss_control_region.name = "wjets_ss_cr"
        inclusive_ss_control_region._variable = None

        inclusive_ss_control_region.cuts.remove("m_t")
        inclusive_ss_control_region.cuts.remove("os")

        inclusive_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # initialize root objects and systematics
        root_objects = []
        systematic._WandQCD_systematics = []

        # for extrapolation factors: only W MC is needed
        for category in [
                high_mt_os_control_region, low_mt_os_control_region,
                high_mt_ss_control_region, low_mt_ss_control_region,
                inclusive_os_control_region, inclusive_ss_control_region
        ]:
            s = Systematic(category=category,
                           process=self._w_process,
                           analysis=systematic.analysis,
                           era=self.era,
                           variation=systematic.variation,
                           mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        # for yields in high mt control regions: data and other bg processes needed
        for process in [self._data_process] + self._bg_processes:
            for category in [
                    high_mt_os_control_region, high_mt_ss_control_region
            ]:
                s = Systematic(category=category,
                               process=process,
                               analysis=systematic.analysis,
                               era=self.era,
                               variation=systematic.variation,
                               mass=125)
                systematic._WandQCD_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects

        # for signal region shape
        s = Systematic(category=signal_region,
                       process=self._w_process,
                       analysis=systematic.analysis,
                       era=self.era,
                       variation=systematic.variation,
                       mass=125)
        systematic._WandQCD_systematics.append(s)
        s.create_root_objects()
        root_objects += s.root_objects

        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_WandQCD_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _WandQCD_systematics needed for WandQCD estimation.",
                systematic.name)
            raise Exception

        # Sort shapes and counts
        wjets_mc_shape = None
        wjets_high_mt_ss_cr_counts = {}
        wjets_high_mt_os_cr_counts = {}
        wjets_low_mt_os_cr_count = None
        wjets_low_mt_ss_cr_count = None
        wjets_os_cr_count = None
        wjets_ss_cr_count = None
        for s in systematic._WandQCD_systematics:
            s.do_estimation()
            if s.category.name.endswith("for_wjets_mc"):
                wjets_mc_shape = s.shape
            elif s.category.name.endswith("wjets_high_mt_ss_cr"):
                wjets_high_mt_ss_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_os_cr"):
                wjets_high_mt_os_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_low_mt_os_cr"):
                wjets_low_mt_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_low_mt_ss_cr"):
                wjets_low_mt_ss_cr_count = s.shape
            elif s.category.name.endswith("wjets_os_cr"):
                wjets_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_ss_cr"):
                wjets_ss_cr_count = s.shape

        # Determine extrapolation factors
        R_ss_to_os = wjets_os_cr_count.result / wjets_ss_cr_count.result

        wjets_integral_low_mt_os = wjets_low_mt_os_cr_count.result
        wjets_integral_high_mt_os = wjets_high_mt_os_cr_counts.pop(
            self._w_process.name).result
        logger.debug("Integral of WJets MC in low mt OS region: %s",
                     str(wjets_integral_low_mt_os))
        logger.debug("Integral of WJets MC in high mt OS region: %s",
                     str(wjets_integral_high_mt_os))

        R_high_to_low_mt_os = wjets_integral_low_mt_os / wjets_integral_high_mt_os
        R_high_to_low_mt_ss = wjets_low_mt_ss_cr_count.result / wjets_high_mt_ss_cr_counts.pop(
            self._w_process.name).result
        logger.debug("WJets SS to OS extrapolation factor: %s",
                     str(R_ss_to_os))
        logger.debug("WJets high to low mt os extrapolation factor: %s",
                     str(R_high_to_low_mt_os))
        logger.debug("WJets high to low mt ss extrapolation factor: %s",
                     str(R_high_to_low_mt_ss))

        # Determine yields in wjets CRs
        logger.debug(
            "Data yield in ss high mt region: %s",
            str(wjets_high_mt_ss_cr_counts[self._data_process.name].result))
        high_mt_ss_yield = wjets_high_mt_ss_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_ss_cr_counts.values()])
        sum_mc = sum([s.result for s in wjets_high_mt_ss_cr_counts.values()])
        logger.debug("MC yield to be subtracted: %s", str(sum_mc))
        for name, s in wjets_high_mt_ss_cr_counts.items():
            logger.debug(name + " : " + str(s.result / sum_mc))

        logger.debug(
            "Data yield in os high mt region: %s",
            str(wjets_high_mt_os_cr_counts[self._data_process.name].result))
        high_mt_os_yield = wjets_high_mt_os_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_os_cr_counts.values()])
        sum_mc = sum([s.result for s in wjets_high_mt_os_cr_counts.values()])
        logger.debug("MC yield to be subtracted: %s", str(sum_mc))
        for name, s in wjets_high_mt_os_cr_counts.items():
            logger.debug(name + " : " + str(s.result / sum_mc))

        logger.debug("WJets + QCD yield in ss high mt region: %s",
                     str(high_mt_ss_yield))
        logger.debug("WJets + QCD yield in os high mt region: %s",
                     str(high_mt_os_yield))

        # Derive and normalize final shape
        logger.debug("WJets MC yield in signal region: %s",
                     str(wjets_integral_low_mt_os))
        sf = R_ss_to_os * (
            high_mt_os_yield -
            self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor
            ) / wjets_integral_high_mt_os
        estimated_yield = R_high_to_low_mt_os * R_ss_to_os * (
            high_mt_os_yield -
            self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor)
        logger.debug("WJets Estimated yield in signal region: %s",
                     str(estimated_yield))
        logger.debug("Scale WJets by %s", str(sf))
        wjets_shape = copy.deepcopy(wjets_mc_shape)
        wjets_shape.result.Scale(sf)

        # Rename root object accordingly
        wjets_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        wjets_shape.replace_negative_entries_and_renormalize(tolerance=100.5)

        return wjets_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError


class QCDEstimationWithW(EstimationMethod):
    def __init__(
            self,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            w_process,
            qcd_ss_to_os_extrapolation_factor,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel):
        super(QCDEstimationWithW, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=
            get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign=None)
        self._bg_processes = bg_processes
        self._data_process = data_process
        self._w_process = w_process
        self._qcd_ss_to_os_extrapolation_factor = qcd_ss_to_os_extrapolation_factor

    def create_root_objects(self, systematic):

        # create category for WJets and QCD shape estimation in the qcd control region
        qcd_control_region = copy.deepcopy(systematic.category)
        qcd_control_region.name = (qcd_control_region.name +
                                   "_ss_for_qcd").lstrip(self.channel.name +
                                                         "_")

        qcd_control_region.cuts.remove("os")

        qcd_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W yield estimation
        high_mt_ss_control_region = copy.deepcopy(systematic.category)
        high_mt_ss_control_region.name = "wjets_high_mt_ss_cr"
        high_mt_ss_control_region._variable = None

        high_mt_ss_control_region.cuts.remove("m_t")
        high_mt_ss_control_region.cuts.remove("os")

        high_mt_ss_control_region.cuts.add(Cut("mt_1>70", "m_t"))
        high_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W high mt to low mt extrapolation factor
        high_mt_os_control_region = copy.deepcopy(
            systematic.category)  # this one also used for W yield estimation
        high_mt_os_control_region.name = "wjets_high_mt_os_cr"
        high_mt_os_control_region._variable = None

        high_mt_os_control_region.cuts.remove("m_t")

        high_mt_os_control_region.cuts.add(Cut("mt_1>70", "m_t"))

        low_mt_os_control_region = copy.deepcopy(systematic.category)
        low_mt_os_control_region.name = "wjets_low_mt_os_cr"
        low_mt_os_control_region._variable = None

        low_mt_ss_control_region = copy.deepcopy(systematic.category)
        low_mt_ss_control_region.name = "wjets_low_mt_ss_cr"
        low_mt_ss_control_region._variable = None

        low_mt_ss_control_region.cuts.remove("os")

        low_mt_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # create control regions for W ss to os extrapolation factor
        inclusive_os_control_region = copy.deepcopy(systematic.category)
        inclusive_os_control_region.name = "wjets_os_cr"
        inclusive_os_control_region._variable = None

        inclusive_os_control_region.cuts.remove("m_t")

        inclusive_ss_control_region = copy.deepcopy(systematic.category)
        inclusive_ss_control_region.name = "wjets_ss_cr"
        inclusive_ss_control_region._variable = None

        inclusive_ss_control_region.cuts.remove("m_t")
        inclusive_ss_control_region.cuts.remove("os")

        inclusive_ss_control_region.cuts.add(Cut("q_1*q_2>0", "ss"))

        # initialize root objects and systematics
        root_objects = []
        systematic._WandQCD_systematics = []

        # for extrapolation factors: only W MC is needed
        for category in [
                high_mt_os_control_region, low_mt_os_control_region,
                high_mt_ss_control_region, low_mt_ss_control_region,
                inclusive_os_control_region, inclusive_ss_control_region
        ]:
            s = Systematic(category=category,
                           process=self._w_process,
                           analysis=systematic.analysis,
                           era=self.era,
                           variation=systematic.variation,
                           mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        # for yields in high mt control regions: data and other bg processes needed
        for process in [self._data_process] + self._bg_processes:
            for category in [
                    high_mt_os_control_region, high_mt_ss_control_region
            ]:
                s = Systematic(category=category,
                               process=process,
                               analysis=systematic.analysis,
                               era=self.era,
                               variation=systematic.variation,
                               mass=125)
                systematic._WandQCD_systematics.append(s)
                s.create_root_objects()
                root_objects += s.root_objects

        # for Wjets and QCD shape
        for process in [self._data_process, self._w_process
                        ] + self._bg_processes:
            s = Systematic(category=qcd_control_region,
                           process=process,
                           analysis=systematic.analysis,
                           era=self.era,
                           variation=systematic.variation,
                           mass=125)
            systematic._WandQCD_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects

        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_WandQCD_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _WandQCD_systematics needed for WandQCD estimation.",
                systematic.name)
            raise Exception

        # Sort shapes and counts
        qcd_control_region_shapes = {}
        wjets_high_mt_ss_cr_counts = {}
        wjets_high_mt_os_cr_counts = {}
        wjets_low_mt_os_cr_count = None
        wjets_low_mt_ss_cr_count = None
        wjets_os_cr_count = None
        wjets_ss_cr_count = None
        for s in systematic._WandQCD_systematics:
            s.do_estimation()
            if s.category.name.endswith("ss_for_qcd"):
                qcd_control_region_shapes[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_ss_cr"):
                wjets_high_mt_ss_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_high_mt_os_cr"):
                wjets_high_mt_os_cr_counts[s.process.name] = s.shape
            elif s.category.name.endswith("wjets_low_mt_os_cr"):
                wjets_low_mt_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_low_mt_ss_cr"):
                wjets_low_mt_ss_cr_count = s.shape
            elif s.category.name.endswith("wjets_os_cr"):
                wjets_os_cr_count = s.shape
            elif s.category.name.endswith("wjets_ss_cr"):
                wjets_ss_cr_count = s.shape

        # Determine extrapolation factors
        R_ss_to_os = wjets_os_cr_count.result / wjets_ss_cr_count.result

        wjets_integral_low_mt_ss = wjets_low_mt_ss_cr_count.result
        wjets_integral_high_mt_ss = wjets_high_mt_ss_cr_counts.pop(
            self._w_process.name).result

        R_high_to_low_mt_os = wjets_low_mt_os_cr_count.result / wjets_high_mt_os_cr_counts.pop(
            self._w_process.name).result
        R_high_to_low_mt_ss = wjets_integral_low_mt_ss / wjets_integral_high_mt_ss

        # Determine yields in wjets CRs
        high_mt_ss_yield = wjets_high_mt_ss_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_ss_cr_counts.values()])

        high_mt_os_yield = wjets_high_mt_os_cr_counts.pop(
            self._data_process.name).result - sum(
                [s.result for s in wjets_high_mt_os_cr_counts.values()])

        # Derive and normalize final shape for QCD
        wjets_shape = qcd_control_region_shapes.pop(self._w_process.name)
        logger.debug("WJets MC yield in qcd control region: %s",
                     str(wjets_integral_low_mt_ss))
        sf = (high_mt_os_yield -
              self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                  R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor
              ) / wjets_integral_high_mt_ss
        estimated_yield = R_high_to_low_mt_ss * (
            high_mt_os_yield -
            self._qcd_ss_to_os_extrapolation_factor * high_mt_ss_yield) / (
                R_ss_to_os - self._qcd_ss_to_os_extrapolation_factor)
        logger.debug("WJets Estimated yield in qcd control region: %s",
                     str(estimated_yield))
        logger.debug("Scale WJets by %s", str(sf))
        wjets_shape.result.Scale(sf)
        wjets_shape._result.Write()

        qcd_shape = copy.deepcopy(
            qcd_control_region_shapes.pop(self._data_process.name))
        qcd_shape.result.Add(wjets_shape.result, -1.0)
        for sh in qcd_control_region_shapes.values():
            qcd_shape.result.Add(sh.result, -1.0)
        # Saving QCD shape in ss control region
        qcd_ss_shape = copy.deepcopy(qcd_shape)
        ss_category_name = ""
        for s in systematic._WandQCD_systematics:
            if s.category.name.endswith("ss_for_qcd"):
                ss_category_name = s.category._name
        qcd_ss_shape.name = systematic.name.replace(systematic.category._name,
                                                    ss_category_name)
        qcd_ss_shape._result.Write()

        # Rescale QCD shape for signal region
        qcd_shape.result.Scale(self._qcd_ss_to_os_extrapolation_factor)

        # Rename root object accordingly
        qcd_shape.name = systematic.name

        # Replace negative entries by zeros and renormalize shape
        qcd_shape.replace_negative_entries_and_renormalize(tolerance=100.5)

        return qcd_shape

    # Data-driven estimation, no associated files and weights
    def get_files(self):
        raise NotImplementedError

    def get_weights(self):
        raise NotImplementedError
