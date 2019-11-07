

# -*- coding: utf-8 -*-

from cutstring import *
from estimation_methods import EstimationMethod, SStoOSEstimationMethod, ABCDEstimationMethod, SumUpEstimationMethod, NewFakeEstimationMethodLT, NewFakeEstimationMethodTT
from estimation_methods_2016 import DataEstimation as DataEstimation2016
from estimation_methods_2016 import WEstimationWithQCD as WEstimationWithQCD2016
from estimation_methods_2016 import QCDEstimationWithW as QCDEstimationWithW2016
from estimation_methods_2016 import ggH_htxs, qqH_htxs
from systematics import *
from era import log_query
from process import *


def get_triggerweight_for_channel(channel):
    weight = Weight("1.0", "triggerweight")

    singleMC = "singleTriggerMCEfficiencyWeightKIT_1"
    crossMCL = "crossTriggerMCEfficiencyWeight_1"
    MCTau_1 = "((byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5 && byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)*crossTriggerMCEfficiencyWeight_vloose_MVAv2_1 + (byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)*crossTriggerMCEfficiencyWeight_tight_MVAv2_1)"
    MCTau_2 = MCTau_1.replace("_1","_2")

    if "mt" in channel:
        trig_sL = "(trg_singlemuon_27 || trg_singlemuon_24)"
        trig_X = "(pt_1 > 21 && pt_1 < 25 && trg_crossmuon_mu20tau27)"

        # Eff = Eff(singleL)*(1 - Eff(xTau)) + Eff(xL)*Eff(xTau)
        #MuTauMC = "*".join([trig_sL,singleMC,"(1-"+trig_X+"*"+crossMCL+")"])+"+"+"*".join([trig_X,crossMCL,MCTau_2])
        #MuTauData = MuTauMC.replace("MC","Data")
        #MuTau = "("+MuTauData+")/("+MuTauMC+")"

        MuTauMC = "*".join([trig_sL, singleMC]) + "+" + "*".join([trig_X, crossMCL, MCTau_2])
        MuTauData = MuTauMC.replace("MC","Data")
        MuTau = "("+MuTauData+")/("+MuTauMC+")"
        weight = Weight(MuTau,"triggerweight")

    elif "et" in channel:
        trig_sL = "(trg_singleelectron_35 || trg_singleelectron_32 || trg_singleelectron_27)"
        trig_X = "(pt_1>25 && pt_1<28 && trg_crossele_ele24tau30)"

        # Eff = Eff(singleL)*(1 - Eff(xTau)) + Eff(xL)*Eff(xTau)
        #ElTauMC = "*".join([trig_sL,singleMC,"(1-"+trig_X+"*"+crossMCL+")"])+"+"+"*".join([trig_X,crossMCL,MCTau_2])
        #ElTauData = ElTauMC.replace("MC","Data")
        #ElTau = "("+ElTauData+")/("+ElTauMC+")"

        ElTauMC = "*".join([trig_sL, singleMC]) + "+" + "*".join([trig_X, crossMCL, MCTau_2])
        ElTauData = ElTauMC.replace("MC","Data")
        ElTau = "("+ElTauData+")/("+ElTauMC+")"
        weight = Weight(ElTau,"triggerweight")

    elif "tt" in channel:
        DiTauMC = "*".join([MCTau_1,MCTau_2])
        DiTauData = DiTauMC.replace("MC","Data")
        DiTau = "("+DiTauData+")/("+DiTauMC+")"
        weight = Weight(DiTau,"triggerweight")

    elif "em" in channel:
        weight = Weight(
            "(trigger_23_data_Weight_2*trigger_12_data_Weight_1*(trg_muonelectron_mu23ele12==1)+trigger_23_data_Weight_1*trigger_8_data_Weight_2*(trg_muonelectron_mu8ele23==1) - trigger_23_data_Weight_2*trigger_23_data_Weight_1*(trg_muonelectron_mu8ele23==1 && trg_muonelectron_mu23ele12==1))/(trigger_23_mc_Weight_2*trigger_12_mc_Weight_1*(trg_muonelectron_mu23ele12==1)+trigger_23_mc_Weight_1*trigger_8_mc_Weight_2*(trg_muonelectron_mu8ele23==1) - trigger_23_mc_Weight_2*trigger_23_mc_Weight_1*(trg_muonelectron_mu8ele23==1 && trg_muonelectron_mu23ele12==1))",
            "trigger_lepton_sf")

    elif "mm" in channel:
        weight = Weight(
            "singleTriggerDataEfficiencyWeightKIT_1/singleTriggerMCEfficiencyWeightKIT_1",
            "trigger_lepton_sf")

    return weight


def get_singlelepton_triggerweight_for_channel(channel):
    weight = Weight("1.0","triggerweight")

    MCTau_1 = "((byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5 && byMediumIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)*crossTriggerMCEfficiencyWeight_medium_MVAv2_1 + (byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5)*crossTriggerMCEfficiencyWeight_tight_MVAv2_1)"
    MCTau_2 = MCTau_1.replace("_1","_2")

    if "mt" in channel or "et" in channel:
        weight = Weight("singleTriggerDataEfficiencyWeightKIT_1/singleTriggerMCEfficiencyWeightKIT_1","triggerweight")
    elif "tt" in channel:
        DiTauMC = "*".join([MCTau_1,MCTau_2])
        DiTauData = DiTauMC.replace("MC","Data")
        DiTau = "("+DiTauData+")/("+DiTauMC+")"
        weight = Weight(DiTau,"triggerweight")

    return weight


def get_tauByIsoIdWeight_for_channel(channel):
    # WPs: VLoose 0.88, Loose 0.89, Medium 0.89, Tight 0.89, VTight 0.86, VVTight 0.84. Currently used: SR mt,et Tight; SR tt Tight, anti-iso CR tt Medium; VVLoose is used for SF estimation and therefore not listed here.
    # Source: https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
    weight = Weight("1.0","taubyIsoIdWeight")
    if "mt" in channel or "et" in channel:
        weight = Weight("((gen_match_2 == 5)*0.89 + (gen_match_2 != 5))", "taubyIsoIdWeight")
    elif "tt" in channel:
        weight = Weight("((gen_match_1 == 5)*0.89 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.89 + (gen_match_2 != 5))", "taubyIsoIdWeight")
    return weight


def get_eleHLTZvtxWeight_for_channel(channel):
    weight = Weight("1.0","eleHLTZvtxWeight")
    if "et" in channel:
        weight = Weight("(trg_singleelectron_35 || trg_singleelectron_32 || trg_singleelectron_27 || trg_crossele_ele24tau30)*0.991 + (!(trg_singleelectron_35 || trg_singleelectron_32 || trg_singleelectron_27 || trg_crossele_ele24tau30))*1.0", "eleHLTZvtxWeight")
    return weight


class DataEstimation(DataEstimation2016):
    pass


class WEstimationWithQCD(WEstimationWithQCD2016):
    pass


class QCDEstimationWithW(QCDEstimationWithW2016):
    pass


class QCDEstimation_SStoOS_MTETEM(SStoOSEstimationMethod):
    def __init__(self,
            era,
            directory,
            channel,
            bg_processes,
            data_process,
            friend_directory=None,
            extrapolation_factor=1.0,
            qcd_weight=Weight("1.0","qcd_Weight"),
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(QCDEstimation_SStoOS_MTETEM, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            friend_directory=friend_directory,
            data_process=data_process,
            extrapolation_factor=extrapolation_factor,
            qcd_weight = qcd_weight
            )


class QCDEstimation_ABCD_TT_ISO2(ABCDEstimationMethod):
    def __init__(self,
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
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(QCDEstimation_ABCD_TT_ISO2, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            friend_directory=friend_directory,
            data_process=data_process,
            AC_cut_names=[ # cuts applied in AC, which should be removed in the BD control regions
                "tau_2_iso",
            ],
            BD_cuts=[      # cuts to be applied instead of cuts removed above
                Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5", "tau_2_iso"),
                Cut("byMediumIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                    "tau_2_iso_loose"),
            ],
            AB_cut_names=[ # cuts applied in AB, which should be removed in the CD control regions
                "os"
            ],
            CD_cuts=[      # cuts to be applied instead of cuts removed above
                Cut("q_1*q_2>0", "ss")
            ]
        )


class QCDEstimation_ABCD_TT_ISO2_TRANSPOSED(ABCDEstimationMethod):
    def __init__(self,
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
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(QCDEstimation_ABCD_TT_ISO2_TRANSPOSED, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            friend_directory=friend_directory,
            data_process=data_process,
            AB_cut_names=[ # cuts applied in AB, which should be removed in the CD control regions
                "tau_2_iso"
            ],
            CD_cuts=[      # cuts to be applied instead of cuts removed above
                Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5", "tau_2_iso"),
                Cut("byMediumIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                    "tau_2_iso_loose"),
            ],
            AC_cut_names=[ # cuts applied in AC, which should be removed in the BD control regions
                "os"
            ],
            BD_cuts=[      # cuts to be applied instead of cuts removed above
                Cut("q_1*q_2>0", "ss")
            ]
        )


class QCDEstimation_ABCD_TT_ISO1(ABCDEstimationMethod):
    def __init__(self,
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
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(QCDEstimation_ABCD_TT_ISO1, self).__init__(
            name="QCD",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            bg_processes=bg_processes,
            friend_directory=friend_directory,
            data_process=data_process,
            AC_cut_names=[ # cuts applied in AC, which should be removed in the BD control regions
                "tau_1_iso"
            ],
            BD_cuts=[      # cuts to be applied instead of cuts removed above
                Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5", "tau_1_iso"),
                Cut("byLooseIsolationMVArun2v1DBoldDMwLT_1>0.5",
                    "tau_1_iso_loose")
            ],
            AB_cut_names=[ # cuts applied in AB, which should be removed in the CD control regions
                "os"
            ],
            CD_cuts=[      # cuts to be applied instead of cuts removed above
                Cut("q_1*q_2>0", "ss")
            ]
        )


class VVEstimation(EstimationMethod):
    def __init__(self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(VVEstimation, self).__init__(
            name="VV",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel.name),
            # self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("118.7*(abs(crossSectionPerEventWeight - 63.21) < 0.01) + crossSectionPerEventWeight*(abs(crossSectionPerEventWeight - 63.21) > 0.01)", "crossSectionPerEventWeight"),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(WW|ZZ|WZ)$",  # Query for Di-Boson samples
            "data": False,
            "generator": "^pythia8",
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process": "ST",  # Query for Single-Top samples (newer v2)
            "data": False,
            "scenario": "^PU2017$",
            "version": "v2",
            "generator": "powheg\-pythia8",
            "campaign": self._mc_campaign
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        query = {
            "process": "ST",  # Query for Single-Top samples (newer pileup mixing)
            "data": False,
            "scenario": "^PU2017newpmx$",
            "generator": "powheg\-pythia8",
            "campaign": self._mc_campaign
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class VVLEstimation(VVEstimation):
    def __init__(self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(VVEstimation, self).__init__(
            name="VVL",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

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
        return Cuts(Cut("%s && %s"%(emb_veto,ff_veto), "vv_emb_and_ff_veto"))

class VVTEstimation(VVEstimation):
    def __init__(self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(VVEstimation, self).__init__(
            name="VVT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

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
    def __init__(self,
            era,
            directory,
            channel,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(VVEstimation, self).__init__(
            name="VVJ",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "(gen_match_2 == 6 && gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            ct = "0.0 == 1.0"
        return Cuts(Cut(ct, "vv_fakes"))

class EWKZEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(EWKZEstimation, self).__init__(
            name="EWKZ",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel.name),
            # self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^EWKZ2Jets.",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class DYJetsToLLEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal", atNLO=False,
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
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        if self.atNLO:
            z_stitching_weight = Weight("((genbosonmass >= 50.0) * 2.9688e-05 + (genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight)","z_stitching_weight") # xsec_NNLO [pb] = 2075.14*3, N_inclusive_NLO = 209694262, xsec_NNLO/N_inclusive_NLO = 2.9688e-05; fraction of negative events in 'generatorWeight'
        else:
            z_stitching_weight = Weight("((genbosonmass >= 50.0)*6.3654e-05*((npartons == 0 || npartons >= 5)*1.0 + (npartons == 1)*0.1743 + (npartons == 2)*0.3556 + (npartons == 3)*0.2273 + (npartons == 4)*0.2104) + (genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight)","z_stitching_weight")
              # xsec_NNLO [pb] = 2075.14*3, N_inclusive = 97800939, xsec_NNLO/N_inclusive = 6.3654e-05 [pb] weights: [1.0, 0.1743347690195873, 0.3556373947627093, 0.22728901609456784, 0.21040417678899315]
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            z_stitching_weight,

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("zPtReweightWeight", "zPtReweightWeight"),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        queryM10 = {
            "process": "DYJetsToLL_M10to50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "ext1",
            "version": "v2"
        }
        queryM50_inclusive_2_3jet = {
            "process": "DY(|2|3)JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "version": "v1"
        }
        queryM50_1jet_v1 = {
            "process": "DY1JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "^$",
            "version": "v1"
        }
        queryM50_1jet_ext1_v2 = {
            "process": "DY1JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8",
            "extension": "ext1",
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
            files = self.era.datasets_helper.get_nicks_with_query(queryM50_inclusive_2_3jet) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM50_1jet_v1) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM50_1jet_ext1_v2) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM50_4jet) + \
                    self.era.datasets_helper.get_nicks_with_query(queryM10) + \
                    self.era.datasets_helper.get_nicks_with_query(queryEWKZ)
            log_query(self.name, queryM10, files)
        return self.artus_file_names(files)


class ZTTEstimation(DYJetsToLLEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal", atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        self.atNLO = atNLO
        name = "ZTTNLO" if self.atNLO else "ZTT"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

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
        return Cuts(Cut(tt_cut, "ztt_cut"))

class ZJEstimation(DYJetsToLLEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal", atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        self.atNLO = atNLO
        name = "ZJNLO" if self.atNLO else "ZJ"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "gen_match_2 == 6"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            ct = "0 == 1"
        return Cuts(Cut(ct, "dy_fakes"))


class ZLEstimation(DYJetsToLLEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal", atNLO=False,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        self.atNLO = atNLO
        name = "ZLNLO" if self.atNLO else "ZL"
        super(DYJetsToLLEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

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
        return Cuts(Cut("%s && %s"%(emb_veto,ff_veto), "dy_emb_and_ff_veto"))


class ZTTEmbeddedEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(ZTTEmbeddedEstimation, self).__init__(
            name="EMB",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            friend_directory=friend_directory,
            directory=directory,
            channel=channel,
            mc_campaign=None)

    def get_weights(self):
        if self.channel.name in ["mt"]:
            return Weights(
                Weight("generatorWeight",
                       "simulation_sf"),
                Weight("muonEffTrgWeight*muonEffIDWeight_1*muonEffIDWeight_2", "scale_factor"),
                Weight("idWeight_1*(trigger_24_27_Weight_1*(pt_1>25)+((0.81*(pt_1>=21 && pt_1<22) + 0.82*(pt_1>=22 && pt_1<23) + 0.83*(pt_1>=23))*(pt_1<25)))*isoWeight_1", "lepton_sf"),
                Weight("(pt_1>25)+(pt_1 >= 21 && pt_1<25)*((pt_2>=20 && pt_2<25)*0.12714+(pt_2>=25 && pt_2<30)*0.46930+0.71983*(pt_2>=30 && pt_2<35) + 0.75209*(pt_2>=35 && pt_2<40) + 0.78164*(pt_2>=40 && pt_2<45) + 0.83241*(pt_2>=45 && pt_2<50) + 0.86694*(pt_2>=50 && pt_2<60) + 0.89966*(pt_2>=60 && pt_2<80) + 0.88534*(pt_2>=80 && pt_2<100) + 0.90095*(pt_2>=100 && pt_2<150) + 0.84402*(pt_2>=150 && pt_2<200) + (pt_2>=200))","tau_leg_weight"),
                Weight("(gen_match_2==5)*0.97+(gen_match_2!=5)", "emb_tau_id"),
                Weight("gen_match_1==4 && gen_match_2==5","emb_veto"),
                Weight("embeddedDecayModeWeight", "decayMode_SF"))
        elif self.channel.name in ["et"]:
            return Weights(
                Weight("generatorWeight",
                       "simulation_sf"),
                Weight("muonEffTrgWeight*muonEffIDWeight_1*muonEffIDWeight_2", "scale_factor"),
                Weight("(pt_1>28)+(pt_1<28)*(crossTriggerDataEfficiencyWeight_tight_MVAv2_2*(abs(eta_1)>=1.5)+((1.29079*(pt_2>=30 && pt_2<35) + 1.06504*(pt_2>=35 && pt_2<40) + 0.93972*(pt_2>=40 && pt_2<45) + 0.91923*(pt_2>=45 && pt_2<50) + 0.89598*(pt_2>=50 && pt_2<60) + 0.90597*(pt_2>=60 && pt_2<80) + 0.88761*(pt_2>=80 && pt_2<100) + 0.90210*(pt_2>=100 && pt_2<150) + 0.84939*(pt_2>=150 && pt_2<200) + (pt_2>=200))*(abs(eta_1)<1.5)))","tau_leg_weight"),
                Weight("(pt_1>28)+(pt_1<28)*(crossTriggerDataEfficiencyWeight_1*(abs(eta_1)>=1.5)+((0.39*(pt_1>=25 && pt_1<26) + 0.46*(pt_1>=26 && pt_1<27) + 0.48*(pt_1>=27 && pt_1<28))*(abs(eta_1)<1.5)))","lepton_leg_weight"),
                Weight("idWeight_1*((pt_1>28)*(trigger_27_32_35_Weight_1*(abs(eta_1) < 1.5) + singleTriggerDataEfficiencyWeightKIT_1*(abs(eta_1)>=1.5))+(pt_1<28))*isoWeight_1", "lepton_sf"),
                Weight("(gen_match_2==5)*0.97+(gen_match_2!=5)", "emb_tau_id"),
                Weight("gen_match_1==3 && gen_match_2==5","emb_veto"),
                Weight("embeddedDecayModeWeight", "decayMode_SF"))
        elif self.channel.name == "tt":
            return Weights(
                Weight("generatorWeight",
                       "simulation_sf"),
                Weight("muonEffTrgWeight*muonEffIDWeight_1*muonEffIDWeight_2", "scale_factor"),
                Weight("(0.18321*(pt_1>=30 && pt_1<35) + 0.53906*(pt_1>=35 && pt_1<40) + 0.63658*(pt_1>=40 && pt_1<45) + 0.73152*(pt_1>=45 && pt_1<50) + 0.79002*(pt_1>=50 && pt_1<60) + 0.84666*(pt_1>=60 && pt_1<80) + 0.84919*(pt_1>=80 && pt_1<100) + 0.86819*(pt_1>=100 && pt_1<150) + 0.88206*(pt_1>=150 && pt_1<200) + (pt_1>=200))","tau1_leg_weight"),
                Weight("(0.18321*(pt_2>=30 && pt_2<35) + 0.53906*(pt_2>=35 && pt_2<40) + 0.63658*(pt_2>=40 && pt_2<45) + 0.73152*(pt_2>=45 && pt_2<50) + 0.79002*(pt_2>=50 && pt_2<60) + 0.84666*(pt_2>=60 && pt_2<80) + 0.84919*(pt_2>=80 && pt_2<100) + 0.86819*(pt_2>=100 && pt_2<150) + 0.88206*(pt_2>=150 && pt_2<200) + (pt_2>=200))","tau2_leg_weight"),
                Weight("((gen_match_1==5)*0.97+(gen_match_1!=5))*((gen_match_2==5)*0.97+(gen_match_2!=5))", "emb_tau_id"),
                Weight("gen_match_1==5 && gen_match_2==5","emb_veto"),
                Weight("embeddedDecayModeWeight", "decayMode_SF"))
        elif self.channel.name == "em":
            return Weights(
                Weight("1.043*generatorWeight", "simulation_sf"), # 1.043 for event loss in ntuple production, to be removed with new ntuples
                Weight("(gen_match_1==3 && gen_match_2==4)", "emb_gen_match"),
                Weight("muonEffTrgWeight*muonEffIDWeight_1*muonEffIDWeight_2", "scale_factor"),
                Weight("0.99*trackWeight_1*trackWeight_2*idWeight_1*isoWeight_1*idWeight_2*looseIsoWeight_2",
                       "idiso_lepton_sf"),
                Weight("(trigger_23_data_Weight_2*trigger_12_data_Weight_1*(trg_muonelectron_mu23ele12==1)+trigger_23_data_Weight_1*trigger_8_data_Weight_2*(trg_muonelectron_mu8ele23==1) - trigger_23_data_Weight_2*trigger_23_data_Weight_1*(trg_muonelectron_mu8ele23==1 && trg_muonelectron_mu23ele12==1))/(trigger_23_embed_Weight_2*trigger_12_embed_Weight_1*(trg_muonelectron_mu23ele12==1)+trigger_23_embed_Weight_1*trigger_8_embed_Weight_2*(trg_muonelectron_mu8ele23==1) - trigger_23_embed_Weight_2*trigger_23_embed_Weight_1*(trg_muonelectron_mu8ele23==1 && trg_muonelectron_mu23ele12==1))",
                       "trigger_lepton_sf"))


    def get_files(self):
        query = {"process": "Embedding2017(B|C|D|E|F)", "embedded": True}
        if self.channel.name == "mt":
            query["campaign"] = "MuTauFinalState"
            query["scenario"] = ".*v2"
        elif self.channel.name == "et":
            query["campaign"] = "ElTauFinalState"
            query["scenario"] = ".*v2"
        elif self.channel.name == "tt":
            query["campaign"] = "TauTauFinalState"
            query["scenario"] = ".*(v2|v3)"
        elif self.channel.name == "em":
            query["campaign"] = "ElMuFinalState"
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

#    def get_cuts(self):
#        ztt_genmatch_cut = Cut("1 == 1", "ztt_genmatch")
#        if self.channel.name in ["mt", "et"]:
#            ztt_genmatch_cut = Cut("gen_match_2==5", "ztt_genmatch")
#        elif self.channel.name == "tt":
#            ztt_genmatch_cut = Cut("(gen_match_1==5) && (gen_match_2==5)",
#                                   "ztt_genmatch")
#        elif self.channel.name == "em":
#            ztt_genmatch_cut = Cut("(gen_match_1>2) && (gen_match_2>3)",
#                                   "ztt_genmatch")
#        return Cuts(ztt_genmatch_cut)
    def get_cuts(self):
        return Cuts(Cut("((gen_match_1>2 && gen_match_1<6) &&  (gen_match_2>2 && gen_match_2<6))", "dy_genuine_tau"))


class ZttEmbeddingEstimation_ScaledToMC(EstimationMethod):
    def __init__(self, era, directory, channel, embedding_process,
                 ttbar_tautau_mc_process, z_tautau_mc_process, folder='nominal'):
        super(ZttEmbeddingEstimation_ScaledToMC, self).__init__(
            name="Ztt",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            mc_campaign=None)
        self._embedding_process = copy.deepcopy(embedding_process)
        self._ttbar_tautau_mc_process = copy.deepcopy(ttbar_tautau_mc_process)
        self._z_tautau_mc_process = copy.deepcopy(z_tautau_mc_process)

    def create_root_objects(self, systematic):
        yield_category = copy.deepcopy(systematic.category)
        yield_category._variable = None

        shape_category = copy.deepcopy(systematic.category)
        shape_category._name += "_unscaled"

        root_objects = []
        systematic._embedding_systematics = []
        shape_systematic = Systematic(
            category=shape_category,
            process=self._embedding_process,
            analysis=systematic.analysis,
            era=self.era,
            variation=systematic.variation,
            mass=125)
        systematic._embedding_systematics.append(shape_systematic)
        shape_systematic.create_root_objects()
        root_objects += shape_systematic.root_objects

        for process in [
                self._embedding_process, self._ttbar_tautau_mc_process,
                self._z_tautau_mc_process
        ]:
            s = Systematic(
                category=yield_category,
                process=process,
                analysis=systematic.analysis,
                era=self.era,
                variation=systematic.variation,
                mass=125)
            systematic._embedding_systematics.append(s)
            s.create_root_objects()
            root_objects += s.root_objects
        return root_objects

    def do_estimation(self, systematic):
        if not hasattr(systematic, "_embedding_systematics"):
            logger.fatal(
                "Systematic %s does not have attribute _embedding_systematics needed for embedding scaled to MC estimation.",
                systematic.name)
            raise Exception

        for s in systematic._embedding_systematics:
            s.do_estimation()

        shapes = [s.shape for s in systematic._embedding_systematics]

        # embedding shape
        embedding_shape = shapes[0]

        # scale factor = MC(TTT + ZTT) yield / embedding yield
        sf = (shapes[2].result + shapes[3].result) / shapes[1].result
        print "Scale factor", sf

        # scaling shape
        embedding_shape.result.Scale(sf)

        embedding_shape.name = systematic.name
        return embedding_shape


class WEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(WEstimation, self).__init__(
            name="W",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            #Weight("numberGeneratedEventsWeight","numberGeneratedEventsWeight"), # to be used only for one inclusive sample
            #Weight("crossSectionPerEventWeight","crossSectionPerEventWeight"), # to be used only for one inclusive sample
            Weight("((0.0010062794683*((npartons <= 0 || npartons >= 5)*1.0 + (npartons == 1)*0.1448 + (npartons == 2)*0.0887 + (npartons == 3)*0.0540 + (npartons == 4)*0.0535)) * (genbosonmass>=0.0) + numberGeneratedEventsWeight * crossSectionPerEventWeight * (genbosonmass<0.0))",
                "wj_stitching_weight"), # xsec_NNLO [pb] = 61526.7, N_inclusive = 61142756, xsec_NNLO/N_inclusive = 0.0010062794683314568 [pb] weights: [1.0, 0.14475938461191984, 0.08870222431187946, 0.05403112360513076, 0.05350733471645885]

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            # self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "W.?JetsToLNu|WGToLNuG",
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


class TTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(TTEstimation, self).__init__(
            name="TT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("(abs(crossSectionPerEventWeight - 380.1) < 0.1)*377.96 + (abs(crossSectionPerEventWeight - 87.31) < 0.1)*88.29 + (abs(crossSectionPerEventWeight - 364.4) < 0.1)*365.35", "crossSectionPerEventWeight"),
            #Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            Weight("topPtReweightWeight", "topPtReweightWeight"),
            self.get_triggerweight_for_channel(self.channel._name),
            # self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "TTTo.*",
            "scenario": "PU2017",
            "dbs" : ".*new_pmx.*",
            "data": False,
            "campaign": self._mc_campaign,
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class TTLEstimation(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(TTEstimation, self).__init__(
            name="TTL",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

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
        return Cuts(Cut("%s && %s"%(emb_veto,ff_veto), "tt_emb_and_ff_veto"))


class TTTEstimation(TTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(TTEstimation, self).__init__(
            name="TTT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

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
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(TTEstimation, self).__init__(
            name="TTJ",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_cuts(self):
        ct = ""
        if "mt" in self.channel.name or "et" in self.channel.name:
            ct = "(gen_match_2 == 6 && gen_match_2 == 6)"
        elif "tt" in self.channel.name:
            ct = "(gen_match_1 == 6 || gen_match_2 == 6)"
        elif "em" in self.channel.name:
            ct = "0 == 1"
        return Cuts(Cut(ct, "tt_fakes"))


class HTTEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="HTT",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            # self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(VBF|GluGlu|Z|W).*HToTauTau_M125",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class HWWEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HWWEstimation, self).__init__(
            name="HWW",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("0.0857883*(abs(numberGeneratedEventsWeight - 2e-06) < 1e-07) + 1.1019558*(abs(numberGeneratedEventsWeight - 2e-06) >= 1e-07)", "crossSectionPerEventWeight"), # TODO intermediate solution for cross-section weight, since wrong in datasets.json

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "(VBF|GluGlu).*HToWWTo2L2Nu_M125",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class ggHWWEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(ggHWWEstimation, self).__init__(
            name="ggHWW125",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("1.1019558", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "GluGlu.*HToWWTo2L2Nu_M125",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class qqHWWEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(qqHWWEstimation, self).__init__(
            name="qqHWW125",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            Weight("0.0857883", "crossSectionPerEventWeight"),

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "VBF.*HToWWTo2L2Nu_M125",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class VHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="VH",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2",)

    def get_cuts(self):
        return Cuts(Cut("(htxs_stage1p1cat>=300)&&(htxs_stage1p1cat<=505)", "htxs_match"))

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
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="WH",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("(htxs_stage1p1cat>=300)&&(htxs_stage1p1cat<=305)", "htxs_match"))

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
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ZH",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut("(htxs_stage1p1cat>=400)&&(htxs_stage1p1cat<=405)", "htxs_match"))

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


# Stage 0 and Stage 1.1 binning of ggH chosen by given name. If no name match, default is stage 0
class ggHEstimation(HTTEstimation):
    htxs_dict = ggH_htxs
    def __init__(self, name,  era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_weights(self):
        weights = super(ggHEstimation, self).get_weights()
        weights.remove("numberGeneratedEventsWeight")
        weights.add(Weight("8.22976e-8", "numberGeneratedEventsWeight"))
        weights.add(Weight("ggh_NNLO_weight", "gghNNLO"))
        weights.add(Weight("1.01", "bbh_inclusion_weight"))
        return weights

    def get_cuts(self):
        return Cuts(Cut(self.htxs_dict.get(self.name, self.htxs_dict["ggH125"]), "htxs_match"))

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


# Stage 0 and Stage 1.1 binning of qqH chosen by given name. If no name match, default is stage 0
class qqHEstimation(HTTEstimation):
    htxs_dict = qqH_htxs
    def __init__(self, name, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name=name,
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_cuts(self):
        return Cuts(Cut(self.htxs_dict.get(self.name, self.htxs_dict["qqH125"]), "htxs_match"))

    def get_files(self):
        query = {
            "process": "(^VBFHToTauTau.*125.*|^W(minus|plus)HToTauTau.*125.*|^ZHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class SUSYggHEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, mass, contribution, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(SUSYggHEstimation, self).__init__(
            name="_".join(["gg"+contribution,str(mass)]),
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")
        self.mass = mass
        self.contribution = contribution

    def get_weights(self):
        contribution_weight = "1.0"
        if self.contribution in ["A_i", "A_t", "A_b", "H_i", "H_t", "H_b", "h_i", "h_t", "h_b"]:
            contribution_weight = "gg%s_weight"%self.contribution

        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight(contribution_weight, "contributionWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            #Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"), # not needed, since should be at 1.0

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            # self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^SUSYGluGluToHToTauTau_M{MASS}$".format(MASS=self.mass),
            "data": False,
            "campaign": self._mc_campaign
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class SUSYbbHEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, mass, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(SUSYbbHEstimation, self).__init__(
            name="_".join(["bbH",str(mass)]),
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIFall17MiniAODv2")
        self.mass = mass

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            Weight("numberGeneratedEventsWeight",
                   "numberGeneratedEventsWeight"),
            #Weight("crossSectionPerEventWeight", "crossSectionPerEventWeight"), # not needed, since should be at 1.0

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            # self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("prefiringweight", "prefireWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^SUSYGluGluToBBHToTauTau_M{MASS}$".format(MASS=self.mass),
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "amcatnlo-pythia8",
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class bbH120Estimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="bbH120",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_files(self):
        query = {
           "process": "(^SUSYGluGluToBBHToTauTau.*120$)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class bbH130Estimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="bbH130",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_files(self):
        query = {
            "process": "(^SUSYGluGluToBBHToTauTau.*130$)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class ttHEstimation(HTTEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ttH",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIFall17MiniAODv2")

    def get_files(self):
        query = {
            "process": "(^ttHToTauTau.*125.*)",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "powheg\-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

class FakeEstimationLT(DataEstimation2016):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(DataEstimation2016, self).__init__(
            name="fakes",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
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


'''class AntiFakeEstimationLT(DataEstimation2016):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(DataEstimation2016, self).__init__(
            name="antifakes",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def get_weights(self):
        return Weights(Weight("(1.0-ff2_nom)", "anti_fake_factor"))

    def create_root_objects(self, systematic):
        aiso_systematic = copy.deepcopy(systematic)
        aiso_systematic.category.cuts.remove("tau_iso")
        aiso_systematic.category.cuts.add(
            Cut(
                "byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5&&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                "tau_aiso"))
        return super(AntiFakeEstimationLT,
                     self).create_root_objects(aiso_systematic)


class DataFakeEstimationLT(DataEstimation2016):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(DataEstimation2016, self).__init__(
            name="datafakes",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign=None)
        self._channel = channel

    def create_root_objects(self, systematic):
        aiso_systematic = copy.deepcopy(systematic)
        aiso_systematic.category.cuts.remove("tau_iso")
        aiso_systematic.category.cuts.add(
            Cut(
                "byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5&&byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                "tau_aiso"))
        return super(DataFakeEstimationLT,
                     self).create_root_objects(aiso_systematic)


class NewFakeEstimationLT(SumUpEstimationMethod):
    def __init__(self, era, directory, channel, friend_directory=None):
        super(NewFakeEstimationLT, self).__init__(
            name="fakes",
            folder="nominal",
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            factors=[1.0,-1.0],
            processes=[
                Process(
                    "dataFakes",
                    DataFakeEstimationLT(
                        era,
                        directory,
                        channel,
                        friend_directory=friend_directory)),
                Process(
                    "antiFakes",
                    AntiFakeEstimationLT(
                        era,
                        directory,
                        channel,
                        friend_directory=friend_directory))
            ])'''


class FakeEstimationTT(DataEstimation2016):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(DataEstimation2016, self).__init__(
            name="fakes",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
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


class NewFakeEstimationLT(NewFakeEstimationMethodLT):
    def __init__(self,
            era,
            directory,
            channel,
            nofake_processes,
            data_process,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(NewFakeEstimationLT, self).__init__(
            name="jetFakes",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
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
    def __init__(self,
            era,
            directory,
            channel,
            nofake_processes,
            data_process,
            friend_directory=None,
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(NewFakeEstimationTT, self).__init__(
            name="jetFakes",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
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
            fakeWeightstring="(0.5*ff1_nom*(byTightIsolationMVArun2017v2DBoldDMwLT2017_1<0.5)+0.5*ff2_nom*(byTightIsolationMVArun2017v2DBoldDMwLT2017_2<0.5))")
