

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
    weight = Weight("1.0","triggerweight")

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
        weight = Weight("(crossTriggerMCWeight_1*(crossTriggerMCWeight_1<10 && crossTriggerMCWeight_1>0.1)+(crossTriggerMCWeight_1>10 || crossTriggerMCWeight_1<0.1))*(pt_1<25) + (trigger_24_27_Weight_1*(pt_1>25))","triggerweight")

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
        weight = Weight("(crossTriggerMCWeight_1*(crossTriggerMCWeight_1<10)+(crossTriggerMCWeight_1>10))*(pt_1<33)+((pt_1>=33)*trigger_32_35_Weight_1)","triggerweight")

    elif "tt" in channel:
        DiTauMC = "*".join([MCTau_1,MCTau_2])
        DiTauData = DiTauMC.replace("MC","Data")
        DiTau = "("+DiTauData+")/("+DiTauMC+")"
        weight = Weight(DiTau,"triggerweight")

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
        weight = Weight("((gen_match_2 == 5)*0.90 + (gen_match_2 != 5))", "taubyIsoIdWeight")
    elif "tt" in channel:
        weight = Weight("((gen_match_1 == 5)*0.90 + (gen_match_1 != 5))*((gen_match_2 == 5)*0.90 + (gen_match_2 != 5))", "taubyIsoIdWeight")
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
            folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            qcd_weight=Weight("1.0","qcd_Weight")):
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
            qcd_weight=qcd_weight)

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            # self.get_eleHLTZvtxWeight_for_channel(self.channel.name),
            Weight("118.7*(abs(crossSectionPerEventWeight - 63.21) < 0.01) + crossSectionPerEventWeight*(abs(crossSectionPerEventWeight - 63.21) > 0.01)", "crossSectionPerEventWeight"),

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
            "generator": "powheg\-pythia8",
            "campaign": self._mc_campaign
         }
        files += self.era.datasets_helper.get_nicks_with_query(query)

        # query = {
        #     "process": "ST",  # Query for Single-Top samples (newer pileup mixing)
        #     "data": False,
        #     "scenario": "^PU2017newpmx$",
        #     "generator": "powheg\-pythia8",
        #     "campaign": self._mc_campaign
        # }
        # files += self.era.datasets_helper.get_nicks_with_query(query)
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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            # self.get_eleHLTZvtxWeight_for_channel(self.channel.name),

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
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
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
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_weights(self):
        if self.atNLO:
            z_stitching_weight = Weight("((genbosonmass >= 50.0) * 3.3847e-05 + (genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight)","z_stitching_weight") # xsec_NNLO [pb] = 2075.14*3, N_inclusive_NLO = 183925831, xsec_NNLO/N_inclusive_NLO = 3.3847e-05; fraction of negative events in 'generatorWeight'
        else:
            z_stitching_weight = Weight("((genbosonmass >= 50.0)*0.00005754202*((npartons == 0 || npartons >= 5)*1.0 + (npartons == 1)*0.194267667208 + (npartons == 2)*0.21727746547 + (npartons == 3)*0.26760465744 + (npartons == 4)*0.294078683662) + (genbosonmass < 50.0)*numberGeneratedEventsWeight*crossSectionPerEventWeight)", "z_stitching_weight")
              # xsec_NNLO [pb] = 2075.14*3, N_inclusive = 100194597,  xsec_NNLO/N_inclusive = 0.00005754202 [pb] weights: [1.0, 0.194267667208, 0.21727746547, 0.26760465744, 0.294078683662]
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
            Weight("zPtReweightWeight", "zPtReweightWeight"),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        queryM10 = {
            "process": "DYJetsToLL_M10to50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        queryM50_inclusive = {
            "process": "DYJetsToLLM50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        queryM50_1jet = {
            "process": "DY1JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        queryM50_2jet = {
            "process": "DY2JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
                    }
        queryM50_3jet = {
            "process": "DY3JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        queryM50_4jet = {
            "process": "DY4JetsToLL_M50",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph\-pythia8"
        }
        queryEWKZ = {
            "process": "EWKZ2Jets_ZToLL_M50",
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
            files = self.era.datasets_helper.get_nicks_with_query(queryM50_inclusive)
            files += self.era.datasets_helper.get_nicks_with_query(queryM50_1jet)
            files += self.era.datasets_helper.get_nicks_with_query(queryM50_2jet)
            files += self.era.datasets_helper.get_nicks_with_query(queryM50_3jet)
            files += self.era.datasets_helper.get_nicks_with_query(queryM50_4jet)
            files += self.era.datasets_helper.get_nicks_with_query(queryM10)
            files += self.era.datasets_helper.get_nicks_with_query(queryEWKZ)

        log_query(self.name, queryM50_inclusive, files)
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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
                Weight("isoWeight_1*idWeight_1*((pt_1>=25)*(trigger_24_27_Weight_1*(trigger_24_27_Weight_1<2.0)+(trigger_24_27_Weight_1>2.0))+(pt_1<25)*(crossTriggerEmbeddedWeight_2*crossTriggerEmbeddedWeight_1))", "lepton_sf"),
                Weight("(gen_match_2==5)*0.97+(gen_match_2!=5)", "emb_tau_id"),
                Weight("embeddedDecayModeWeight", "decayMode_SF"),
                Weight("gen_match_1==4 && gen_match_2==5","emb_veto"))
        elif self.channel.name in ["et"]:
            return Weights(
                Weight("generatorWeight",
                       "simulation_sf"),
                Weight("muonEffTrgWeight*muonEffIDWeight_1*muonEffIDWeight_2", "scale_factor"),
                Weight("idWeight_1*((crossTriggerEmbeddedWeight_2*(crossTriggerEmbeddedWeight_2<1)+(crossTriggerEmbeddedWeight_2>1))*(crossTriggerEmbeddedWeight_1*(crossTriggerEmbeddedWeight_1<10)+(crossTriggerEmbeddedWeight_1>10))*(pt_1<33)+(pt_1>=33)*trigger_32_35_Weight_1)*isoWeight_1", "lepton_sf"),
                Weight("(gen_match_2==5)*0.97+(gen_match_2!=5)", "emb_tau_id"),
                Weight("gen_match_1==3 && gen_match_2==5","emb_veto"),
                Weight("embeddedDecayModeWeight", "decayMode_SF"))
        elif self.channel.name == "tt":
            return Weights(
                Weight("generatorWeight",
                       "simulation_sf"),
                Weight("muonEffTrgWeight*muonEffIDWeight_1*muonEffIDWeight_2", "scale_factor"),
                Weight("(triggerWeight_1+triggerWeight_2)/2.0","tau1_leg_weight"),
                #~ Weight("(0.104*(pt_2>=30 && pt_2<35) + 0.519*(pt_2>=35 && pt_2<40) + 0.682*(pt_2>=40 && pt_2<45) + 0.766*(pt_2>=45 && pt_2<50) + 0.786*(pt_2>=50 && pt_2<60) + 0.804*(pt_2>=60 && pt_2<80) + 0.735*(pt_2>=80 && pt_2<100) + 0.730*(pt_2>=100 && pt_2<150) + 0.683*(pt_2>=150 && pt_2<200) + (pt_2>=200))","tau2_leg_weight"),
                Weight("((gen_match_1==5)*0.97+(gen_match_1!=5))*((gen_match_2==5)*0.97+(gen_match_2!=5))", "emb_tau_id"),
                Weight("gen_match_1==5 && gen_match_2==5","emb_veto"),
                Weight("embeddedDecayModeWeight", "decayMode_SF"))
        elif self.channel.name == "em":
            return Weights(
                Weight("generatorWeight", "simulation_sf"),
                Weight("muonEffTrgWeight*muonEffIDWeight_1*muonEffIDWeight_2", "scale_factor"),
                Weight("idWeight_1*isoWeight_1*idWeight_2*isoWeight_2",
                       "idiso_lepton_sf"),
                Weight("gen_match_1==3 && gen_match_2==4","emb_veto"))#,
                # Weight("(trigger_23_data_Weight_2/trigger_23_embed_Weight_2)*(pt_2>24)+(trigger_8_data_Weight_2/trigger_8_embed_Weight_2)*(pt_2<24)+(trigger_12_data_Weight_1/trigger_12_embed_Weight_1)*(pt_1<24)+(trigger_23_data_Weight_1/trigger_23_embed_Weight_1)*(pt_1<24)",
                    #    "trigger_lepton_sf"))


    def get_files(self):
        query = {"process": "Embedding2018", "embedded": True}
        if self.channel.name == "mt":
            query["campaign"] = "MuTauFinalState"
        elif self.channel.name == "et":
            query["campaign"] = "ElTauFinalState"
            query["process"] = "Embedding2018"
        elif self.channel.name == "tt":
            query["campaign"] = "TauTauFinalState"
        elif self.channel.name == "em":
            query["process"] = "Embedding2018"
            query["campaign"] = "ElMuFinalState"
        files = self.era.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return self.artus_file_names(files)

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
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_weights(self):
        return Weights(
            # MC related weights
            Weight("generatorWeight", "generatorWeight"),
            # Weight("numberGeneratedEventsWeight","numberGeneratedEventsWeight"), # to be used only for one inclusive sample
            #Weight("crossSectionPerEventWeight","crossSectionPerEventWeight"), # to be used only for one inclusive sample
            Weight("((0.00092600048*((npartons <= 0 || npartons >= 5)*1.0 + (npartons == 1)*0.1647043928 + (npartons == 2)*0.128547226623 + (npartons == 3)*0.0767138313139 + (npartons == 4)*0.0631529545476)) * (genbosonmass>=0.0) + numberGeneratedEventsWeight * crossSectionPerEventWeight * (genbosonmass<0.0))",
                 "wj_stitching_weight"), # xsec_NNLO [pb] = 61526.7, N_inclusive = 66443486, xsec_NNLO/N_inclusive = 0.00092600048 [pb] weights: [1.0, 0.1647043928, 0.128547226623, 0.0767138313139, 0.0631529545476]

            # Weights for corrections
            Weight("puweight", "puweight"),
            Weight("idWeight_1*idWeight_2","idweight"),
            Weight("isoWeight_1*isoWeight_2","isoweight"),
            Weight("trackWeight_1*trackWeight_2","trackweight"),
            self.get_triggerweight_for_channel(self.channel._name),
            # self.get_singlelepton_triggerweight_for_channel(self.channel.name),
            Weight("eleTauFakeRateWeight*muTauFakeRateWeight", "leptonTauFakeRateWeight"),
            self.get_tauByIsoIdWeight_for_channel(self.channel.name),
            # self.get_eleHLTZvtxWeight_for_channel(self.channel.name),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "WJetsToLNu",
            #"process": "WJetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files = self.era.datasets_helper.get_nicks_with_query(query)
        query = {
            "process": "W1JetsToLNu",
            #"process": "WJetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)
        query = {
            "process": "W2JetsToLNu",
            #"process": "WJetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)
        query = {
            "process": "W3JetsToLNu",
            #"process": "WJetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)
        query = {
            "process": "W4JetsToLNu",
            #"process": "WJetsToLNu",
            "data": False,
            "campaign": self._mc_campaign,
            "generator": "madgraph-pythia8"
        }
        files += self.era.datasets_helper.get_nicks_with_query(query)
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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            # self.get_eleHLTZvtxWeight_for_channel(self.channel.name),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "TTTo.*",
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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            # self.get_eleHLTZvtxWeight_for_channel(self.channel.name),

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
            mc_campaign="RunIIAutumn18MiniAOD")

    # def get_cuts(self):
    #     return Cuts(Cut("(htxs_stage1p1cat>=300)&&(htxs_stage1p1cat<=404)", "htxs_match"))

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
            mc_campaign="RunIIAutumn18MiniAOD")

    # def get_cuts(self):
    #     return Cuts(Cut("(htxs_stage1p1cat>=300)&&(htxs_stage1p1cat<=304)", "htxs_match"))

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
            mc_campaign="RunIIAutumn18MiniAOD")

    # def get_cuts(self):
    #     return Cuts(Cut("(htxs_stage1p1cat>=400)&&(htxs_stage1p1cat<=404)", "htxs_match"))

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


class ggHEstimation(HTTEstimation):
    htxs_dict = ggH_htxs
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
            mc_campaign="RunIIAutumn18MiniAOD")

    # def get_weights(self):
    #     weights = super(ggHEstimation, self).get_weights()
    #     # weights.remove("numberGeneratedEventsWeight")
    #     # weights.add(Weight("8.22976e-8", "numberGeneratedEventsWeight"))
    #     weights.add(Weight("ggh_NNLO_weight", "gghNNLO"))
    #     weights.add(Weight("1.01", "bbh_inclusion_weight"))
    #     return weights

    # def get_cuts(self):
    #     return Cuts(Cut("(htxs_stage1p1cat>=101)&&(htxs_stage1p1cat<=111)", "htxs_match"))

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
            mc_campaign="RunIIAutumn18MiniAOD")

    # def get_cuts(self):
    #     return Cuts(Cut("(htxs_stage1p1cat>=201)&&(htxs_stage1p1cat<=205)", "htxs_match"))

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


class ggHEstimation_VBFTOPO_JET3VETO(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_VBFTOPO_JET3VETO",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    # def get_cuts(self):
    #     return Cuts(Cut("htxs_stage1p1cat==101", "htxs_match"))


class ggHEstimation_VBFTOPO_JET3(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_VBFTOPO_JET3",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==102", "htxs_match"))


class ggHEstimation_0J(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_0J",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    # def get_cuts(self):
    #     return Cuts(Cut("htxs_stage1p1cat==103", "htxs_match"))


class ggHEstimation_1J_PTH_0_60(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_1J_PTH_0_60",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==104", "htxs_match"))


class ggHEstimation_1J_PTH_60_120(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_1J_PTH_60_120",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==105", "htxs_match"))


class ggHEstimation_1J_PTH_120_200(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_1J_PTH_120_200",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==106", "htxs_match"))


class ggHEstimation_1J_PTH_GT200(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_1J_PTH_GT200",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==107", "htxs_match"))


class ggHEstimation_GE2J_PTH_0_60(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_GE2J_PTH_0_60",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==108", "htxs_match"))


class ggHEstimation_GE2J_PTH_60_120(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_GE2J_PTH_60_120",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==109", "htxs_match"))


class ggHEstimation_GE2J_PTH_120_200(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_GE2J_PTH_120_200",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==110", "htxs_match"))


class ggHEstimation_GE2J_PTH_GT200(ggHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="ggH_GE2J_PTH_GT200",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==111", "htxs_match"))


class qqHEstimation_VBFTOPO_JET3VETO(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="qqH_VBFTOPO_JET3VETO",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==201", "htxs_match"))


class qqHEstimation_VBFTOPO_JET3(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="qqH_VBFTOPO_JET3",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==202", "htxs_match"))


class qqHEstimation_VH2JET(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="qqH_VH2JET",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==203", "htxs_match"))


class qqHEstimation_REST(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="qqH_REST",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==204", "htxs_match"))


class qqHEstimation_PTJET1_GT200(qqHEstimation):
    def __init__(self, era, directory, channel, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(HTTEstimation, self).__init__(
            name="qqH_PTJET1_GT200",
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            friend_directory=friend_directory,
            channel=channel,
            mc_campaign="RunIIAutumn18MiniAOD")

    def get_cuts(self):
        return Cuts(Cut("htxs_stage1p1cat==205", "htxs_match"))


class SUSYggHEstimation(EstimationMethod):
    def __init__(self, era, directory, channel, mass, friend_directory=None, folder="nominal",
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,):
        super(SUSYggHEstimation, self).__init__(
            name="_".join(["ggH",str(mass)]),
            folder=folder,
            get_triggerweight_for_channel=get_triggerweight_for_channel,
            get_singlelepton_triggerweight_for_channel=get_singlelepton_triggerweight_for_channel,
            get_tauByIsoIdWeight_for_channel=get_tauByIsoIdWeight_for_channel,
            get_eleHLTZvtxWeight_for_channel=get_eleHLTZvtxWeight_for_channel,
            era=era,
            directory=directory,
            channel=channel,
            friend_directory=friend_directory,
            mc_campaign="RunIIAutumn18MiniAOD")
        self.mass = mass

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
            # self.get_eleHLTZvtxWeight_for_channel(self.channel.name),

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
            mc_campaign="RunIIAutumn18MiniAOD")
        self.mass = mass

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
            # self.get_eleHLTZvtxWeight_for_channel(self.channel.name),

            # Data related scale-factors
            self.era.lumi_weight)

    def get_files(self):
        query = {
            "process": "^SUSYGluGluToBBHToTauTau_M{MASS}$".format(MASS=self.mass),
            "data": False,
            "campaign": self._mc_campaign
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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
            mc_campaign="RunIIAutumn18MiniAOD")

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
