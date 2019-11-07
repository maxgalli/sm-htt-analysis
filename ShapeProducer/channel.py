# -*- coding: utf-8 -*-

from cutstring import Cuts, Cut
import logging
logger = logging.getLogger(__name__)
"""
"""


class Channel(object):
    @property
    def cuts(self):
        return self._cuts

    @property
    def name(self):
        return self._name


class EESM(Channel):
    def __init__(self):
        self._name = "ee"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.1 && iso_2<0.1", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("(trg_singleelectron==1 && pt_1>26 && pt_2>26)",
                "trg_singleelectron"))


class MMSM(Channel):
    def __init__(self):
        self._name = "mm"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15 && iso_2<0.15", "muon_iso"), Cut(
                "q_1*q_2<0", "os"),
            Cut("(trg_singlemuon==1 && pt_1>25 && pt_2>25)", "trg_singlemuon"))


class MMSM2016(Channel):
    def __init__(self):
        self._name = "mm"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15 && iso_2<0.15", "muon_iso"), Cut(
                "q_1*q_2<0", "os"),
            Cut("m_vis > 50","m_vis_cut"),
            Cut("(pt_1 > 23 && trg_singlemuon==1)&&(0<1)", "trg_selection"))


class MMSM2017(Channel):
    def __init__(self):
        self._name = "mm"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15 && iso_2<0.15", "muon_iso"), Cut(
                "q_1*q_2<0", "os"),
            Cut("m_vis > 50","m_vis_cut"),
            Cut("(trg_singlemuon_27==1 || trg_singlemuon_24==1)", "trg_selection"))


class MMSM2018(Channel):
    def __init__(self):
        self._name = "mm"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("iso_1<0.15 && iso_2<0.15", "muon_iso"), Cut(
                "q_1*q_2<0", "os"),
            Cut("m_vis > 50","m_vis_cut"),
            Cut("(trg_singlemuon_27==1 || trg_singlemuon_24==1)", "trg_selection"))


class MT(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonTight3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2v1DBoldDMwLT_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("trg_singlemuon==1", "trg_singlemuon"))


class MTMSSM2017(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonTight3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && ((trg_singlemuon_27 == 1) || (trg_singlemuon_24 == 1) || (pt_1 < 25 && trg_crossmuon_mu20tau27 == 1))",
                "trg_selection"))


class MTSM2017(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonTight3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && ((trg_singlemuon_27 == 1) || (trg_singlemuon_24 == 1) || (pt_1 < 25 && trg_crossmuon_mu20tau27 == 1))",
                "trg_selection"))

class MTSM2018(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonTight3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && ((trg_singlemuon_27 == 1) || (trg_singlemuon_24 == 1)) || (pt_1 < 25 && (trg_crossmuon_mu20tau27_hps == 1 || trg_crossmuon_mu20tau27 == 1))",
                "trg_selection"))

class MTSM2016(Channel):
    def __init__(self):
        self._name = "mt"
        self._cuts = Cuts(
            Cut("flagMETFilter==1", "met_filter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonTight3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_2>0.5","againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "muon_iso"), 
            Cut("q_1*q_2<0", "os"),
            Cut("mt_1<50", "m_t"),
            Cut("pt_2>30 && ((trg_singlemuon == 1) || (trg_mutaucross == 1 && pt_1 <= 23 && pt_1 > 20))","trg_selection")
            )

class ET(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronTightMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_iso"),
            Cut("iso_1<0.1", "ele_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("trg_singleelectron==1", "trg_singleelectron"))


class ETMSSM2017(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronTightMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && pt_1 > 25 && (((trg_singleelectron_35 == 1) || (trg_singleelectron_32 == 1) || ((trg_singleelectron_27 == 1))) || (abs(eta_1)>1.5 && isEmbedded)) || (pt_1>25 && pt_1<28 && pt_2>35 && ((isEmbedded && (abs(eta_1)>1.5)) || (trg_crossele_ele24tau30 == 1)))",
                "trg_selection"))


class ETSM2017(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronTightMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && pt_1 > 25 && (((trg_singleelectron_35 == 1) || (trg_singleelectron_32 == 1) || ((trg_singleelectron_27 == 1))) || (abs(eta_1)>1.5 && isEmbedded)) || (pt_1>25 && pt_1<28 && pt_2>35 && ((isEmbedded && (abs(eta_1)>1.5)) || (trg_crossele_ele24tau30 == 1)))",
                "trg_selection"))

class ETSM2018(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronTightMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_iso"),
            Cut("iso_1<0.15", "ele_iso"), Cut("q_1*q_2<0", "os"),
            Cut("pt_2>30 && ((trg_singleelectron_35 == 1) || (trg_singleelectron_32 == 1) || (pt_1>25 && pt_1<33 && pt_2>35 && (trg_crossele_ele24tau30_hps == 1 || trg_crossele_ele24tau30 == 1)))",
                "trg_selection"))

class ETSM2016(Channel):
    def __init__(self):
        self._name = "et"
        self._cuts = Cuts(
            Cut("flagMETFilter==1", "met_filter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_2>0.5", "againstMuonDiscriminator"),
            Cut("againstElectronTightMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_iso"),
            Cut("iso_1<0.1", "ele_iso"), 
            Cut("q_1*q_2<0", "os"), 
            Cut("mt_1<50", "m_t"),
            Cut("(trg_singleelectron==1 && pt_1>26 && pt_2>30)","trg_singleelectron"))


class TT(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_1>0.5 && againstMuonLoose3_2>0.5",
                "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_1>0.5 && againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5", "tau_1_iso"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_2_iso"),
            Cut("q_1*q_2<0", "os"), Cut("trg_doubletau==1", "trg_doubletau"))


class TTSM2016(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("flagMETFilter==1", "met_filter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_1>0.5 && againstMuonLoose3_2>0.5",
                "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_1>0.5 && againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5", "tau_1_iso"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5", "tau_2_iso"),
            Cut("q_1*q_2<0", "os"), 
            Cut("pt_tt>50", "pt_h"),
            Cut("trg_doubletau==1", "trg_doubletau"))


class TTMSSM2017(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_1>0.5 && againstMuonLoose3_2>0.5",
                "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_1>0.5 && againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5",
                "tau_1_iso"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                "tau_2_iso"), Cut("q_1*q_2<0", "os"),
            Cut("(trg_doubletau_35_tightiso_tightid == 1) || (trg_doubletau_40_mediso_tightid == 1) || (trg_doubletau_40_tightiso == 1)",
                "trg_selection"))


class TTSM2017(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_1>0.5 && againstMuonLoose3_2>0.5",
                "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_1>0.5 && againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5",
                "tau_1_iso"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                "tau_2_iso"), Cut("q_1*q_2<0", "os"),
            Cut("(trg_doubletau_35_tightiso_tightid == 1) || (trg_doubletau_40_mediso_tightid == 1) || (trg_doubletau_40_tightiso == 1)",
                "trg_selection"))

class TTSM2018(Channel):
    def __init__(self):
        self._name = "tt"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("againstMuonLoose3_1>0.5 && againstMuonLoose3_2>0.5",
                "againstMuonDiscriminator"),
            Cut("againstElectronVLooseMVA6_1>0.5 && againstElectronVLooseMVA6_2>0.5",
                "againstElectronDiscriminator"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_1>0.5",
                "tau_1_iso"),
            Cut("byTightIsolationMVArun2017v2DBoldDMwLT2017_2>0.5",
                "tau_2_iso"), Cut("q_1*q_2<0", "os"),
            Cut("(((!isMC && run>=317509) || isMC) && (trg_doubletau_35_mediso_hps == 1)) || (!isMC && (run<317509) && ((trg_doubletau_35_tightiso_tightid == 1) || (trg_doubletau_40_mediso_tightid == 1) || (trg_doubletau_40_tightiso == 1)))",
                "trg_selection"))

class EM(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("trg_muonelectron==1", "trg_muonelectron"))

class EMSM2016(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("nbtag==0 && pZetaMissVis>-35 && mTdileptonMET<60","dzeta"),
            Cut("pt_2>10 && ((pt_1>13 && pt_2>24 && trg_muonelectron_mu23ele12 == 1) || (pt_1>24 && pt_2>10 && trg_muonelectron_mu8ele23 == 1))",
                "trg_selection"))

class EMMSSM2017(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("pt_2>10 && ((trg_muonelectron_mu23ele12 == 1) || (trg_muonelectron_mu8ele23 == 1))",
                "trg_selection"))

class EMSM2017(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("pt_2>10 && ((trg_muonelectron_mu23ele12 == 1) || (trg_muonelectron_mu8ele23 == 1))",
                "trg_selection"))

class EMSM2018(Channel):
    def __init__(self):
        self._name = "em"
        self._cuts = Cuts(
            Cut("flagMETFilter == 1", "METFilter"),
            Cut("extraelec_veto<0.5", "extraelec_veto"),
            Cut("extramuon_veto<0.5", "extramuon_veto"),
            Cut("dilepton_veto<0.5", "dilepton_veto"),
            Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
            Cut("q_1*q_2<0", "os"),
            Cut("(trg_muonelectron_mu23ele12 == 1 && pt_1>13 && pt_2 > 24) || (trg_muonelectron_mu8ele23 == 1 && pt_1>24 && pt_2>10)",
                "trg_selection"))

class PU(Channel):
    def __init__(self):
        self._name = "pu"
        self._cuts = Cuts()


# collection of channels an analysis can be ran on
class Channels(object):
    def __init__(self, name):
        self._name = name
        self._channels = []

    def add(self, channel):
        self._channels.append(channel)

    @property
    def name(self):
        return self._name
