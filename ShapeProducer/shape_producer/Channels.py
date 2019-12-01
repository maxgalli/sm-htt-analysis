import ROOT
from CutScheduler import CutScheduler
from Cutstring import Cut
from Cutstring import CutVec


class Channel(CutScheduler):
    def __init__(self,
            name,
            dataframe,
            cuts,
            era = None,
            model = None,
            weights = None):
        self._name = name
        self._era = era
        self._model = model
        CutScheduler.__init__(self, dataframe, cuts, weights)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, a_name):
        self._name = a_name

    @property
    def era(self):
        return self._era

    @era.setter
    def era(self, a_era):
        self._era = a_era

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, a_model):
        self._model = a_model

def create_channels(name, era, model, dataframe):
    if name == "em":
        if era == "2016":
            cut_vec = CutVec(Cut("flagMETFilter == 1", "METFilter"),
                    Cut("extraelec_veto<0.5", "extraelec_veto"),
                    Cut("extramuon_veto<0.5", "extramuon_veto"),
                    Cut("dilepton_veto<0.5", "dilepton_veto"),
                    Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
                    Cut("q_1*q_2<0", "os"),
                    Cut("nbtag==0 && pZetaMissVis>-35 && mTdileptonMET<60","dzeta"),
                    Cut("pt_2>10 && ((pt_1>13 && pt_2>24 && trg_muonelectron_mu23ele12 == 1) || (pt_1>24 && pt_2>10 && trg_muonelectron_mu8ele23 == 1))","trg_selection"))
        elif era == "2017":
            if model == "MSSM":
                cut_vec = CutVec(Cut("flagMETFilter == 1", "METFilter"),
                        Cut("extraelec_veto<0.5", "extraelec_veto"),
                        Cut("extramuon_veto<0.5", "extramuon_veto"),
                        Cut("dilepton_veto<0.5", "dilepton_veto"),
                        Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
                        Cut("q_1*q_2<0", "os"),
                        Cut("pt_2>10 && ((trg_muonelectron_mu23ele12 == 1) || (trg_muonelectron_mu8ele23 == 1))",
                            "trg_selection"))
            else:
                cut_vec = CutVec(Cut("flagMETFilter == 1", "METFilter"),
                        Cut("extraelec_veto<0.5", "extraelec_veto"),
                        Cut("extramuon_veto<0.5", "extramuon_veto"),
                        Cut("dilepton_veto<0.5", "dilepton_veto"),
                        Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
                        Cut("q_1*q_2<0", "os"),
                        Cut("pt_2>10 && ((trg_muonelectron_mu23ele12 == 1) || (trg_muonelectron_mu8ele23 == 1))",
                            "trg_selection"))
        else:
            cut_vec = CutVec(Cut("flagMETFilter == 1", "METFilter"),
                    Cut("extraelec_veto<0.5", "extraelec_veto"),
                    Cut("extramuon_veto<0.5", "extramuon_veto"),
                    Cut("dilepton_veto<0.5", "dilepton_veto"),
                    Cut("iso_1<0.15", "ele_iso"), Cut("iso_2<0.2", "muon_iso"),
                    Cut("q_1*q_2<0", "os"),
                    Cut("(trg_muonelectron_mu23ele12 == 1 && pt_1>13 && pt_2 > 24) || (trg_muonelectron_mu8ele23 == 1 && pt_1>24 && pt_2>10)", "trg_selection"))

    channel = Channel(name, dataframe, cut_vec, era, model)
    return channel

