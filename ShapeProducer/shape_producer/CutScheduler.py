import ROOT

class CutScheduler:
    def __init__(self, dataframe, cuts = None, weights = None):
        self._in_rdf = dataframe
        self._cuts = cuts

    @property
    def in_rdf(self):
        return self._in_rdf

    @in_rdf.setter
    def in_rdf(self, a_dataframe):
        self._in_rdf = a_dataframe

    @property
    def cuts(self):
        return self._cuts

    @cuts.setter
    def cuts(self, a_cuts):
        self._cuts = a_cuts

    def apply_cuts(self):
        internal_rdf = self._in_rdf
        for cut in self._cuts:
            self._out_rdf = internal_rdf.Filter(cut.string)
            internal_rdf = self.out_rdf
        return self._out_rdf

    @property
    def out_rdf(self):
        return self._out_rdf
