import ROOT

class Cut:
    def __init__(self, expression, name):
        self._expression = expression
        self._name = name

    @property
    def expression(self):
        return self._expression

    @expression.setter
    def expression(self, a_expression):
        self._expression = a_expression

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, a_name):
        self._name = a_name

class CutVec(list):
    def __init__(self, *args):
        for cut in args:
            self.append(cut)
