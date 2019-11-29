import ROOT

class Cut():
    def __init__(self, string, name):
        self._name = name
        self._string = string

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, a_name):
        self._name = a_name

    @property
    def string(self):
        return self._string

    @string.setter
    def string(self, a_string):
        self._string = a_string
