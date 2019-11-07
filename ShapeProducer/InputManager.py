import ROOT
from ROOT import TChain

def CreateTChain(chain_name = "HTTInputChain", *args):
    chain = TChain(chain_name)
    for input_file in args:
        chain.Add(input_file)
    return chain
