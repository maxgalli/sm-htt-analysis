import ROOT
from ROOT import TChain
import os

def GetListOfSubTDirectoryFiles(root_file):
    sub_tdirf = []
    for key in root_file.GetListOfKeys():
        sub_dir_name = key.GetName()
        sub_tdirf.append(root_file.GetDirectory(sub_dir_name))
    return sub_tdirf

def CreateTChainFromPath(tree_name, *args):
    # Parameters:
    # every parameter is a path to a directory
    chain = TChain()
    for path_to_dir in args:
        dir_content = os.listdir(path_to_dir)
        for subdir_or_file in dir_content:
            if '.root' in subdir_or_file:
                '''
                case 1:
                .
                ├── one.root
                └── two.root
                '''
                file_name = '/'.join([path_to_dir, subdir_or_file])
            else:
                '''
                case 2:
                .
                ├── subdir_one
                │   └── one.root
                └── subdir_two
                    └── two.root
                '''
                path_to_subdir = '/'.join([path_to_dir, subdir_or_file])
                file_name = '/'.join([path_to_dir, subdir_or_file, os.listdir(path_to_subdir)[0]])
                if '.root' not in file_name:
                    raise AttributeError(file_name)

            file_object = ROOT.TFile(file_name)
            file_content = GetListOfSubTDirectoryFiles(file_object)
            for sub_tdirf in file_content:
                chain.Add('/'.join([file_name, sub_tdirf.GetName(), tree_name]))

    return chain
