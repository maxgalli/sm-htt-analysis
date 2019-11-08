import ROOT
from ROOT import TChain
import os

def CreateTChainFromPath(*args):
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
                chain.Add('/'.join([path_to_dir, subdir_or_file]))
            else:
                '''
                case 2:
                .
                ├── subdir_one
                │   └── one.root
                └── subdir_two
                    └── two.root
                '''
                file_name = '/'.join([path_to_dir, subdir_or_file, os.listdir(path_to_dir)[0]])
                if '.root' in file_name:
                    chain.Add(file_name)
                else:
                    pass # should raise Error??
    return chain
