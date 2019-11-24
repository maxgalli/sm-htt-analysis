#include <sys/types.h>
#include <dirent.h>

// ShapeProducer headers
#include "Utils.h"

// ROOT headers
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"

// Read inside a directory and put all the files contained
// in a vector of strings, after appending them to the full
// path.
// Args:
// - name: reference to a string containing the name of the
// directory;
// - v: string vector filled with the full path to the files
// inside the directory called "name".
void read_directory(const std::string&, stringvec&);

// Read inside a ROOT file and put all the TDFs contained in
// a vector of TDFs.
// Args:
// - root_file: pointer to a TFile object
// - v: TDirectoryFile vector filled with the TDFs inside root_file
void read_root_file(TFile*, tdfvec&);

// Variadic function: takes a fixed argument (const char* tree_name)
// which is the name of the tree contained in the ROOT files and a
// variable number of templated arguments, which are strings containing
// the path to the directories with inside other directories containing
// only one ROOT file each (see sketch below).
template<typename...Args>
TChain* CreateTChainFromPath(const char*, Args...);
