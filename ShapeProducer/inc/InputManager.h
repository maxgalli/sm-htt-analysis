#include <sys/types.h>
#include <dirent.h>

#include "Utils.h"

// ROOT headers
#include "TTree.h"
#include "TChain.h"

// Read inside a directory and put all the files contained
// in a vector of strings, after appending them to the full
// path.
// Args:
// - name: reference to a string containing the name of the
// directory;
// - v: string vector filled with the full path to the files
// inside the directory called "name".
void read_directory(const std::string& name, stringvec& v) {
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        const char* last_char = &name.back();
        const char* file_name = dp->d_name;
        std::string full_path;
        // Exclude "." and ".."
        if (strcmp(file_name, ".") == 0 || strcmp(file_name, "..") == 0)
            continue;
        // Don't add "/" if it is already the last character in the
        // directory name
        if (strcmp(last_char, "/") == 0) {
            full_path = name + file_name;
        } else {
            // Add "/" if it is not the last character in the
            // directory name
            full_path = name + "/" + file_name;
        }
        v.push_back(full_path);
    }
    closedir(dirp);
}

// Read inside a ROOT file and put all the TDFs contained in
// a vector of TDFs.
// Args:
// - root_file: pointer to a TFile object
// - v: TDirectoryFile vector filled with the TDFs inside root_file
void read_root_file(TFile* root_file, tdfvec& v) {
    TList* keys = root_file->GetListOfKeys();
    for (auto key : *keys) {
        const char* tdf_name = key->GetName();
        TDirectoryFile* tdf = (TDirectoryFile*)root_file->GetDirectory(tdf_name);
        v.push_back(tdf);
    }
}

namespace shape_producer {
    namespace InputManager {
        // Variadic function: takes a fixed argument (const char* tree_name)
        // which is the name of the tree contained in the ROOT files and a
        // variable number of templated arguments, which are strings containing
        // the path to the directories with inside other directories containing
        // only one ROOT file each (see sketch below).
        template<typename...Args>
            TChain* CreateTChainFromPath(const char* tree_name, Args... path_names){
                TChain* f_chain = new TChain(tree_name);
                TFile* root_file;
                std::vector<Args...> path_names_list = {path_names...};
                // Loop on the paths passed as arguments in the shell scripts
                for (auto path : path_names_list) {
                    /*
                    .path
                    ├── subdir_one
                    │   └── one.root
                    └── subdir_two
                        └── two.root
                    */
                    stringvec rf_containers_vec;
                    read_directory(path, rf_containers_vec);
                    // Loop on the subdirs containing one ROOT file each
                    for (auto rf_container : rf_containers_vec) {
                        stringvec rfs_vec;
                        tdfvec tdfs_vec;
                        read_directory(rf_container, rfs_vec);
                        // Only one ROOT file in the subdir
                        const char* file_path = &rfs_vec[0][0];
                        root_file = new TFile(file_path);
                        read_root_file(root_file, tdfs_vec);
                        char full_tree_path[512];
                        // Loop on the TDFs
                        for (auto tdf : tdfs_vec) {
                            strncpy(full_tree_path, file_path, sizeof(full_tree_path));
                            strncat(full_tree_path, "/", 1);
                            strncat(full_tree_path, tdf->GetName(), sizeof(full_tree_path) - strlen(full_tree_path) - 1);
                            strncat(full_tree_path, "/", 1);
                            strncat(full_tree_path, tree_name, sizeof(full_tree_path) - strlen(full_tree_path) - 1);
                            f_chain->Add(full_tree_path);
                        }
                    }
                }
                return f_chain;
            }
    }
}

