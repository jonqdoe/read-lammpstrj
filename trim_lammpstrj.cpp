#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <algorithm>
using namespace std;

#include "lammpstrj.h"
#define Dim 3

FILE *fopen_mkdir(char *path, char *mode);

int trim_lammpstrj(const char* name, float border, bool sep_files) {

    std::string str;

    std::vector<int> accepted_molecules(nsites);
    std::vector<int> rejected_molecules(nsites);

    string filename(name);
    string filename_extended;

    bool rejected=false;

    int accepted_ns = 0;
    int mol_count = 0;

    for (int t = 0; t < nframes; t++) {

        accepted_molecules.resize(0);
        rejected_molecules.resize(0);
        accepted_ns = 0;
        mol_count = 0; 

        for (int ns = 0; ns < nsites; ns++) {
            if (ns==0) {mol_count = 1;}
            else if ( (mol[ns] == mol[ns-1] && rejected)){
            continue;
            }
            else if (mol[ns] == mol[ns-1] && not rejected){
                mol_count++;
            }
            else if (mol[ns] != mol[ns-1] && not rejected){
                accepted_molecules.push_back(mol[ns-1]);
                accepted_ns += mol_count; 
                mol_count = 1; 
                rejected=false;
            }
            else if( mol[ns] != mol[ns-1] && rejected ){
                mol_count = 1;
                rejected=false;
            }


            for (int j = 0; j < Dim; j++) {
                if (abs(xt[t][ns][j] - xlo[t][j]) < abs(border) 
                    || abs(xt[t][ns][j] - xhi[t][j]) < abs(border)  ){
                        rejected=true;
                        rejected_molecules.push_back(mol[ns]);
                        break;
                    }
            }
            // if (rejected) rejected_molecules.push_back(mol[ns]);
        }

        FILE *opt;
        if (sep_files == true){
            filename_extended = "data/" + filename + "_" + to_string(t);
            opt = fopen_mkdir((char*) filename_extended.c_str(),"w");
        } else {
            filename_extended = "data/" + filename;
            if (t==0)
                opt = fopen_mkdir((char*) filename_extended.c_str(),"w");
        }

        fprintf(opt,"ITEM: TIMESTEP\n%d\n", t);
        fprintf(opt,"ITEM: NUMBER OF ATOMS\n%d\n", accepted_ns);
        fprintf(opt,"ITEM: BOX BOUNDS pp pp pp\n");
        
        for ( int j=0 ; j<3 ; j++ ) {
        fprintf(opt,"%f %f\n", xlo[t][j],  xhi[t][j]);
        }

        fprintf(opt,"ITEM: ATOMS id type mol x y z\n");
        

        int new_id = 1;
        int new_mol_id = 0;
        int last_mol_id = 0;
        for (int ns = 0; ns < nsites; ns++){
            if (std::find(accepted_molecules.begin(), accepted_molecules.end(),mol[ns])
            != accepted_molecules.end()){
                // If the first particle for the frame
                if (new_id==0){
                    last_mol_id = mol[ns];
                    new_mol_id++;
                }
                // If not the first particle
                else{
                    // Is the molecule id different
                    if (last_mol_id != mol[ns]){
                        new_mol_id++;
                        last_mol_id = mol[ns];
                    }
                }
                fprintf(opt,"%d %d %d",new_id++,type[ns],new_mol_id);
                for (int j = 0; j < Dim; j++) {
                    fprintf(opt," %f", xt[t][ns][j]);
                }
                fprintf(opt,"\n");
            }
            

        }

        if (sep_files == true){
            fclose(opt);
        } else {
            if (t==nframes) fclose(opt);
        }
    }
}