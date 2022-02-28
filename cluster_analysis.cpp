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
double pbc_mdr2(std::vector<double> ri,std::vector<double> rj,std::vector<double> L );
double pbc_mdr2(std::vector<double> ri,std::vector<double> rj,std::vector<double> dr,std::vector<double> L );

int cluster_analysis(std::vector<int> tp, float cutoff) {

    double cutoff_2 = cutoff*cutoff;
    FILE *opt = fopen_mkdir("data/clusters.txt","w");

    fprintf(opt,"time cluster_id particle_id molecule_id\n");
    
    vector<int> atom_og(nsites);

    int particles_of_type_tp=0;
    for (int i = 0; i<nsites; i++){
        if (std::find(tp.begin(), tp.end(), type[i]) != tp.end()) {
            atom_og.at(particles_of_type_tp) = i;
            particles_of_type_tp++;
        }
    }
    atom_og.resize(particles_of_type_tp);

    vector<int> atom_id;
    vector<double> dr(3);

    for (int t = 0; t < nframes; t++) {
        
        std::vector<vector<int>> clusters{};
        clusters.reserve(nsites);

        atom_id = atom_og;

        while (atom_id.size() > 0){
            vector<int> cluster_temp { atom_id.at(0) }; 
            atom_id.erase(atom_id.begin());

            for (int i = 0; i < cluster_temp.size(); i++){
                for (int j = 0; j < atom_id.size(); j++){
                    if (
                        pbc_mdr2(xt[t][cluster_temp.at(i)], xt[t][atom_id.at(j)], dr, L[t]) 
                        <= cutoff_2){
                        cluster_temp.push_back(atom_id.at(j));
                        atom_id.erase(atom_id.begin()+j);
                        j--;
                    }
                }
            }
            clusters.push_back(cluster_temp);
        }

        printf("\rTimestep %d out of %d", t, nframes); fflush( stdout );

        for (int i = 0; i < clusters.size(); i++){
            for (const int& particle: clusters.at(i)){
                fprintf(opt,"%d %d %d %d\n",t, i, particle, mol[particle]);
            }
        }

    }

    fclose(opt);

}
