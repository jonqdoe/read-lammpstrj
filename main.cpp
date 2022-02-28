#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <complex>
#include "stdlib.h"
#include <algorithm>
#include "fftw.h"
#include "io_utils.h"
#include "mesh_globals.h"
using namespace std;

#include "lammpstrj.h"

#define PI 3.141592654
#define PI2 6.2831853071

int read_lammpstrj(const char*, const int, const int) ;
void connect_molecules(vector<vector<vector<double>>>&, vector<int>, int, int, 
    vector<vector<double>>) ;
void remove_pbc_time_jumps(vector<vector<vector<double>>>&, int, int, vector<vector<double>>) ;
void lc_order(vector<vector<vector<double>>>, int, int, vector<int>, int, int);
void calc_msd(vector<vector<vector<double>>>, const vector<int>, const int, const int, vector<vector<double>> );
void calc_rdf(double, string, string);
void nlist_init(void);
void van_hove(vector<vector<vector<double>>>, int, int, int, double, vector<vector<double>> );
int trim_lammpstrj(const char* name, float border, bool sep_files);
int cluster_analysis(std::vector<int> tp, float border);


int main( int argc, char* argv[] ) {

  if ( argc < 5 ) {
    cout << "Usage: ./postproc-lammpstrj [input.lammpstrj] [first frame] [last frame] [calc_type]..." << endl;
    exit(1);
  }

  string calc_type(argv[4]);

  fr1 = stoi(argv[2]);
  fr2 = stoi(argv[3]);

  frs = read_lammpstrj(argv[1], fr1, fr2);

  if ( frs != ( fr2 - fr1 ) ) {
    cout << "Mismatch in frames read and input!" << endl;
    exit(1);
  }
  
  cout << frs << " frames read" << endl;
  nlist_init();


  if ( calc_type == "RDF" ) {

  /*#include "nl_globals.h"

    void make_nlist(int, vector<vector<double>>, vector<double>, double);
    make_nlist(nsites, xt[2], L[2], 2.75);
    for ( int i=0 ; i<2500 ; i+=500 )
      cout << i << " " << neigh_ct[i] << endl;
    exit(1);
*/
    if ( argc < 8 ) {
        cout << "Usage: postproc-lammpstrj [input.lammpstrj] [first frame index] [last frame index] ";
        cout << "RDF [dr_bin] [type1] [type2]" << endl;
        exit(1);
    }

    string tp1(argv[6]);
    string tp2(argv[7]);
    calc_rdf( stod(argv[5]), tp1, tp2 );

  } // RDF calculation

  else if ( calc_type == "MSD" ) {
    cout << "MSD optional arguments: " << endl;
    cout << "sitemax [integer], max number of sites to use in MSD, ignoring some sites at the end of the position array" << endl;
    cout << "type [integer], particle type on which to calculate MSD." << endl;


    // Set default list of sites to calculate MSD for
    int ns_msd = nsites;

    vector<int> site_list(ns_msd,0);
    for ( int i=0 ; i<ns_msd ; i++ ) {
      site_list[i] = i;
    }



    if ( argc > 5 ) {
      string opts(argv[5]);
      if ( opts == "sitemax" ) {
        ns_msd = atoi(argv[6]);
      }

      else if ( opts == "type" ) {
        if ( !has_type ) {
          cout << "ERROR: TYPE NOT INCLUDED IN LAMMPSTRJ FILE" << endl;
          exit(1);
        }

        ns_msd = 0;
        int msd_typ = atoi(argv[6]);
        for ( int i=0 ; i<nsites; i++ ) {
          if (type[i] == msd_typ) {
            site_list[ns_msd] = i;
            ns_msd++;
          }
        }
        cout << "Found " << ns_msd << " sites of type " << msd_typ << endl;
        if ( ns_msd == 0 ) {
          cout << "ERROR: FOUND NO SITES OF THE GIVEN TYPE!" << endl;
          exit(1);
        }
      }// opts == type
    }// Extra arguments (args = 6)
    

    remove_pbc_time_jumps(xt, nsites, frs, L);

    calc_msd(xt, site_list, ns_msd, frs, L);
  } // MSD calculation

  else if ( calc_type == "VAN-HOVE" ) {
    if ( argc < 7 ) {
      cout << "Usage: VAN-HOVE [delt] [dx_bin] [optional: site_max]" << endl;
      exit(1);
    }

    double dx_bin = atof( argv[6] );
    int delt = atoi( argv[5] );

    int sitemax = nsites ;
    if ( argc == 8 )
      sitemax = atoi( argv[7] );

    van_hove( xt, sitemax, frs, delt, dx_bin, L);

  }

  else if ( calc_type == "LC_ORDER" ) {
    if ( argc < 7 ) {
      cout << "Usage: postproc-lammpstrj [input.lammpstrj] [first frame index] [last frame index] ";
      cout << "LC_ORDER [LC type] [sites per LC]" << endl;
      cout << "NOTE: the LC type is the type listed in the lammpstrj file, it is not shifted to be zero indexed." << endl;
      exit(1);
    }

    if ( !has_type || !has_mol ) {
      cout << "Requires both molecule id and atom type in input file!" << endl;
      return 1 ;
    }
 
    connect_molecules(xt, mol, nsites, frs, L ) ;
 
    int lc_type = stoi( argv[5] ) ;
    int per_lc = stoi(argv[6] ) ;
 
    lc_order(xt, nsites, frs, type, lc_type, per_lc ) ;
  } // LC_ORDER calculation
  else if ( calc_type == "FT_ANALYSIS" ) {
    #include "mesh_globals.h"
    if ( argc < 8 ) {
      cout << "Usage: postproc-lammpstrj [input.lammpstrj] [first frame index] [last frame index] ";
      cout << "FT_ANALYSIS PER_FRAME_SQ_FLAG[true/false] [Nx] [Ny] [Nz]\n";
      cout << "Required arguments:\n";
      cout << "  [type1] [type2] or \"all\"\n";
      cout << "    all outputs a grid denisty for each type\n";
      cout << "Optional arguments:\n";
      cout << "  k2_cutoff [cutoff] \n" << endl;
      exit(1);
    }

    if (std::string(argv[5]) == "false" || ("0" == std::string(argv[5])))
      per_frame_sq_flag = false;
    else if (std::string(argv[5]) == "true" || ("1" == std::string(argv[5])))
      per_frame_sq_flag = true;
    else {
      cout << "PER_FRAME_SQ_FLAG must be either false, true, 0 or 1" << endl;
      exit(1);
    }

    Nx[0] = atoi(argv[6]);
    Nx[1] = atoi(argv[7]);
    Nx[2] = atoi(argv[8]);
    ML = Nx[0] * Nx[1] * Nx[2];
    M = ML;

    MPI_Init(&argc,&argv);
    k2_cutoff = 0.0f;

    ntypes = *max_element(std::begin(type), std::end(type)); 
    for (int i = 9; i <argc; i++){
      if (strcmp("all", argv[i]) == 0) {
        unique_types = type;
      } else if (string(argv[i]) == "k2_cutoff"){
          k2_cutoff = atof(argv[i+1]);
          i += 1;
          if (k2_cutoff < 0){
            cout << "Cannot have negative cutoff"<<endl;
            exit(1);
          }
      } else if (atoi(argv[i]) > ntypes){
          cout << "Please make sure the selected molecule type exists in the dump file."<< endl;
          exit(1);
      } else if (atoi(argv[i]) < 0){
          cout << "Please enter a non-negative molecule type." << endl;
          exit(1);
      } else {
          unique_types.push_back(atoi(argv[i]));
        }
    }

    std::sort(unique_types.begin(), unique_types.end());
    vector<int>::iterator ip = std::unique(unique_types.begin(), unique_types.end());
    unique_types.resize(std::distance(unique_types.begin(), ip));

    sq_routine();

    MPI_Finalize();

  } else if ( calc_type == "TRIM_TRAJ" ) {
    
    if ( argc < 7 ) {
      cout << "Usage: postproc-lammpstrj [input.lammpstrj] [first frame index] [last frame index] ";
      cout << "TRIM_TRAJ [output basename]\n" << endl;

      cout << "Required arguments:\n";
      cout << "  border [distance] \n";
      cout << "Optional arguments:\n";
      cout << "  [cutoff] \n" << endl;
      exit(1);
    }
    std::string output_file = string(argv[5]);

    bool separate_files = false;
    float border = 0.0f;
    for (int i = 6; i <argc; i++){
      if (string(argv[i]) == "border"){
        border = atof(argv[i+1]);
        i += 1;
      } else if (string(argv[i]) == "sep_files_flag"){
        i++;
        if (std::string(argv[i]) == "false" || ("0" == std::string(argv[i])))
          separate_files = false;
        else if (std::string(argv[i]) == "true" || ("1" == std::string(argv[i])))
          separate_files = true;
        else {
          cout << "sep_files_flag must be either false, true, 0 or 1" << endl;
          exit(1);
        }
        i += 1;
      }
    }

    if (border == 0.0f) cout << "Warning: Border not defined." <<endl;
    trim_lammpstrj(output_file.c_str(), border,separate_files);

  } else if ( calc_type == "CLUSTER") {
    if ( argc < 7 ) {
      cout << "Usage: postproc-lammpstrj [input.lammpstrj] [first frame index] [last frame index] ";
      cout << "CLUSTER " << endl;
      cout << "Required arguments:\n";
      cout << "  [type1] [type2] or \"all\"\n";
      cout << "    Cluster analysis is based on all types used specified\n";
      cout << "  cutoff [distance] or \"all\"\n";
      exit(1);
    }

    connect_molecules(xt, mol, nsites, frs, L ) ;

    float border=0.0f;

    ntypes = *max_element(std::begin(type), std::end(type)); 
    for (int i = 5; i <argc; i++){
      if (strcmp("all", argv[i]) == 0) {
        unique_types = type;
      } else if (string(argv[i]) == "cutoff"){
        border = atof(argv[i+1]);
        i += 1;
      } else if (atoi(argv[i]) > ntypes){
          cout << "Please make sure the selected molecule type exists in the dump file."<< endl;
          exit(1);
      } else if (atoi(argv[i]) < 0){
          cout << "Please enter a non-negative molecule type." << endl;
          exit(1);
      } else {
          unique_types.push_back(atoi(argv[i]));
        }
    }

    std::sort(unique_types.begin(), unique_types.end());
    vector<int>::iterator ip = std::unique(unique_types.begin(), unique_types.end());
    unique_types.resize(std::distance(unique_types.begin(), ip));

    cluster_analysis(unique_types, border);

  } else {
    cout << "Not a valid analysis type.\n";
    cout << "Please check your spelling or consider implementing this feature into the code." << endl;
    exit(1);
  }

  return 0;
}

