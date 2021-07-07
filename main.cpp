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
void calc_msd(vector<vector<vector<double>>>, const int, const int, vector<vector<double>> );
void calc_rdf(double, string, string);
void nlist_init(void);
void van_hove(vector<vector<vector<double>>>, int, int, int, double, vector<vector<double>> );


int main( int argc, char* argv[] ) {

  if ( argc < 4 ) {
    cout << "Usage: ./postproc-lammpstrj [input.lammpstr] [first frame] [last frame] [calc_type]..." << endl;
    exit(1);
  }

  fr1 = stoi(argv[2]);
  fr2 = stoi(argv[3]);

  frs = read_lammpstrj(argv[1], fr1, fr2);

  if ( frs != ( fr2 - fr1 ) ) {
    cout << "Mismatch in frames read and input!" << endl;
    exit(1);
  }
  
  cout << frs << " frames read" << endl;
  nlist_init();


  string calc_type(argv[4]);


  if ( calc_type == "RDF" ) {
    #include "nl_globals.h"

    void make_nlist(int, vector<vector<double>>, vector<double>, double);
    make_nlist(nsites, xt[2], L[2], 2.75);
    for ( int i=0 ; i<2500 ; i+=500 )
      cout << i << " " << neigh_ct[i] << endl;
    exit(1);

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
    cout << "sitemax [integer], maximum number of sites to use in MSD calculation, ignoring some sites at the end of the position array" << endl;

    if ( argc == 6 )
      nsites = atoi( argv[5] );

    remove_pbc_time_jumps(xt, nsites, frs, L);

    calc_msd(xt, nsites, frs, L);
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
      cout << "FT_ANALYSIS PER_FRAME_SQ_FLAG [Nx] [Ny] [Nz] [type1] [type2]... \n" << endl;
      exit(1);
    }

    if (std::string(argv[5]) == "false" || ("0" == std::string(argv[5])))
      per_frame_sq_flag = false;
    else if (std::string(argv[5]) == "true" || ("1" == std::string(argv[5])))
      per_frame_sq_flag = true;
    else {
      cout << "PER_FRAME_SQ_FLAG must be either FALSE, TRUE, 0 or 1" << endl;
      exit(1);
    }

    Nx[0] = atoi(argv[6]);
    Nx[1] = atoi(argv[7]);
    Nx[2] = atoi(argv[8]);
    ML = Nx[0] * Nx[1] * Nx[2];
    M = ML;

    MPI_Init(&argc,&argv);

    ntypes = *max_element(std::begin(type), std::end(type)); 
    if (strcmp("all", argv[9]) == 0) {
      unique_mol_id = type;
    } else {
        for (int i = 9; i <argc; i++){
          if (atoi(argv[i]) > ntypes){
            cout << "Please make sure the selected molecule type exists in the dump file."<< endl;
            exit(1);
          }
          else if (atoi(argv[i]) < 0){
            cout << "Please enter a non-negative molecule type." << endl;
            exit(1);
          }
          else {
            unique_mol_id.push_back(atoi(argv[i]));
          }
      }
    }
    std::sort(unique_mol_id.begin(), unique_mol_id.end());
    vector<int>::iterator ip = std::unique(unique_mol_id.begin(), unique_mol_id.end());
    unique_mol_id.resize(std::distance(unique_mol_id.begin(), ip));

    sq_routine();

    MPI_Finalize();

  }
  else {
    cout << "Not a valid analysis command.\n";
    cout << "Please check your spelling or consider implementing this feature into the code." << endl;
    exit(1);
  }

  return 0;
}

