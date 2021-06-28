#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
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


int main( const int argc, const char* argv[] ) {

  if ( argc < 4 ) {
    cout << "Usage: ./postproc-lammpstrj [input.lammpstr] [first frame] [last frame] [calc_type]..." << endl;
    exit(1);
  }

  int fr1 = stoi(argv[2]);
  int fr2 = stoi(argv[3]);

  int frs = read_lammpstrj(argv[1], fr1, fr2);

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

  return 0;
}

