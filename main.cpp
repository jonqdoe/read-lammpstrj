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
void connect_molecules(vector<vector<vector<double>>>&, vector<int>, int, int, vector<vector<double>>) ;
void lc_order(vector<vector<vector<double>>>, int, int, vector<int>, int, int);

int main( const int argc, const char* argv[] ) {

  if ( argc < 6 ) {
    cout << "Usage: postproc-lammpstrj [input.lammpstrj] [first frame index] [last frame index] [LC type] [sites per LC]" << endl;
    cout << "NOTE: the LC type is the type listed in the lammpstrj file, it is not shifted to be zero indexed." << endl;
    exit(1);
  }

  int fr1 = stoi(argv[2]);
  int fr2 = stoi(argv[3]);

  int frs = read_lammpstrj(argv[1], fr1, fr2);

  if ( frs != ( fr2 - fr1 ) ) {
    cout << "Mismatch in frames read and input!" << endl;
    exit(1);
  }


  if ( !has_type || !has_mol ) {
    cout << "REquires both molecule id and atom type in input file!" << endl;
    return 1 ;
  }

  connect_molecules(xt, mol, nsites, frs, L ) ;

  int lc_type = stoi( argv[4] ) ;
  int per_lc = stoi(argv[5] ) ;

  lc_order(xt, nsites, frs, type, lc_type, per_lc ) ;

  return 0;
}

