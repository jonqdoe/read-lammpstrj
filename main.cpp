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
void remove_pbc_time_jumps(vector<vector<vector<double>>>, int, int, vector<vector<double>>) ;
void calc_msd(vector<vector<vector<double>>>, const int, const int, vector<vector<double>> );


int main( const int argc, const char* argv[] ) {

  if ( argc < 4 ) {
    cout << "Usage: postproc-lammpstrj [input.lammpstrj] [first frame index] [last frame index]" << endl;
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

  remove_pbc_time_jumps(xt, nsites, frs, L);

  calc_msd(xt, nsites, frs, L);

  return 0;
}

