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
void calc_rdf( double, string, string ) ;

int main( const int argc, const char* argv[] ) {

  if ( argc < 7 ) {
    cout << "Usage: postproc-lammpstrj [input.lammpstrj] [first frame index] [last frame index] [dr_bin] [type1] [type2]" << endl;
    exit(1);
  }

  int fr1 = stoi(argv[2]);
  int fr2 = stoi(argv[3]);

  nframes = read_lammpstrj(argv[1], fr1, fr2);

  if ( nframes != ( fr2 - fr1 ) ) {
    cout << "Mismatch in frames read and input!" << endl;
    exit(1);
  }

  string tp1(argv[5]);
  string tp2(argv[6]);
  calc_rdf( stod(argv[4]), tp1, tp2 );

  return 0;
}

