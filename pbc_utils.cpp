#include <vector>

using namespace std ;
#include "lammpstrj.h"


// Compute rij = ri - rj, returns |rij|^2
double pbc_mdr2( vector<double> ri, vector<double> rj, vector<double> dr, vector<double> L ) {

  double mdr2 = 0.0;

  for ( int j=0 ; j<3 ; j++ ) {
    dr[j] = ri[j] - rj[j] ;
    if ( dr[j] > L[j] ) dr[j] -= L[j] * 0.5 ;
    else if ( dr[j] < -L[j] ) dr[j] += L[j] * 0.5 ;

    mdr2 += dr[j] * dr[j] ;
  }

  return mdr2 ;
}




