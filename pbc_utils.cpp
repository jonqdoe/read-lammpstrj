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


// Puts molecules all in the same periodic image
// This assumes the molecules will not stretch across half the box length
void connect_molecules(vector<vector<vector<double>>> xt, vector<int> mol, 
    int ns, int frs, vector<double> L ) {

  vector<double> dr(3);
  for ( int t=0 ; t<frs ; t++ ) {
    int cur_mol = mol[0] ;
    int first_id = 0 ;

    for ( int i=1 ; i<ns ; i++ ) {
 
      // Check if we're starting a new molecule
      if ( mol[i] != cur_mol ) {
        cur_mol = mol[i] ;
        first_id = i;
      }

      // Not on a new molecule
      else {
        for ( int j=0 ; j<3 ; j++ ) {
          dr[j] = xt[t][i][j] - xt[t][first_id][j] ;
          
          if ( dr[j] > 0.5 * L[j] ) dr[j] -= L[j] ;
          else if ( dr[j] < -0.5 * L[j] ) dr[j] += L[j] ;

          xt[t][i][j] = xt[t][first_id][j] + dr[j] ;
        }
      }
    }// i=0:ns

  }// t=0:frs

}

// Connects particles' trajectories in time //
void remove_pbc_time_jumps( vector<vector<vector<double>>> xt, int ns, int frs, vector<double> L ) {

  vector<double> dr(3);

  for ( int t=0 ; t<frs-1 ; t++ ) {
    for ( int i=0 ; i<ns ; i++ ) {
      for ( int j=0 ; j<3 ; j++ ) {
        dr[j] = xt[t][i][j] - xt[t+1][i][j] ;

        while ( dr[j] > L[j]*0.5 ) dr[j] -= L[j] ;
        while ( dr[j] < -L[j]*0.5) dr[j] += L[j] ;

        xt[t+1][i][j] = xt[t][i][j] - dr[j] ;
      }
    }// i=0:ns
  }// t=0:frs-1

}

