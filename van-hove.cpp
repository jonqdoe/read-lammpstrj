#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <cmath>
#include "log_space.h"

using namespace std;

double pbc_mdr2(vector<double>, vector<double>, vector<double>, vector<double>);

void van_hove( vector<vector<vector<double>>> xt,  // Positions as function of time
    int nsites,                              // Number of particles to consider
    int nfr,                                 // Number of frames in trajectory
    int delt,                                // Lag time for van Hove
    double dx_bin,                           // Bin size for van Hove
    vector<vector<double>> L ) {                   // Box dimensions with time


  vector<double> dr(3);
  double dx_max = 250.0;
  int nbins = int( dx_max / dx_bin );

  double *gr = new double[nbins];
  double *gx = new double[nbins];
  double *gy = new double[nbins];
  double *gz = new double[nbins];
  double norm = 0.0;


  for ( int t=0 ; t<nfr-delt ; t++ ) {

    for ( int i=0 ; i<nsites ; i++ ) {

      // double mdr2 = pbc_mdr2(xt[t][i], xt[t+delt][i], dr, L[t] ) ;

      double mdr2 = 0.0;
      for ( int j=0 ; j<3 ; j++ ) {
        dr[j] = xt[t+delt][i][j] - xt[t][i][j] ;
        mdr2 += dr[j] * dr[j];
      }

      double mdr = sqrt(mdr2);

      int rbin = int( mdr / dx_bin );
      if ( rbin < nbins ) gr[rbin] += 1.0;

      rbin = int( fabs(dr[0]) / dx_bin );
      if ( rbin < nbins ) gx[rbin] += 1.0;

      rbin = int( fabs(dr[1]) / dx_bin );
      if ( rbin < nbins ) gy[rbin] += 1.0;

      rbin = int( fabs(dr[2]) / dx_bin );
      if ( rbin < nbins ) gz[rbin] += 1.0;

      norm += 1.0;

    }// i=0:nsites

  }// t=0:nfr-delt


  char nm[25];
  sprintf(nm, "Gs-dt%d.dat", delt);
  ofstream otp(nm);

  for ( int i=0 ; i<nbins ; i++ ) {
    double r = double(i) * dx_bin ;

    gr[i] *= 1.0 / norm;
    gx[i] *= 1.0 / norm;
    gy[i] *= 1.0 / norm;
    gz[i] *= 1.0 / norm;

    if ( gr[i] > 0.0 ) 
      otp << r << " " << gr[i] << " " << gx[i] << " " << gy[i] << " " << gz[i] << endl;

  }

  otp.close() ;
}
