#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
using namespace std;

#include "log_space.h"

void log_space(int,int,int);
double pbc_mdr2(vector<double>, vector<double>, vector<double>, vector<double>);

void calc_msd( vector<vector<vector<double>>> xt, int nsites, int nframes, vector<vector<double>> L ) {

  log_space(1,nframes,nframes/2) ;

  vector<double> msd(n_times,0), dr(3);

  for ( int ti = 0; ti < n_times ; ti++ ) {

    int delt = times[ti] ;

    for ( int t=0 ; t<nframes-delt ; t++ ) {

      for ( int i=0 ; i<nsites ; i++ ) {

        double mdr2 = pbc_mdr2(xt[t][i], xt[t+delt][i], dr, L[t] ) ;

        msd[ti] += mdr2 ;

      }// i=0:nsites

    }// t=0:nframes-delt

  }// ti=0:n_times

  ofstream otp("msd.dat");

  for ( int ti=0 ; ti<n_times ; ti++ ) {
    int delt = times[ti] ;
    
    msd[ti] *= ( 1.0 / double(nsites * (nframes-delt) ) ) ;

    otp << ti << " " << msd[ti] << endl;

  }

  otp.close() ;
}
