#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <cmath>
using namespace std;

#include "lammpstrj.h"

#define PI 3.141592654
#define PI2 6.2831853071

double pbc_mdr2(vector<double>, vector<double>, vector<double>, vector<double>);

void calc_rdf( double drbin, string type1, string type2 ) {

  cout << endl;
  cout << "NOTE: calc_rdf.cpp assumes 3D simulation cell!!!" << endl;
  cout << "      a new code is needed to calculate 2D g(r)"  << endl << endl;

  // Flags for the atom type 
  int tp1 = -1, ns1 ;
  int tp2 = -1, ns2 ;

  // Set up type 1 for g(r)
  if ( type1 == "all" )
    ns1 = nsites ;
  else {
    tp1 = stoi(type1) ;
    ns1 = 0;
    for ( int i=0 ; i<nsites ; i++ ) {
      if ( type[i] == tp1 )
        ns1++;
    }
  }

  // Set up type 2 for g(r)
  if ( type2 == "all" )
    ns2 = nsites ;
  else {
    tp2 = stoi(type2) ;
    ns2 = 0;
    for ( int i=0 ; i<nsites ; i++ ) {
      if ( type[i] == tp2 )
        ns2++;
    }
  }

  cout << "Found " << ns1 << " sites in group 1, " << ns2 << " in group 2" << endl;

  // Find the minimum box dimension over the whole simulation
  double bx_min = 23489101243.4 ;
  double vol ;
  for ( int t=0 ; t<nframes ; t++ ) 
    for ( int j=0 ; j<3 ; j++ )
      if ( L[t][j] < bx_min ) 
        bx_min = L[t][j] ;


  // Number of bins in g(r)
  int n_r_bins = int( bx_min / 2.0 / drbin ) + 1;
  cout << "Allocating " << n_r_bins << " bins for g(r)" << endl;

  // allocate g(r)
  vector<double> gr(n_r_bins);

  vector<double> dr(3);
  double max_dr2 = bx_min * bx_min / 4.0 ;

  // Main calculation loop for g(r)
  for ( int t=0 ; t<nframes ; t++ ) {
    for ( int i=0 ; i<nsites ; i++ ) {

      if ( tp1 >= 0 && type[i] != tp1 ) 
        continue ;

      for ( int j=0 ; j<nsites ; j++ ) {
        if ( tp2 >= 0 && type[j] != tp2 )
          continue ;
        if ( i == j )
          continue ;

        double mdr2 = pbc_mdr2( xt[t][i], xt[t][j], dr, L[t] ) ;
        if ( mdr2 >= max_dr2 )
          continue ;

        double mdr = sqrt(mdr2);
        int r_bin = int( mdr / drbin );

        if ( r_bin < 0 || r_bin >= n_r_bins ) {
          cout << "bin error!" << " r_bin = " << r_bin << endl;
          exit(1);
        }

        gr[r_bin] += 1.0 ;

      }// j=i+1:nsites
    }// i=i:nsites-1
  }// t=0:frs


  double avg_vol = 0.0 ;
  for ( int t=0 ; t<nframes ; t++ ) {
    double vol = 1.0 ;
    for ( int j=0 ; j<3 ; j++ ) 
      vol *= L[t][j] ;
    avg_vol += vol ;
  }
  avg_vol = avg_vol / double(nframes) ;

  double norm1 = avg_vol / double( ns1 * ns2 * nframes ) ;

  for ( int i=0 ; i<n_r_bins ; i++ ) {
    double v1 = double(pow(i,3));
    double v2 = double(pow(i+1,3));

    double shell_vol = (v2 - v1) * drbin*drbin*drbin * 4.0 * PI / 3.0 ;

    //gr[i] = shell_vol ;
    gr[i] *= norm1 / shell_vol ;
  }// i=0:n_r_bins


  ofstream otp("gr.dat");

  for ( int i=0 ; i<n_r_bins-1 ; i++ ) 
    otp << (double(i)+0.5) * drbin << " " << gr[i] << endl;

  otp.close();

}
