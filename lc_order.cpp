#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>


using namespace std ;
vector<vector<double>> make_tensor( vector<double>) ;


void lc_order( vector<vector<vector<double>>> xt, int ns, int frs, vector<int> type, int lc_type, int ns_per_lc) {

  vector<int> lc_inds(ns,-1);

  // Figure out how many LC sites we have
  int ns_tot_lc = 0;

  for ( int i=0 ; i<ns; i++ ) {
    if (type[i] == lc_type) {
      lc_inds[ ns_tot_lc ] = i;
      ns_tot_lc++;
    }
  }


  // Check for consistency
  if ( ns_tot_lc % ns_per_lc != 0 ) {
    cout << "Number of LC sites not divisible by sites per LC molceule. Check input values" << endl;
    exit(1);
  }

  int n_lc_molecs = ns_tot_lc / ns_per_lc ;

  // Allocate unit orientation vector, tensor order parameter
  // Note that S
  vector<vector<double>> u(n_lc_molecs, vector<double>(3)) ;

  cout << "Calculated there are " << n_lc_molecs << " total mesogens in the box" << endl;


  for ( int t=0 ; t<frs ; t++ ) {

    vector<vector<double>> S_avg(3, vector<double>(3,0.0)) ;

    for ( int i=0 ; i<n_lc_molecs ; i++ ) {

      // Find the indices of the two particles at the ends
      // of the current LC
      int i1 = lc_inds[ i * ns_per_lc ];
      int i2 = lc_inds[ (i+1) * ns_per_lc - 1 ];


      // Calculate the vector between the two ends of mesogen
      vector<double> ui(3) ;
      double mdr2 = 0.0;
      for ( int j=0 ; j<3 ; j++ ) {
        ui[j] = xt[t][i1][j] - xt[t][i2][j] ;
        mdr2 += ui[j] * ui[j] ;
      }

      // Normalize it so its a unit vector
      double mdr = sqrt(mdr2);
      for ( int j=0 ; j<3 ; j++ ) 
        ui[j] *= 1.0/mdr ;


      // Declare and calculate S
      vector<vector<double>> S(3, vector<double>(3)) ;
      S = make_tensor(ui);


      // Accumulate the average of S
      for ( int j=0; j<3 ; j++ )
        for ( int k=0 ; k<3 ; k++ )
          S_avg[j][k] += S[j][k] ;


    }// i=0:n_lc_molecs


    cout << "Frame " << t << ":" << endl;
    for ( int j=0 ; j<3 ; j++ ) {
      for ( int k=0 ; k<3 ; k++ ) {
        S_avg[j][k] *= 1.0 / double(n_lc_molecs) ;
        cout << S_avg[j][k] << " " ;
      }
      cout << endl;
    }
    cout << endl;

  }// t=0:frs

}





vector<vector<double>> make_tensor( vector<double> ui ) {

  vector<vector<double>> S(3,vector<double>(3)) ;
  for ( int i=0 ; i<3 ; i++ ) 
    for ( int j=0 ; j<3 ; j++ )
      S[i][j] = ui[i] * ui[j] - ( i==j ? 1.0 : 0.0 ) / 3.0 ;
  return S;

}
