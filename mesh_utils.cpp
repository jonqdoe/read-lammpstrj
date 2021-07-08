#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std ;

#define MESHUTILS
#define Dim 3
#include "mesh_globals.h"

void lagrange_get_weights( double , double , double* ) ;
void spline_get_weights( double , double , double* ) ;
void add_charge( int ) ;
void charge_grid(void);


int t;

void sq_routine(){

    fft_init(Nx);

    allocate_grid_variables();

    for (auto type : unique_mol_id){

      for (auto &x: tmp1){ x = 0.0f; }
      for (auto &x: tmp2){ x = 0.0f; }

      for (t=0; t<frs; t++) {
        printf("\rTimestep %d out of %d", t, frs); fflush( stdout );

        V = 1; 
        vector<double> &Lloc = L.at(t);
        for (int j=0 ; j<Dim ; j++ ) {
          V *= Lloc.at(j) ;
          dx[j] = Lloc.at(j)/Nx[j]; 
          Lh[j] = 0.5 * Lloc[j] ;
        }
        gvol = V / double( M ) ;
        
        charge_grid();


        const vector<double> &tmp_rho = rho.at(type);

        fftw_fwd( tmp_rho.data(), tmp1.data(), M ) ;

        for (int i=0 ; i<M ; i++ ){
          tmp1.at(i) = tmp1.at(i) * conj(tmp1.at(i)) ;
        }
        
        if (per_frame_sq_flag == true){
          std::string strdata = "sq.time"+to_string(t+fr1)+".type"+to_string(type);
          write_kspace_data(strdata.c_str(), tmp1.data());
        }
        for (int i = 0; i < tmp2.size(); i++){
        tmp2.at(i) += tmp1.at(i);
        }
      }
      
        for (int i = 0; i < tmp2.size(); i++){
          tmp2.at(i) /= frs;
        }
      std::string strdata = "sq,avg."+to_string(frs)+".type"+to_string(type);
      write_kspace_data(strdata.c_str(), tmp2.data());
    }
}

void charge_grid(void){

  int i, j;
  
  // for ( i=0 ; i<M ; i++ ) 
  //   for ( j=0 ; j<ntypes ; j++ ) 
  //     rho[j][i] = 0.0 ;
  for (vector<double> &x: rho)
    fill(x.begin(), x.end(), 0.0f);


  ////////////////////////////////////////////////////
  // Add segments to each processors density fields //
  ////////////////////////////////////////////////////

  int id=9;
  j=0;
  int t=0;
  // cout<<x[t][id][0]<<'\t'<<dx[j]<<endl;
  for ( i=0 ; i<nsites; i++ ) {
    if ( type.at(i) != -1 ){
     add_segment( i ) ;
     }
  }

}

///////////////////////////////////
// Add all particles to the grid //
///////////////////////////////////
// void fill_grid( ) {
//   int i, id;
//   for ( i=0 ; i<M ; i++ ) 
//     rhoc[i] = 0.0 ;
//   for ( i=0 ; i<ncharge ; i++ )  {
//     id = chg_list[i] ;
//     add_particle( id ) ;
//   }
// }


///////////////////////////////////////////////////////////////////
// Taken directly from appendix in Petersen JCP V103 3668 (1995) //
///////////////////////////////////////////////////////////////////
void lagrange_get_weights( double dx , double H , double *W ) {
  double norm, dx2, dx4, H2, H4;
  
  if ( pmeorder == 0 ) 
    W[0] = 1. ;

  else if ( pmeorder == 2 ) {
    norm = 2.*H*H ;
    dx2 = dx * dx ;
    W[0] = (dx2 - H * dx ) / norm ;
    W[1] = ( -2.*dx2 + 2*H*H ) / norm ;
    W[2] = (dx2 + H * dx ) / norm ;
  }

  else if ( pmeorder == 3 ) {
    norm = 48.*H*H*H ;
    dx2 = dx * dx ;
    W[0] = ( -8.*dx2*dx + 12.*H*dx2 + 2.*H*H*dx - 3.*H*H*H ) / norm ;
    W[1] = ( 24.*dx2*dx - 12.*H*dx2 - 54.*H*H*dx + 27.*H*H*H ) / norm ;
    W[2] = ( -24.*dx2*dx - 12.*H*dx2 + 54.*H*H*dx + 27.*H*H*H ) / norm ;
    W[3] = ( 8.*dx2*dx + 12.*H*dx2 - 2.*H*H*dx - 3.*H*H*H ) / norm ;
  }

  else if ( pmeorder == 4 ) {
    H2 = H * H;
    H4 = H2 * H2 ;
    dx2 = dx*dx;
    dx4 = dx2*dx2;

    norm = 24.0 * H4;

    W[0] = ( dx4 - 2.*H*dx2*dx - H2*dx2 + 2*H2*H*dx ) / norm ;
    W[1] = ( -4.*dx4 + 4.*H*dx2*dx + 16.*H2*dx2 - 16.*H2*H*dx ) / norm ;
    W[2] = ( 6.*dx4 - 30.*H2*dx2 + 24.*H4 ) / norm ;
    W[3] = ( -4.*dx4 - 4.*H*dx2*dx + 16.*H2*dx2 + 16.*H2*H*dx ) / norm ;
    W[4] = ( dx4 + 2.*H*dx2*dx - H2*dx2 - 2*H2*H*dx ) / norm ;

  }

  else if ( pmeorder == 5 ) {
    // I think there is a bug in here! //
    cout << "Possible bug in lagrange interpolation with pmeorder==5\n" ;
    H2 = H * H ;
    H4 = H2 * H2 ;

    dx2 = dx * dx;
    dx4 = dx2 * dx2 ;

    norm = 3840. * H4 * H;

    W[0] = ( -32.*dx4*dx + 80.*H*dx4 + 80.*H2*dx2*dx - 200.*H2*H*dx2
             - 18.*H4*dx + 45.*H4*H ) / norm ;
    W[1] = ( 160.*dx4* - 240.*H*dx4 - 1040.*H2*dx2*dx + 1560.*H2*H*dx2
           + 250.*H4*dx - 375.*H4*H ) / norm ;
    W[2] = ( -320.*dx4*dx + 160.*H*dx4 + 2720.*H2*dx2*dx - 1360.*H2*H*dx2
           - 4500.*H4*dx + 2250*H4*H ) / norm ;

    W[2] = ( 320.*dx4*dx + 160.*H*dx4 - 2720.*H2*dx2*dx - 1360.*H2*H*dx2
           + 4500.*H4*dx + 2250*H4*H ) / norm ;
    W[1] = ( -160.*dx4* - 240.*H*dx4 + 1040.*H2*dx2*dx + 1560.*H2*H*dx2
           - 250.*H4*dx - 375.*H4*H ) / norm ;
    W[0] = ( 32.*dx4*dx + 80.*H*dx4 - 80.*H2*dx2*dx - 200.*H2*H*dx2
             + 18.*H4*dx + 45.*H4*H ) / norm ;
  }


  else {
    cout << "PME order is " << pmeorder << endl;
    cout << "PME not set up for this interpolation order!" << endl;
    exit(1);
  }

}


void spline_get_weights( double dx , double H , double *W ) {
  double sx = dx / H ;
  
  double sx2, sx3, sx4, sx5, sx6;

  if ( pmeorder == 0 ) 
    W[0] = 1. ;

  else if ( pmeorder == 1 ) {
    W[0] = ( 1. - 2. * sx ) / 2.0 ;
    W[1] = ( 1. + 2. * sx ) / 2.0 ;
  }

  else if ( pmeorder == 2 ) {
    sx2 = sx * sx ;

    W[0] = (1. - 4. * sx + 4.*sx2 ) / 8. ;
    W[1] = (3. - 4. * sx2 ) / 4. ;
    W[2] = (1. + 4. * sx + 4.*sx2 ) / 8. ;

  }

  else if ( pmeorder == 3 ) {
    sx2 = sx * sx ;
    sx3 = sx2 * sx ;

    W[0] = ( 1. - 6.*sx + 12.*sx2 - 8.*sx3 ) / 48. ;
    W[1] = ( 23. - 30.*sx - 12.*sx2 + 24.*sx3 ) / 48. ;
    W[2] = ( 23. + 30.*sx - 12.*sx2 - 24.*sx3 ) / 48. ;
    W[3] = ( 1. + 6.*sx + 12.*sx2 + 8.*sx3 ) / 48. ;
  }

  else if ( pmeorder == 4 ) {
    sx2 = sx * sx ;
    sx3 = sx2 * sx2 ;
    sx4 = sx2 * sx2 ;

    W[0] = (1. - 8.*sx + 24.*sx2 - 32.*sx3 + 16.*sx4 ) / 384.;
    W[1] = (19. - 44.*sx + 24.*sx2 + 16.*sx3 - 16.*sx4 ) / 96. ;
    W[2] = (115. - 120.*sx2 + 48.*sx4 ) / 192. ;
    W[3] = (19. + 44.*sx + 24.*sx2 - 16.*sx3 - 16.*sx4 ) / 96. ;
    W[4] = (1. + 8.*sx + 24.*sx2 + 32.*sx3 + 16.*sx4 ) / 384.;
  }

  else if ( pmeorder == 5 ) {
    sx2 = sx * sx ;
    sx3 = sx2 * sx2 ;
    sx4 = sx2 * sx2 ;
    sx5 = sx4 * sx ;

    W[0] = (1. - 10.*sx + 40.*sx2 - 80.*sx3 + 80.*sx4 - 32.*sx5 ) / 3840.;
    W[1] = (237. - 750.*sx + 840.*sx2 - 240.*sx3 - 240.*sx4 + 160.*sx5 ) / 3840.;
    W[2] = (841. - 770.*sx - 440.*sx2 + 560.*sx3 + 80.*sx4 - 160.*sx5 ) / 1920. ;
    W[3] = (841. + 770.*sx - 440.*sx2 - 560.*sx3 + 80.*sx4 + 160.*sx5 ) / 1920. ;
    W[4] = (237. + 750.*sx + 840.*sx2 + 240.*sx3 - 240.*sx4 - 160.*sx5 ) / 3840.;
    W[5] = (1. + 10.*sx + 40.*sx2 + 80.*sx3 + 80.*sx4 + 32.*sx5 ) / 3840.;

  }


  else
    {
      cout << "P3M not set up for this interpolation order!" <<endl;
      exit(10);
    }

}




// Stacks vector x into 1D array index in [ 0, M-1 ]
int me_stack(int x[3]) {
  return  (x[0] + (x[1] + x[2]*Nx[1])*Nx[0] );
}



// Receives index id in [0 , M[b] ] and makes array
// nn[Dim] in [ 0 , Nx[Dim] ]
void me_unstack(int id, int nn[3]) {

  nn[2] = id/Nx[1]/Nx[0];

  nn[1] = id/Nx[0] - nn[2]*Nx[1];
  
  nn[0] = id - (nn[1] + nn[2]*Nx[1])*Nx[0];
  
}



double me_get_k(int id, double *k ) {

  double kmag = 0.0;

  int i, n[3];
  if (t==frs) t--;
  const vector<double> &box = L.at(t);

  me_unstack(id, n);

  if ( double(n[0]) < double(Nx[0]) / 2. )
   k[0] = 2*PI*double(n[0])/box[0];
  else
   k[0] = 2*PI*double(n[0]-Nx[0])/box[0];

  if ( double(n[1]) < double(Nx[1]) / 2. )
   k[1] = 2*PI*double(n[1])/box[1];
  else
   k[1] = 2*PI*double(n[1]-Nx[1])/box[1];

  if ( double(n[2]) < double(Nx[2]) / 2. )
    k[2] = 2*PI*double(n[2])/box[2];
  else
    k[2] = 2*PI*double(n[2]-Nx[2])/box[2];

  for (i=0; i<3; i++)
    kmag += k[i]*k[i];

  return kmag;

}




//////////////////////////////////////////////
// Adds the segment associated with particle //
// "id" to the PME grid using Lagrange      //
// interpolation scheme. From JCP V103 3668 //
//////////////////////////////////////////////
void add_segment( int id ) {


  
  int j, g_ind[Dim] , ix, iy, iz, nn[Dim] , Mindex, grid_ct; 
  //double **W , gdx , W3;
  
  // double gdx , W3;

  double **W , gdx , W3;
  
  W = ( double** ) calloc( Dim , sizeof( double* ) );

  std::vector<std::vector<double>> &x = xt.at(t);

  // std::cout<<t<<std::endl;
  // std::cout<<L[t][0]<<std::endl;
  // std::cout<<L[t][1]<<std::endl;
  // std::cout<<L[t][2]<<std::endl;

  ///////////////////////////////////////////////
  // First, determine the relevant weights for //
  // all grid points in all directions.        //
  ///////////////////////////////////////////////

  for ( j=0 ; j<Dim ; j++ ) {

   W[j] = ( double* ) calloc( pmeorder+1 , sizeof( double ) );

    // Distance to nearest grid point if even //
    if ( pmeorder % 2 == 0 ) {
      g_ind[j] = int( ( x.at(id).at(j) + 0.5 * dx[j] ) / dx[j] ) ;

      gdx = x.at(id).at(j) - double( g_ind[j] ) * dx[j] ;
    }
 

    // Distance to nearest mid-point between grid points if odd //
    else {
      g_ind[j] = int( ( x.at(id).at(j)  ) / dx[j] ) ;
      if ( g_ind[j] >= Nx[j] ){
      g_ind[j] = Nx[j]-1;
      }
      gdx = x.at(id).at(j)  - ( double( g_ind[j] ) + 0.5 ) * dx[j] ;
    }


    /////////////////////////////////////////
    // Get the weights for each grid point //
    /////////////////////////////////////////

      spline_get_weights( gdx , dx[j] , W[j] );

  }//for ( j=0 ; j<3...


  ////////////////////////////////////////////////////
  // Assign the weights to all relevant grid points //
  ////////////////////////////////////////////////////
  grid_ct = 0 ;
  
  ///////////////////////////////////////////
  // 3D version of particle-to-mesh scheme //
  ///////////////////////////////////////////
  if ( Dim == 3 ) {
    for ( ix = 0 ; ix < pmeorder+1 ; ix++ ) {
      // cout<<g_ind[0]<<"\t"<<ix<<'\t'<<( pmeorder/2 + pmeorder % 2 )<<endl; 
      nn[0] = g_ind[0] + ix - ( pmeorder/2 + pmeorder % 2 );
      
      if ( nn[0] < 0 ) nn[0] += Nx[0] ;
      else if ( nn[0] >= Nx[0] ) nn[0] -= Nx[0] ;
  
      for ( iy = 0 ; iy < pmeorder+1 ; iy++ ) {
  
        nn[1] = g_ind[1] + iy - ( pmeorder/2 + pmeorder % 2 ) ;
        
        if ( nn[1] < 0 ) nn[1] += Nx[1] ;
        else if ( nn[1] >= Nx[1] ) nn[1] -= Nx[1] ;
  
        for ( iz = 0 ; iz < pmeorder+1 ; iz++ ) {
  
          nn[2] = g_ind[2] + iz - ( pmeorder/2 + pmeorder % 2 ) ;
          
          if ( nn[2] < 0 ) nn[2] += Nx[2] ;
          else if ( nn[2] >= Nx[2] ) nn[2] -= Nx[2] ;
  
  
          // stack() returns index in [0,M]
            Mindex = me_stack( nn ) ;
  
          if ( Mindex >= M ) {
            char nm[80] ;
            sprintf(nm, "%d Index = %d out of range, particle %d %lf %lf %lf\n" , 
                t, Mindex, id , x[id][0], x[id][1], x[id][2] ) ;
            cout << nm << endl;

            // die(nm) ;
            exit(10);
          }
          // if (id==M-2)
        // cout<<"I was playing possu"<<endl;
  
          W3 = W[0][ix] * W[1][iy] * W[2][iz] / gvol ;
  
          
          rho.at(type.at(id)).at(Mindex) += W3; 

  
          grid_inds[ id ][ grid_ct ] = Mindex ;
          grid_W[ id ][ grid_ct ] = W3 ;
  
          grid_ct++ ;
  
        }
      }
    }
  }


  ///////////////////////////
  // Free allocated memory //
  ///////////////////////////
 for ( j=0 ; j<Dim ; j++ ) 
   free( W[j] ) ;

 free(W) ;

}// End lagrange add charge


void allocate_grid_variables(){
  
  pmeorder=1;
  tmp1.resize(M);
  tmp2.resize(M);
  grid_per_particle = pow(pmeorder+1,3);

  grid_W.resize(nsites,std::vector<double>(grid_per_particle));
  rho.resize(ntypes+1,std::vector<double>(M));
  grid_inds.resize(nsites,std::vector<int>(grid_per_particle));

}
