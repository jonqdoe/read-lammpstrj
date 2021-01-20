#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <iostream>

using namespace std ;

#define NLUTILS
#include "nl_globals.h"

void cell_death( const char* ) ;
int cell_stack( vector<int> ) ;
void cell_unstack( int, vector<int> ) ;




// This routine should be called once only //
void nlist_init() {
    
  cell_allocate_flag = 0 ;

}





// Primary nlist routine. Variables:
// ns: number of sites in x
// x[ns][3]: x-y-z position array
// L[3]: Box dimensions
// rc: cut-off to be used for neighbor list.
void make_nlist( int ns,  vector<vector<double>> x, vector<double> L, double rc ) {

  int i, j, k, id2 ;

  // Place particles in the box //
  for ( i=0 ; i<ns ; i++ ) {
    for ( j=0 ; j<3 ; j++ ) {
      while ( x[i][j] > L[j] ) x[i][j] -= L[j] ;
      while ( x[i][j] < 0.0 ) x[i][j] += L[j] ;
    }
  }


  if ( !cell_allocate_flag )
    n_cells.resize(3);


  // Determine cell size //
  tot_cells = 1 ;
  for ( j=0 ; j<3 ; j++ ) {
    n_cells[j] = int( L[j] / rc ) ;

    if ( n_cells[j] < 3 ) 
      n_cells[j] = 3 ;
 
    tot_cells *= n_cells[j] ;

    cell_sz[j] = L[j] / double( n_cells[j] ) ;
  }

  

  
  // Allocate cell memory //
  if ( !cell_allocate_flag ) {
    // Assumes that the maximum local density in a cell is 
    // less than 15 times the box average
    printf("Allocating cell memory..."); fflush(stdout) ;
    CELL_MAX = int( double(ns) / double(tot_cells) * 15.0 ) ;

    part_cells.resize(ns);
    cell_n.resize(tot_cells,0);
    cell_inds.resize(tot_cells, vector<int>(CELL_MAX));

//    part_cells = ( int* ) calloc( ns , sizeof( int ) ) ;
//    cell_n = ( int* ) calloc( tot_cells , sizeof( int ) ) ;
//    cell_inds = ( int** ) calloc( tot_cells, sizeof( int* ) ) ;
//    for ( i=0 ; i<tot_cells ; i++ ) {
//      cell_n[i] = 0 ;
//      cell_inds[i] = ( int* ) calloc( CELL_MAX , sizeof( int ) ) ;
//    }
    printf("done!\n") ; fflush(stdout) ;

    cout << "Allocated roughly " << (ns+tot_cells*(1+CELL_MAX)*sizeof(int))/1.0E6 << "MB" << endl;
  }
  else
    for ( i=0 ; i<tot_cells; i++ )
      cell_n[i] = 0 ;



  // Place particles in cells //
  int c_ind ;
  vector<int> cid(3);
  for ( i=0 ; i<ns ; i++ ) {
    for ( j=0 ; j<3 ; j++ ) {
      cid[j] = int( x[i][j] / cell_sz[j] ) ;
      if ( cid[j] < 0 || cid[j] >= n_cells[j] )
        cell_death("Cell Index Error!") ;
    }


    c_ind = cell_stack( cid ) ;
  
    cell_inds[c_ind][ cell_n[ c_ind ] ] = i ;
    cell_n[ c_ind ] += 1 ;

    if ( cell_n[c_ind] == CELL_MAX )
      cell_death( "Cells at max occupancy" ) ;
    
    part_cells[i] = c_ind ;
  }


  // Allocate neighbor list //
  if ( !cell_allocate_flag ) {
    NEIGH_MAX = CELL_MAX ;
    cout << "Allocating neighbor list memory...NEIGH_MAX=" << NEIGH_MAX << endl;
 
    neigh_ct.resize(ns, 0);
    neigh_ind.resize(ns, vector<int>(NEIGH_MAX));

//    neigh_ct = ( int* ) calloc( ns , sizeof( int ) ) ;
//    neigh_ind = ( int** ) calloc( ns , sizeof( int* ) ) ;
//    for ( i=0 ; i<ns ; i++ ) {
//      neigh_ct[i] = 0 ;
//      neigh_ind[i] = ( int* ) calloc( NEIGH_MAX , sizeof( int ) ) ;
//    }

    cell_allocate_flag = 1 ;
    cout << "Neighborlist allocated!" << endl;
    cout << "Allocated roughly " << (ns+ns*(NEIGH_MAX)*sizeof(int))/1.0E6 << "MB" << endl;
  }

  for ( i=0 ; i<ns ; i++ )
    neigh_ct[i] = 0 ;



  int icell, ncell;

  vector<int> nc(3), ic(3), did(3);
  vector<double> dr(3);

  double mdr2, rc2 = rc * rc ;
  for ( i=0 ; i<ns ; i++ ) {
    icell = part_cells[i] ;

    cell_unstack( icell , ic ) ;

    // Loop over 27 neighbors of i's cells //
    for ( did[0] = -1 ; did[0] <= 1 ; did[0]++ ) {
      for ( did[1] = -1 ; did[1] <= 1 ; did[1]++ ) {
        for ( did[2] = -1 ; did[2] <= 1 ; did[2]++ ) {

          for ( j=0 ; j<3 ; j++ ) {
            nc[j] = ic[j] + did[j] ;
            if ( nc[j] < 0 )
              nc[j] = n_cells[j] - 1 ;
            else if ( nc[j] >= n_cells[j] )
              nc[j] = 0 ;
          }

          ncell = cell_stack( nc ) ;

          // Loop over particles in ncell 
          for ( k=0 ; k<cell_n[ ncell ] ; k++ ) {
            id2 = cell_inds[ ncell ][ k ] ;

            if ( id2 == i )
              continue ;

            mdr2 = 0.0 ;
            for ( j=0 ; j<3 ; j++ ) {
              dr[j] = x[i][j] - x[id2][j] ;
              if ( dr[j] > 0.5 * L[j] ) dr[j] -= L[j] ;
              else if ( dr[j] < -0.5 * L[j] ) dr[j] += L[j] ;

              mdr2 += dr[j] * dr[j] ;
            }

            if ( mdr2 > rc2 )
              continue ;

            //if ( i > 6600 ) 
            //  cout << i << "ct, id2: " << neigh_ct[i] << " " << id2 << endl;
            neigh_ind[i][ neigh_ct[i] ] = id2 ;
            neigh_ct[i] += 1 ;

            if ( neigh_ct[i] >= NEIGH_MAX )
              cell_death( "Neighbor list hit NEIGH_MAX" ) ;

          }// Loop over cell occupants

        }// Loops over neighboring cells
      }
    }

  }// Loop over particles








}


void cell_death( const char *msg ) {
  printf( "%s\n" , msg ) ;
  exit(1);
}

// Stacks vector x into 1D array index in [ 0, tot_cells ]
int cell_stack(vector<int> x) {
  return  (x[0] + (x[1] + x[2]*n_cells[1])*n_cells[0] );
}




// Receives index id in [0 , M[b] ] and makes array
// nn[Dim] in [ 0 , n_cells[Dim] ]
void cell_unstack(int id, vector<int> nn) {

  nn[2] = id/n_cells[1]/n_cells[0];

  nn[1] = id/n_cells[0] - nn[2]*n_cells[1];

  nn[0] = id - (nn[1] + nn[2]*n_cells[1])*n_cells[0];

}

