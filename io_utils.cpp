// #include "globals.h"

#include <sys/stat.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "io_utils.h"
#include "mesh_globals.h"

extern int t;
FILE *fopen_mkdir(char*, char*);

void write_kspace_data( const char *lbl , std::complex<double> *kdt ) {
  int i, j , n[3] ;
  FILE *otp ;
  double kv[3], k2 ;
  
  char nm[80] ;

  sprintf( nm, "data/%s.dat" , lbl ) ;
  char *mm;
  mm = (char*) nm;

  char *mode= strdup("w");
  char *hehe;
  hehe = (char*) "w";
  otp = fopen_mkdir(mm,hehe);

  for ( i=1 ; i<M ; i++ ) {
    me_unstack( i , n) ;

    k2 = me_get_k( i , kv ) ;
    if ((k2 >= k2_cutoff) && (k2_cutoff > 0) ) continue; 

    fprintf( otp , "%d %1.5e %1.5e %1.5e %1.5e\n" , t, abs(kdt[i]), sqrt(k2), 
       real(kdt[i]) , imag(kdt[i]) ) ;
  }

  fclose( otp ) ;


}

void write_grid_data( const char *nm , double *dat ) {

  int i, j, nn[3] ;
  FILE *otp ;
  double r[3] ;
  otp = fopen( nm , "w" ) ;
  double dx[3];

  for ( i=0 ; i<M ; i++ ) {
    me_unstack( i , nn ) ;
    
    for ( j=0 ; j<3; j++ )
      fprintf( otp , "%lf " , double(nn[j]) * dx[j] ) ;
    
    fprintf( otp , "%1.16e \n" , dat[i] ) ;
  }

  fclose( otp ) ;
}


void rek_mkdir(char *path)
{
  char *sep = strrchr(path, '/' );
  if(sep != NULL) {
    *sep = 0;
    rek_mkdir(path);
    *sep = '/';
  }
  if( mkdir(path,0755) && errno != EEXIST )
    printf("error while trying to create '%s'\n%m\n",path ); 
}


FILE *fopen_mkdir( char *path, char *mode )
{
    char *sep = strrchr(path, '/' );
    if(sep ) { 
       char *path0 = strdup(path);
       path0[ sep - path ] = 0;
       rek_mkdir(path0);
       free(path0);
    } 
    return fopen(path,mode);
}
