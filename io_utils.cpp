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
  // sprintf( nm, "%s.p.dat" , lbl ) ;

  sprintf( nm, "data/%s.%d.dat" , lbl,t ) ;
  char *mm;
  mm = (char*) nm;

  // char *mode = new char[100];
  // sprintf(mode,"w");
  char *mode= strdup("w");
  char *hehe;
  hehe = (char*) "w";
  otp = fopen_mkdir(mm,hehe);




  // fprintf( otp ,"step %d ,post_spin_dt= %lf\n", t, (step-sample_wait)*delt );

  // cout<< i <<'\t'<< nn[0]<<'\t'<<nn[1]<<'\t'<<nn[2]<<"\t"<<M<<endl;
  for ( i=1 ; i<M ; i++ ) {
    me_unstack( i , n) ;

  // cout<< i <<'\t'<< nn[0]<<'\t'<<nn[1]<<'\t'<<nn[2]<<"\t"<<M<<endl;
    k2 = me_get_k( i , kv ) ;
    if (k2 >= 10) continue;
  // cout<< i <<'\t'<< nn[0]<<'\t'<<nn[1]<<'\t'<<nn[2]<<"\t"<<M<<endl;
    // cout<< i<<"\t"<<kv<<'\t'<<k2<<endl;
    // cout<<i<<endl;
    // cout<<k2<<endl;

    // for ( j=0 ; j<Dim ; j++ ) 
    //   fprintf( otp , "%lf " , kv[j] ) ;

    fprintf( otp , "%d %1.5e %1.5e %1.5e %1.5e\n" , t, abs(kdt[i]), sqrt(k2), 
       real(kdt[i]) , imag(kdt[i]) ) ;
    // fprintf( otp , "%1.5e %1.5e %1.5e\n" ,abs(kdt[i]), (k2), 
    //     real(kdt[i]) ) ;
    //if ( Dim == 2 && nn[0] == Nx[0]-1 )
      //fprintf( otp , "\n" ) ;
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
    
    // if ( Dim == 2 && nn[0] == Nx[0]-1 )
    //   fprintf( otp , "\n" ) ;
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
