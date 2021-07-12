#define FFTW
#include "fftw.h"

// Forward transform 
void fftw_fwd(const double* in, complex<double>* out, int M) {

  int i;

  // Store fft input
  for (i=0; i<M; i++) {
    fmin0[i][0] = double(in[i]);
    fmin0[i][1] = 0.0 ;
  }

  fftw_execute(fwd0);

  double norm = 1.0 / double(M);

  // Store fft output
  for (i=0; i<M; i++) 
    out[i] =( fmot0[i][0] + I * fmot0[i][1] ) * norm ;

}



// Backwards transform and normalization
void fftw_back(complex<double>* in, double* out, int M ) {

  int i;
  // Store input
  for (i=0; i<M; i++) {
    fmin0[i][0] = real(in[i]);
    fmin0[i][1] = imag(in[i]);
  }


  // Perform fft
  fftw_execute(fbk0);
 

  // Store output
  for (i=0; i<M; i++)
    out[i] = fmot0[i][0] ;
  
}

  
int fft_init( int Nx[3] ) {

  int Dim = 3;
  int NxL[3];

  int i, b , total_alloced = 0 ;
  I = complex<double>(0.0,1.0);


  for ( i=0 ; i<3 ; i++ ) 
    NxL[i] = Nx[i] ;

  ptrdiff_t Dm = Dim, Nfp[3], NxLtp, ztp;
  for (i=0; i<Dim; i++)
    Nfp[i] = Nx[Dim-1-i];
 
  int size = fftw_mpi_local_size_many(Dm, Nfp, 1, 0, MPI_COMM_WORLD,
            &NxLtp, &ztp );
 
  NxL[Dim-1] = NxLtp;
  for (i=0; i<Dim-1; i++)
    NxL[i] = Nx[Dim-i-1];
 
  zstart = ztp;
  
  fmin0 = (fftw_complex*) fftw_malloc( size * sizeof(fftw_complex) );
  fmot0 = (fftw_complex*) fftw_malloc( size * sizeof(fftw_complex) );
  
  fwd0 = fftw_mpi_plan_dft(Dim, Nfp, fmin0, fmot0,
      MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE );
  fbk0 = fftw_mpi_plan_dft(Dim, Nfp, fmin0, fmot0,
      MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE );

  ML = 1;
  for (i=0; i<Dim; i++)
    ML *= NxL[i];
  
  total_alloced += size*sizeof(fftw_complex)*2 ;

  return total_alloced ;
}
