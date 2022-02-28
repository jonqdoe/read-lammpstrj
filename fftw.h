#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <fftw3.h>
#include <fftw3-mpi.h>
using namespace std ;

#ifndef FFTW
extern
#endif
int zstart, ML;

#ifndef FFTW
extern
#endif
fftw_complex *fmin0, *fmot0 ;

#ifndef FFTW
extern
#endif
fftw_plan fwd0, fbk0 ;

#ifndef FFTW
extern
#endif
complex<double> I;


#ifndef FFTW_H
#define FFTW_H
int fft_init( int Nx[3] );

void fftw_fwd(const double* in, complex<double>* out, int M);
#endif