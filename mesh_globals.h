#include <vector>
#include <complex>
#include "fftw.h"
#include "lammpstrj.h"
#include "io_utils.h"

#define PI 3.141592654

#ifndef MESHUTILS
extern
#endif
int pmeorder, Nx[3], M, grid_per_particle ;

#ifndef MESHUTILS
extern
#endif
std::vector<std::complex<double>> tmp1 ;

#ifndef MESHUTILS
extern
#endif
std::vector<std::vector<std::complex<double>>> tmp2 ;

#ifndef MESHUTILS
extern
#endif
double dx[3], Lh[3], V, gvol, k2_cutoff;

#ifndef MESHUTILS
extern
#endif
std::vector<std::vector<double>> rho, grid_W;

#ifndef MESHUTILS
extern
#endif
std::vector<int> unique_types;

#ifndef MESHUTILS
extern
#endif
std::vector<std::vector<int>> grid_inds;

#ifndef MESHUTILS
extern
#endif
bool per_frame_sq_flag;

#ifndef MESHUTILS_H
#define MESHUTILS_H
void charge_grid();
void me_unstack(int, int nn[3]);
int me_stack(int nn[3]);
double me_get_k(int, double*);
void sq_routine();
void add_segment( int );
void allocate_grid_variables();

#endif
