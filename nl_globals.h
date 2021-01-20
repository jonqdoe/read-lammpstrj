#define min(A,B) ((A)<(B) ? (A) : (B) )
#define max(A,B) ((A)>(B) ? (A) : (B) )


#ifndef NLUTILS
extern
#endif
double nl_rc, cell_sz[3] ;

#ifndef NLUTILS
extern
#endif
int tot_cells, CELL_MAX, NEIGH_MAX, cell_allocate_flag ;


#ifndef NLUTILS
extern
#endif
vector<vector<int>> neigh_ind, cell_inds;

#ifndef NLUTILS
extern
#endif
vector<int> neigh_ct, cell_n, n_cells, part_cells ;


void make_nlist( int ,  vector<vector<double>>, vector<double>, double  )  ;
void nlist_init( void ) ;
