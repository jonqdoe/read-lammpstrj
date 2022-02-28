/////////////////////////////////////////////////////////
// 10/26/2020                           R. Riggleman   //
// Uses classes to read in lammpstrj files.            //
/////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <vector>
using namespace std;

#define LAMMPSTRJ
#include "lammpstrj.h"

int parse_ATOMS_line(string);

// ASSUMPTIONS IN THIS CODE:
// - 3D positions are printed even if a 2D simulation
// - y and z immediately follow x in the column ordering
// - z is printed even if 2D simulation
// - id, type, mol are all read as ints, everything else as double
// - header is timesteps, number of atoms, box dimensions, atomic stuff
// - compute quantities are all the final columns
// - id, type, mol are constant throughout the simulation

int read_lammpstrj( const char* name, const int frame1, const int lastframe) {

  string word, line, rname ;
  
  id_col = -1;
  mol_col = -1;
  type_col = -1;
  compute_col = -1;
  n_computes = 0;

  int nread = 0, time0, time1; 
  int nstored = 0;

  int frs = lastframe - frame1 ;

  xhi.resize(frs,std::vector<double>(3));
  xlo.resize(frs,std::vector<double>(3));

  string st_name(name);
  ifstream inp(st_name);
  if (not inp.is_open()) {
    cout << "Cannot open file " + st_name << endl;
  }
  while ( !inp.eof() ) {
    // READ TIMESTEP
    getline(inp, line);
    getline(inp, line);

    istringstream iss(line);

    if ( nread == 0 ) {
      string toint ;
      iss >> toint ;
      time0 = stoi(toint);
    }
    else if ( nread == 1 ) {
      string toint ;
      iss >> toint ;
      time1 = stoi(toint);

      delta_frames = time1 - time0;
    }

    // READ NUMBER OF ATOMS
    getline(inp, line);
    getline(inp, line);
    if ( nread == 0 ) {
      nsites = stoi(line);
      xt.resize(frs, vector<vector<double>>(nsites, vector<double>(3)));
      L.resize(frs, vector<double>(3)) ;
    }

    // READ BOX DIMENSIONS
    getline(inp, line) ;

    for ( int j=0 ; j<3 ; j++ ) {
      getline(inp, line) ;
      istringstream bx(line);

      string todub ;
      bx >> todub ;
      double xlotmp = stod(todub);
      bx >> todub ;
      double xhitmp = stod(todub);

      if ( nread >= frame1 ) {
        double index = nread - frame1;
        xhi.at(nread-frame1).at(j)= xhitmp;
        xlo.at(nread-frame1).at(j)= xlotmp;
        L[nread-frame1][j] = xhi.at(nread-frame1).at(j) - xlo.at(nread-frame1).at(j) ;
      }
      
    }


    // read in ITEM: ATOMS line
    getline(inp, line);
    

    // Parse the inputs and reallocate some arrays
    if ( nread == 0 ) {
      n_col = parse_ATOMS_line(line) ;
      cout << "Found " << n_col << " columns and " << n_computes << " computes" << endl;
      if ( has_id ) {
        cout << "id in column " << id_col << endl;
        id.resize(nsites);
      }
      if ( has_mol ) {
        cout << "molID in column " << mol_col << endl;
        mol.resize(nsites);
      }
      if ( has_type ) {
        cout << "type in column " << type_col << endl;
        type.resize(nsites);
      }
      if ( n_computes > 0 ) {
        compute_col = has_id + has_mol + has_type + 3 ;
        cout << "Computes begin in column " << compute_col << endl;
        computes.resize(frs, vector<vector<double>>(nsites, vector<double>(n_computes)));
      }
    }



    ////////////////////////////////
    // Read the main data section //
    ////////////////////////////////
    for ( int i=0 ; i<nsites ; i++ ) {
      getline(inp, line) ;

      if ( nread < frame1 )
        continue ;


      istringstream is2(line);
      string toint, tofloat ;

      int cur_ind = -1;

      for ( int j=0 ; j<n_col ; j++ ) {

        
        if ( j == id_col ) {
          is2 >> toint ;
          cur_ind = stoi(toint)-1;
          if ( nread == frame1 ) id[cur_ind] = stoi(toint);
        }

        else if ( j == mol_col ) {
          is2 >> toint ;
          if ( cur_ind == -1 ) { cout << "index error!" << endl; exit(1); }
          if ( nread == frame1 ) mol[cur_ind] = stoi(toint);
        }

        else if ( j == type_col ) {
          is2 >> toint ; 
          if ( cur_ind == -1 ) { cout << "index error!" << endl; exit(1); }
          if ( nread == frame1 ) type[cur_ind] = stoi(toint);
        }

        else if ( j == x_pos_col ) {
          is2 >> tofloat ;
          if ( cur_ind == -1 ) { cout << "index error!" << endl; exit(1); }
          if ( nread >= frame1 ) xt[nread-frame1][cur_ind][0] = stod(tofloat) ;
        }

        else if ( j == x_pos_col + 1 ) {
          is2 >> tofloat ;
          if ( cur_ind == -1 ) { cout << "index error!" << endl; exit(1); }
          if ( nread >= frame1 ) xt[nread-frame1][cur_ind][1] = stod(tofloat) ;
        }

        else if ( j == x_pos_col + 2) {
          is2 >> tofloat ;
          if ( cur_ind == -1 ) { cout << "index error!" << endl; exit(1); }
          if ( nread >= frame1 ) xt[nread-frame1][cur_ind][2] = stod(tofloat) ;
        }

        else { 
          is2 >> tofloat ;
          if ( cur_ind == -1 ) { cout << "index error!" << endl; exit(1); }
          if ( nread >= frame1 ) computes[nread-frame1][cur_ind][j-n_computes] = stod(tofloat);
        }


      }// j=0:n_col

    }// i=0:nsites

    nread++;
    if ( nread > frame1 )
      nstored++;
    if ( nread == lastframe ) {
      cout << "Breaking before reading " << lastframe << ", kept " << nstored << endl;
      break ;
    }

  }// while ( !inp.eof() )
  inp.close();

  nframes = nstored ;

  return nstored ;
}




int parse_ATOMS_line( string line ) {
  istringstream isa(line);

  int cur_col = 0;

  string word2 ;
  isa >> word2;  // ITEM:
  isa >> word2;  // ATOMS

  while ( isa >> word2 ) {
    if ( word2 == "id" ) {
      has_id = 1;
      id_col = cur_col;
    }

    else if ( word2 == "mol" ) {
      has_mol = 1;
      mol_col = cur_col;
    }

    else if ( word2 == "type" ) {
      has_type = 1;
      type_col = cur_col ;
    }

    else if ( word2 == "x" || word2 == "f_avg_x" ) {
      x_pos_col = cur_col ;
      isa >> word2 ;
      isa >> word2 ;
    }

    else {
      n_computes++;
    }

    cur_col++ ;
  }//while(isa>> word2)
  return cur_col + 2;
}
