#define LOG_SPACE
#include <cmath>
#include <iostream>
#include "log_space.h"

using namespace std ;

void log_space(int low, int high, int NP) {
  high -= 1;

  int i, j;
  double lg_low = log10(double(low));
  double lg_high = log10(double(high));
  double spacing = (lg_high - lg_low) / (double)NP;
 
  cout << lg_low << " " << lg_high << " " << spacing << endl;
  cout << low << " " << high << " " << NP << endl;

  int *temp;
  temp = (int* ) calloc(NP,sizeof(int));
  double arg;

  n_times = 1;

  int last_new = -1;

  for (i=0; i<NP; i++) {
    arg = lg_low + spacing*(double)i;
    arg = pow(10.0,arg);
    temp[i] = (int)arg;

    if (i==0)
      last_new = temp[i];
    else if (temp[i] != last_new) {
      n_times++;
      last_new = temp[i];
    }

  }
 
  times = (int*) calloc(n_times+1,sizeof(int));
  int ct=0;
  for (i=0; i<NP; i++) {
    if (i==0) {
      times[ct] = temp[i];
      ct++;
    }
    else if (times[ct-1] != temp[i]) {
      times[ct] = temp[i];
      ct++;
    }
  }
  
  times[n_times] = high-1;
  n_times++;

}



