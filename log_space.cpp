#define LOG_SPACE
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "log_space.h"

void log_space(int low, int high, int NP) {
  high -= 1;

  int i, j;
  double lg_low = log10(low);
  double lg_high = log10(high);
  double spacing = (lg_high - lg_low) / (double)NP;
  
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



