/*
* @Author: Aakash Gajjar
* @Date:   2017-01-26 23:40:31
* @Last Modified by:   Aakash Gajjar
* @Last Modified time: 2017-03-09 16:53:23
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void fast_dec_bin(int length, unsigned long x, int * c) {
  int i;
  for (i = 0; i < length; i++)
    *(c++) = (x >> i) & 0x1;
}

int fillConfigs(int rows, int cols, FILE *file) {
  unsigned long i;
  int j;
  int *state = calloc(cols , sizeof(state));

  for (i = 0; i < rows; i++) {
    fast_dec_bin(cols, i, state);
    for (j = cols - 1; j >= 0; j--) {
      fprintf(file, "%d", state[j]);
    }
    fprintf(file, "\n");
  }
  return 0;
}

int main() {
  int total_particles;
  int total_states;
  FILE *file;
  clock_t start, end;

  printf("Enter number of particles:");
  scanf("%d", &total_particles);
  printf("Generating states for: %d particles ", total_particles);

  start = clock();

  total_states = pow(2 , total_particles);
  printf("with total %d states\n", total_states);

  file = fopen("states.txt", "w+");
  if (ferror(file)) {
    printf("Error opening file\n");
    return 1;
  }

  fprintf(file, "%d %d\n", total_particles, total_states);
  fillConfigs(total_states, total_particles, file);
  close(file);

  end = clock();
  printf("Which are saved in 'states.txt', it took %f secs.\n", (double)(end - start) / CLOCKS_PER_SEC);

  return 0;
}