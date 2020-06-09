/*
* @Author: Aakash Gajjar
* @Date:   2017-02-16 15:09:50
* @Last Modified by:   Aakash Gajjar
* @Last Modified time: 2017-03-24 00:16:04
* @Subject:
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_log.h>

#define DEBUG 1
#define HAMILTONIAN 1
#define DENSITYMATRIX 1
#define RED_DENSITYMATRIX 1
#define RED_DENSITYMATRIX_EIGEN 1

void init_hamiltonian(int states, gsl_matrix *Hamiltonian) {
  int i, j;
  for (i = 0; i < states; i++)
    for (j = 0; j < states; j++)
      gsl_matrix_set(Hamiltonian, i, j, 0.0);
}

void print_matrix_to_stdout(int states, gsl_matrix *Hamiltonian, const char *msg) {
  int i, j;
  printf("\n====================================\n");
  printf("%s:\n", msg);
  for (i = 0; i < states; i++)
    for (printf("\n"), j = 0; j < states; j++)
      printf("%10.2g ", gsl_matrix_get(Hamiltonian, i, j));
  printf("\n====================================\n");
}

int main() {
  int length;
  int states;
  float Omega, J_interaction;
  int to_read;
  int i, j;
  double value;
  //double *eigen_vector;
  //clock_t start, end;
  FILE *file;

  file = fopen("Hamiltonian.txt", "r");
  if (ferror(file)) {
    fprintf(stderr, "Error opening file Hamiltonian.txt\n");
    exit(1);
  }

  fscanf(file, "%d\n%d\n%f %f\n%d",
         &length, &states,
         &Omega, &J_interaction,
         &to_read);

  gsl_matrix *Hamiltonian = gsl_matrix_alloc(states, states);
  gsl_matrix *density_matrix = gsl_matrix_alloc(states, states);
  gsl_matrix *reduced_density_matrix = gsl_matrix_alloc(states / 2, states / 2);

  init_hamiltonian(states, Hamiltonian);
  while (fscanf(file, "%d %d %lf", &i, &j, &value) == 3) {
    gsl_matrix_set(Hamiltonian, i, j, value);
    gsl_matrix_set(Hamiltonian, j, i, value);
  }
  fclose(file);

#if DEBUG == 1 && HAMILTONIAN == 1
  print_matrix_to_stdout(states, Hamiltonian, "Hamiltonian");
#endif

  gsl_vector *EigenValue_ham = gsl_vector_alloc (states);
  gsl_matrix *EigenVector_ham = gsl_matrix_alloc (states, states);
  gsl_vector *EigenValue_density = gsl_vector_alloc (states);
  gsl_matrix *EigenVector_density = gsl_matrix_alloc (states, states);
  gsl_vector *EigenValue_reduced_density = gsl_vector_alloc (states / 2);
  gsl_matrix *EigenVector_reduced_density = gsl_matrix_alloc (states / 2, states / 2);

  gsl_eigen_symmv_workspace * workspace_ham = gsl_eigen_symmv_alloc (states);
  gsl_eigen_symmv_workspace * workspace_density = gsl_eigen_symmv_alloc (states);
  gsl_eigen_symmv_workspace * workspace_reduced_density = gsl_eigen_symmv_alloc (states / 2);

  //eigen_vector = malloc(states * sizeof(eigen_vector));

  gsl_eigen_symmv (Hamiltonian, EigenValue_ham, EigenVector_ham, workspace_ham);
  //gsl_eigen_symmv_sort (EigenValue_ham, EigenVector_ham, GSL_EIGEN_SORT_ABS_ASC);
  {
    int eigen_i, eigen_j;
    for (i = 0; i < states; i++) {		// Eigenvalues of each state
      for (j = 0; j < states; j++) {	// Get Eigen Values
        printf("H:%f\n", gsl_matrix_get(EigenVector_ham, i, j));
      }
      // Construct Density MAtrix
      for (eigen_i = 0; eigen_i < states; eigen_i++) {
        //printf("%f\n", gsl_matrix_get(EigenVector_ham, i, eigen_i));
        for (eigen_j = eigen_i; eigen_j < states; eigen_j++) {
          // value = eigen_vector[eigen_i] * eigen_vector[eigen_j];
          // gsl_matrix_set(density_matrix, eigen_i, eigen_j, value);
          gsl_matrix_set(density_matrix, eigen_i, eigen_j,
                         gsl_matrix_get(EigenVector_ham, i, eigen_i) *
                         gsl_matrix_get(EigenVector_ham, i, eigen_j));
          gsl_matrix_set(density_matrix, eigen_j, eigen_i,
                         gsl_matrix_get(EigenVector_ham, i, eigen_i) *
                         gsl_matrix_get(EigenVector_ham, i, eigen_j));
        }
      }

#if DEBUG == 1 && DENSITYMATRIX == 1
      value = 0.0;
      for (eigen_i = 0; eigen_i < states; eigen_i++) {
        value += gsl_matrix_get(density_matrix, eigen_i, eigen_i);
      }
      print_matrix_to_stdout(states, density_matrix, "Density Matrix");
      printf("\nTrace: %5.2f", value);
      printf("\n====================================\n");
#endif
      for (eigen_i = 0; eigen_i < states / 2; eigen_i++) {
        for (eigen_j = 0; eigen_j < states / 2; eigen_j++) {
          gsl_matrix_set(reduced_density_matrix, eigen_i, eigen_j,
                         gsl_matrix_get(density_matrix, 2 * eigen_i, 2 * eigen_j));
        }
      }
#if DEBUG == 1 && RED_DENSITYMATRIX == 1
      print_matrix_to_stdout(states / 2, reduced_density_matrix, "Reduced Density Matrix");
#endif
      gsl_eigen_symmv(reduced_density_matrix, EigenValue_reduced_density,
                      EigenVector_reduced_density, workspace_reduced_density);
      gsl_eigen_symmv_sort (EigenValue_reduced_density,
                            EigenVector_reduced_density,
                            GSL_EIGEN_SORT_ABS_ASC); {
        value = 0;
        for (eigen_i = 0; eigen_i < states / 2; eigen_i++) {
          //value += gsl_pow_2(gsl_vector_get(EigenValue_reduced_density, eigen_i));
          printf("\nEN: %.18e", gsl_vector_get(EigenValue_reduced_density, eigen_i));
          // value += gsl_vector_get(EigenValue_reduced_density, eigen_i) *
          //          gsl_sf_log(gsl_vector_get(EigenValue_reduced_density, eigen_i));
        }
        //printf("Entropy:%f\n", value);
        printf("\n");
      }
    }
  }

  return 0;
}