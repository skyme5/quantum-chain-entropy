/*
* @Author: Aakash Gajjar
* @Date:   2017-02-11 22:17:47
* @Last Modified by:   Aakash Gajjar
* @Last Modified time: 2017-03-09 16:47:47
* @Subject:
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define true 1
#define false 0

#define S_EIGEN 0.5
#define S_COEFF 2.0

char *make_file_name(int states, float Omega, float J_interaction, int make_binary_file) {
  char name[] = "Hamiltonian";
  char saperator[] = "_";
  char *interaction = malloc(4 * sizeof(char));
  char *omega_str = malloc(4 * sizeof(char));
  char ext[] = ".txt";
  char *bin = calloc(8, sizeof(char));
  char binary[] = "binary";
  char redable[] = "redable";
  char *states_str = malloc(sizeof(states) * 4);
  char *file_name = calloc(sizeof(name) + sizeof(interaction) + sizeof(omega_str)
                           + sizeof(ext) + sizeof(states_str) + sizeof(bin) + 4,
                           sizeof(char));

  sprintf(states_str, "%d", states);
  strcat(file_name, name);
  strcat(file_name, saperator);
  sprintf(omega_str, "%4.2f", Omega);
  strcat(file_name, omega_str);
  strcat(file_name, saperator);
  sprintf(interaction, "%4.2f", J_interaction);
  strcat(file_name, interaction);
  strcat(file_name, saperator);

  if (make_binary_file) {
    strcat(file_name, binary);
  } else {
    strcat(file_name, redable);
  }

  /*strcat(file_name, states_str);*/
  strcat(file_name, ext);
  return file_name;
}

void printSTDOUT(int states, float **hamiltonian) {
  int i, j;
  for (i = 0; i < states; i++) {
    printf("\n");
    for (j = 0; j < states; j++) {
      printf("%5.2f", hamiltonian[i][j]);
    }
  }
  printf("\n");
}

int non_zero_element(int states, float **hamiltonian) {
  int i, j;
  int count = 0;
  for (i = 0; i < states; i++) {
    for (j = 0; j < states; j++) {
      if (hamiltonian[i][j] != 0)
        count++;
    }
  }
  return count;
}

void printHamiltonianFile(int length, int states, float **hamiltonian, float Omega, float J_interaction, int make_binary_file) {
  FILE *file0, *file1;
  char *file_name = make_file_name(states, Omega, J_interaction, make_binary_file);
  char buffer[256];
  int i, j;
  if (make_binary_file) {
    file0 = fopen(file_name, "wb");
    file1 = fopen("Hamiltonian.bin", "wb");
  } else {
    file0 = fopen(file_name, "w");
    file1 = fopen("Hamiltonian.txt", "w");
  }

  sprintf(buffer, "%d\n%d\n%f %f\n%d\n", length, states, Omega, J_interaction, non_zero_element(states, hamiltonian));
  fwrite(buffer, 1, strlen(buffer), file0);
  fwrite(buffer, 1, strlen(buffer), file1);
  for (i = 0; i < states; i++) {
    printf("\rPrinted Row:%d/%d", i + 1, states);
    for (j = i; j < states; j++) {
      if (hamiltonian[i][j] != 0) {
        if (make_binary_file) {
          fwrite(&i, sizeof (i), 1, file0);
          fwrite(&j, sizeof (j), 1, file1);
          fwrite(&i, sizeof (i), 1, file0);
          fwrite(&j, sizeof (j), 1, file1);
          fwrite(&hamiltonian[i][j], sizeof (hamiltonian[i][j]), 1, file0);
          fwrite(&hamiltonian[i][j], sizeof (hamiltonian[i][j]), 1, file1);
        } else {
          sprintf(buffer, "%d %d %4.2f\n", i, j, hamiltonian[i][j]);
          fwrite(buffer, 1, strlen(buffer), file0);
          fwrite(buffer, 1, strlen(buffer), file1);
        }
      }
    }
  }
  free(file_name);
  fclose(file0);
  fclose(file1);
}

void multiply_Sz_terms(int states, int length, int **basis, float **hamiltonian, float J_interaction) {
  int i, j, k;
  float answer;
  for (i = 0; i < states; i++) {
    for (j = 0; j < length - 1; j += 2) {
      answer = 1.0;
      for (k = j; k <= j + 1; k++) {
        if (basis[i][k] == 1) {
          answer *= -S_EIGEN;
        } else {
          answer *= S_EIGEN;
        }
      }
      hamiltonian[i][i] += answer;
    }
    hamiltonian[i][i] *= (S_COEFF * J_interaction);
  }
}

int get_state_num(int length, int from[], int swapA, int swapB) {
  int i;
  int answer = 0;
  for (i = 0; i < length; i++) {
    if (i == swapA) {
      answer += from[swapB] * pow(2, length - i - 1);
    } else if (i == swapB) {
      answer += from[swapA] * pow(2, length - i - 1);
    } else {
      answer += from[i] * pow(2, length - i - 1);
    }
  }
  return answer;
}

int find_state_num(int length, int find, int * inthis) {
  int i;
  for (i = 0; i < length; i++) {
    if (inthis[i] == find)
      return i;
  }
  fprintf(stderr, "Error cannot find %d\n", find);
  exit(1);
}

int multiply_SS_terms(int states, int length, int **basis, int *basis_num, float **hamiltonian, float J_interaction) {
  int i, j;
  int current_state_num, getting_state_num;
  for (j = 0; j < length - 1; j++) {
    for (i = 0; i < states; i++) {
      if (basis[i][j]^basis[i][j + 1]) {
        current_state_num = i;
        getting_state_num = get_state_num(length, basis[i], j, j + 1);
        getting_state_num = find_state_num(states, getting_state_num, basis_num);
        hamiltonian[current_state_num][getting_state_num] += S_EIGEN * J_interaction;
        hamiltonian[getting_state_num][current_state_num] += S_EIGEN * J_interaction;
      }
    }
  }
  multiply_Sz_terms(states, length, basis, hamiltonian, J_interaction);
  return 0;
}

int Summing_Sz_terms(int states, int length, int **basis, float **hamiltonian, float Omega) {
  int i, j;
  float answer = 0.0;
  for (i = 0; i < states; i++) {
    answer = 0.0;
    for (j = 0; j < length; j++) {
      if (basis[i][j] == 0) {
        answer += -S_EIGEN * Omega;
      } else {
        answer += S_EIGEN * Omega;
      }
    }
    hamiltonian[i][i] += answer;
  }
  return 0;
}

void readBasis(int states, int length, int **basis, int *basis_num, FILE *file) {
  int i, j;
  for (i = 0; i < states; i++) {
    fscanf(file, "%d", &basis_num[i]);
    for (j = 0; j < length; j++) {
      fscanf(file, "%1d", &basis[i][j]);
    }
  }
}

void hamiltonian(int make_binary_file) {
  int i;
  int states;
  int length, minUp;
  float J_interaction;
  float Omega;
  FILE *file;

  file = fopen("states.txt", "r");
  fscanf(file, "%d %d\n%d\n%f %f\n", &length, &minUp, &states, &Omega, &J_interaction);
  /*printf("%d %d %d %f %f\n", length, minUp, states, Omega, J_interaction);*/

  float **hamiltonian = calloc(states , sizeof * hamiltonian);
  int **basis = calloc(states , sizeof * basis);
  int *basis_num = calloc(states, sizeof * basis_num);
  for (i = 0; i < states; i++) {
    hamiltonian[i] = calloc(states, sizeof * (hamiltonian[i]));
    basis[i] = calloc(length, sizeof * (basis[i]));
  }
  printf("Basis and Hamiltonian Matrix Allocated.\n");

  readBasis(states, length, basis, basis_num, file);
  fclose(file);
  printf("Basis Read from file.\n");

  /*printf("Applying ESz.\n");
  Summing_Sz_terms(states, length, basis, hamiltonian, Omega);*/
  printf("Applying SS.\n");
  multiply_SS_terms(states, length, basis, basis_num, hamiltonian, J_interaction);

  printf("Printing Hamiltonian to a file.\n");
  printHamiltonianFile(length, states, hamiltonian, Omega, J_interaction, make_binary_file);
  printSTDOUT(states, hamiltonian);
  free(hamiltonian);
  printf("\rFinished!              ");
}
