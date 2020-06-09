/*
* @Author: Aakash Gajjar
* @Date:   2017-02-15 13:33:38
* @Last Modified by:   Aakash Gajjar
* @Last Modified time: 2017-02-15 14:10:30
* @Subject:
*/

/* Program to print all permutations of a string in sorted order.*/

#include <stdio.h>
#include "chain_states.c"
#include "chain_hamiltonian.c"

int main(int argc, const char *argv[]) {
  int make_binary_file = 0;
  if (argc > 0) {
    make_binary_file = atoi(argv[1]);
  }

  basis();
  hamiltonian(make_binary_file);
  return 0;
}
