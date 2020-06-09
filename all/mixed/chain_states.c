/*
*	This program requires user input of
*		Length of Chain(Total number of prticles)
*		Number of UP Particles
*		Omega
*		The Interation Strength (same in all direction)
*
*	The Output file format will be
*		Length 	UpParticles
*		TotalStates
*		Omega 	InterationStrength
*		States on each new line
*
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#define bool int
#define true 1
#define false 0

int get_states_state_num(int length, char from[]) {
  int i;
  int answer = 0;
  for (i = 0; i < length; i++) {
    answer += (from[i] - 48) * pow(2, length - i - 1);
  }
  return answer;
}


unsigned long long int nCr(int n, int r) {
  if (r > n / 2) {
    r = n - r;
  }
  unsigned long long int ans = 1;
  int i;

  for (i = 1; i <= r; i++) {
    ans *= n - r + i;
    ans /= i;
  }
  return ans;
}

/* Following function is needed for library function qsort().
*/
int compare(const void *a, const void * b) {
  return ( *(char *)a - * (char *)b );
}

/* A utility function two swap two characters a and b
*/
void swap(char* a, char* b) {
  char t = *a;
  *a = *b;
  *b = t;
}

/* This function finds the index of the smallest character
* which is greater than 'first' and is present in str[l..h]
*/
int findCeil(char str[], char first, int l, int h) {
  /*  initialize index of ceiling element*/
  int ceilIndex = l;
  int i;

  /* Now iterate through rest of the elements and find
   * the smallest character greater than 'first'
   */
  for (i = l + 1; i <= h; i++)
    if (str[i] > first && str[i] < str[ceilIndex])
      ceilIndex = i;

  return ceilIndex;
}

/* Print all permutations of str in sorted order*/
void sortedPermutations(char str[], FILE *file, FILE *same) {
  /* Get size of string*/
  int size = strlen(str);

  /* Sort the string in increasing order*/
  qsort(str, size, sizeof( str[0] ), compare);

  /* Print permutations one by one*/
  bool isFinished = false;
  while (!isFinished) {
    /* print this permutation*/
    fprintf(file, "%d %s\n", get_states_state_num(size, str), str);
    fprintf(same, "%d %s\n", get_states_state_num(size, str), str);

    /* Find the rightmost character which is smaller than its next
     * character. Let us call it 'first char'
     */
    int i;
    for (i = size - 2; i >= 0; --i) {
      if (str[i] < str[i + 1]) {
        break;
      }
    }

    /* If there is no such chracter, all are sorted in decreasing order,
     * means we just printed the last permutation and we are done.
     */
    if (i == -1) {
      isFinished = true;
    } else {
      /* Find the ceil of 'first char' in right of first character.
       * Ceil of a character is the smallest character greater than it
       */
      int ceilIndex = findCeil(str, str[i], i + 1, size - 1);

      /* Swap first and second characters*/
      swap(&str[i], &str[ceilIndex]);

      /* Sort the string on right of 'first char'*/
      qsort(str + i + 1, size - i - 1, sizeof(str[0]), compare);
    }
  }
}

void generate_states(int length, int minUp, float Omega, float J_interaction) {
  int i;
  char *str = malloc(sizeof(*str) * length + 1);
  char *name = calloc(20, sizeof(*name));
  char *str_length = calloc(10, sizeof(*str_length));
  char *str_minUp = calloc(10, sizeof(*str_minUp));
  char ch[] = "_";
  char ext[] = ".txt";
  char base[] = "states_";

  sprintf(str_length, "%d", length);
  sprintf(str_minUp, "%d", minUp);

  strcat(name, base);
  strcat(name, str_length);
  strcat(name, ch);
  strcat(name, str_minUp);
  strcat(name, ext);

  /*Assign Data to str
  *	0000....1111
  */
  for (i = 0; i < length; i++) {
    if (i < minUp)
      str[i] = '1';
    else
      str[i] = '0';
  }
  /*Null terminator*/
  str[i] = '\0';

  FILE *file, *same;
  file = fopen(name, "w");
  same = fopen("states.txt", "w");
  if (file == NULL) {
    printf("Error opening file\n");
  }

  /* We should print L & minUp
  * with total number of basis state
  */
  fprintf(file, "%d %d\n%llu\n%f %f\n", length, minUp, nCr(length, minUp), Omega, J_interaction);
  fprintf(same, "%d %d\n%llu\n%f %f\n", length, minUp, nCr(length, minUp), Omega, J_interaction);

  /* call function for permutations.
   * two arguments string: str
   *				   FILE: file
   */
  sortedPermutations(str, file, same);

  fclose(file);
  fclose(same);
}

void basis() {
  int length;
  int minUp;
  float Omega;
  float J_interaction;
  printf("Enter length of chain: ");
  scanf("%d", &length);
  printf("Enter number of spin UP particles: ");
  scanf("%d", &minUp);
  printf("Enter Omega: ");
  scanf("%f", &Omega);
  printf("Enter interaction strength between particles: ");
  scanf("%f", &J_interaction);
  if (minUp > length)
    fprintf(stderr, "Input Error: minUp should be less than length of chain\n");
  generate_states(length, minUp, Omega, J_interaction);
  printf("Configurations of chain are generated.\n");
}