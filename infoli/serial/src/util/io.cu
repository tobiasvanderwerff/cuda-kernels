#include <stdio.h>
#include <stdlib.h>

#include "io.h"

#define NUM_STATES 16

char fpeek(FILE *stream) {
  char c = (char)fgetc(stream);
  ungetc((int)c, stream);

  return c;
}

int readFileLine(FILE *pInFile, mod_prec *iAppArray, int cellCount) {
  char c = fpeek(pInFile);
  if (c == EOF) {
    return 0;
  }

  int buffSize = cellCount * 20;  // more than enough but will do for nao
  char *buffer = (char *)malloc(buffSize * sizeof(char));

  int floatsIgnored = 0;
  int floatsNeeded = cellCount;
  int usefulElementFound = 0;

  // Get one line
  if (fgets(buffer, buffSize, pInFile)) {
    // Convert the ASCII string of one element to a double precision floating
    // point value
    char *strNumber = strtok(buffer, " ");

    // int i will stand for elements already processed, hence it starts at 0
    // even after the first strtok
    int i = 0;

    while ((strNumber != NULL) && (i < cellCount)) {
      if (i >= floatsIgnored) {
        usefulElementFound = 1;
      }

      if (i >= floatsIgnored + floatsNeeded) {
        usefulElementFound = 0;
      }

      if (usefulElementFound) {
        // atof() should change if using integers or fixed point
        iAppArray[i - floatsIgnored] = atof(strNumber);
      }
      strNumber = strtok(NULL, " ");
      i++;
    }

    free(buffer);
    if (i < cellCount) {
      // BUG: if only one element is missing but the line ends in a space, the
      // error is not detected
      printf("Error: Input line doesn't have enough elements, only %d\n", i);
      exit(EXIT_FAILURE);
    }

    return 1;  // success
  } else {
    free(buffer);
    if (!feof(pInFile)) {
      printf("Error: Reading from input file didn't finish successfully\n");
      exit(EXIT_FAILURE);
    }
    return 0;  // end of file
  }
}

void readGCalFromFile(CellState *cellPtr, int cellCount) {
  FILE *fd = fopen("gcal_file.txt", "r");
  for (int i = 0; i < cellCount; i++) {
    fscanf(fd, "%*lf ");
  }
  for (int i = 0; i < cellCount; i++) {
    fscanf(fd, "%lf ", &cellPtr[i].soma.g_CaL);
  }
  fclose(fd);
}

void printState(const CellState *cellPtr, const char *filename, int cellCount) {
  FILE *fd = fopen(filename, "w");
  mod_prec *s = (mod_prec *)malloc(NUM_STATES * sizeof(mod_prec));

  for (int i = 0; i < cellCount; i++) {
    s[0] = cellPtr[i].soma.V_soma;
    s[1] = cellPtr[i].soma.Sodium_m;
    s[2] = cellPtr[i].soma.Potassium_n;
    s[3] = cellPtr[i].soma.Potassium_x_s;
    s[4] = cellPtr[i].soma.Calcium_k;
    s[5] = cellPtr[i].soma.Calcium_l;
    s[6] = cellPtr[i].dend.V_dend;
    s[7] = cellPtr[i].dend.Calcium_r;
    s[8] = cellPtr[i].dend.Potassium_s;
    s[9] = cellPtr[i].dend.Hcurrent_q;
    s[10] = cellPtr[i].dend.Ca2Plus;
    s[11] = cellPtr[i].dend.I_CaH;
    s[12] = cellPtr[i].axon.V_axon;
    s[13] = cellPtr[i].axon.Sodium_m_a;
    s[14] = cellPtr[i].axon.Sodium_h_a;
    s[15] = cellPtr[i].axon.Potassium_x_a;

    for (int j = 0; j < NUM_STATES; j++) {
      fprintf(fd, "%.8lf ", s[j]);
    }
    fprintf(fd, "\n");
  }

  fclose(fd);
}
