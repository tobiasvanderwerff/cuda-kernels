#include <sys/time.h>
#include <time.h>

#include "util.h"

int *allocIntArr(int *pointer, int existing_slots) {
  int new_total_slots = existing_slots + 1;
  int *new_pointer = (int *)realloc(pointer, new_total_slots * sizeof(int));
  return new_pointer;
}

mod_prec *allocModArr(mod_prec *pointer, int existing_slots) {
  int new_total_slots = existing_slots + 1;
  mod_prec *new_pointer = (mod_prec *)realloc(pointer, new_total_slots * sizeof(mod_prec));
  return new_pointer;
}

timestamp_t getTimeStamp() {
  struct timeval now;
  gettimeofday(&now, NULL);
  return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

void cudaSuccessOrExit(cudaError_t err) {
  if (err != cudaSuccess) {
    fprintf(stderr, "CUDA Error: %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}