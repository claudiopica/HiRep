#include <execinfo.h>
#include <stdlib.h>
#include "logger.h"

/* Obtain a backtrace and print it to stdout. */
void print_trace (void) {
  void *array[50];
  char **strings;
  int size, i;

  size = backtrace (array, 50);
  strings = backtrace_symbols (array, size);
  if (strings != NULL)
  {

    lprintf("TRACE",1,"Obtained %d stack frames.\n", size);
    for (i = 0; i < size; i++)
      lprintf("TRACE",1,"%s\n", strings[i]);
  }

  free (strings);
}