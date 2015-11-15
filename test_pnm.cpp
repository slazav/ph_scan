#include <cstdio>
#include <cstdlib>

#include "pnm.h"

/* read and write pnm file */
int
main(int argc, char *argv[]){

  if (argc!=3) {
    fprintf(stderr, "Usage: %s <in> <out>\n", argv[0]);
    exit(1);
  }

  PNM pnm(argv[1]);
  pnm.save(argv[2]);
}
