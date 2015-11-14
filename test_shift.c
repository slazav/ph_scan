#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "pnm.h"

/* find best shift between IR and RGB images */
int
main(int argc, char *argv[]){

  pnm_t *pnm, *ir;
  cnv_t cnv;

  if (argc!=4) {
    fprintf(stderr, "Usage: %s <RGB file> <IR file> <out file>\n", argv[0]);
    exit(1);
  }

  /***************************************************************************/

  pnm = pnm_load(argv[1]);
  ir  = pnm_load(argv[2]);

  cnv = ir_shift(pnm,ir, 1);

  pnm_save(pnm, argv[3]);
  pnm_del(pnm);
  pnm_del(ir);
}
