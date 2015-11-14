#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "pnm.h"

/* reduce dispersion of the IR channel using RGB image
   (without any shift) */
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

  cnv.dx=cnv.dy=0; /* trivial conversion */
  cnv.kx=cnv.ky=1;

  ir_uncorr(pnm, ir, &cnv);

  pnm_save(ir,   argv[3]);
  pnm_del(pnm);
  pnm_del(ir);
}
