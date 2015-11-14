#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "pnm.h"

int
main(int argc, char *argv[]){

  pnm_t *pnm, *ir, *mask;
  cnv_t cnv;
  int i,j;
  double thr=0.06;

  /***************************************************************************/
  /* 1. Process command-line options */
  void usage(void){
    fprintf(stderr, "usage: %s OPTIONS in_file ir_file out1 out2 out3\n", argv[0]);
    fprintf(stderr, "           -T # -- dust threshold 0 .. 1.0 -- default 0.06\n");
    exit(0);
  }
  while (1){
    i=getopt(argc, argv, "T:");
    if (i==-1) break;
    switch(i){
      case 'T': thr=atof(optarg); break;
      default:  usage();
    }
  }
  if (argc-optind!=5) usage();

  /***************************************************************************/

  pnm = pnm_load(argv[optind]);
  ir  = pnm_load(argv[optind+1]);

  cnv.dx=cnv.dy=0;
  cnv.kx=cnv.ky=1;
  cnv = ir_shift(pnm,ir,0);

  ir_uncorr(pnm, ir, &cnv);
  mask = detect_dust1(ir, thr);
  expand_dust(mask);
  interp1(pnm, mask, &cnv);

  pnm_save(ir,   argv[optind+2]);
  pnm_save(mask, argv[optind+3]);
  pnm_save(pnm,  argv[optind+4]);
  pnm_del(pnm);
  pnm_del(ir);
  pnm_del(mask);
}
