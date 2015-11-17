#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "pnm.h"
#include "utils.h"

void usage(void){
  fprintf(stderr, "usage: undust3 OPTIONS in_file ir_file out1 out2 out3\n");
  fprintf(stderr, "           -T # -- dust threshold 0 .. 1.0 -- default 0.06\n");
  fprintf(stderr, "           -S   -- use scaling below the threshold\n");
  exit(0);
}

int
main(int argc, char *argv[]){

  cnv_t cnv;
  int i,j;
  double thr=0.06;
  int sc=0;

  /***************************************************************************/
  /* 1. Process command-line options */
  while (1){
    i=getopt(argc, argv, "T:S");
    if (i==-1) break;
    switch(i){
      case 'T': thr=atof(optarg); break;
      case 'S': sc=1; break;
      default:  usage();
    }
  }
  if (argc-optind!=5) usage();

  /***************************************************************************/

  PNM pnm(argv[optind]);
  PNM  ir(argv[optind+1]);

  cnv = ir_shift(pnm,ir,0);
  if (sc) ir_mult(pnm, ir, cnv, thr);
  ir_uncorr(pnm, ir, cnv);


  PNM mask = detect_dust1(ir, thr);
  expand_dust(mask);
  interp(pnm, mask, cnv);


  ir.save(argv[optind+2]);
  mask.save(argv[optind+3]);
  pnm.save(argv[optind+4]);
}
