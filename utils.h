#include "pnm.h"

typedef struct{ /* 4-parameter conversion */
  double kx,ky,dx,dy;
} cnv_t;


/**************************************************************************/

/* reduce dispersion of the IR channel using RGB image */
int ir_uncorr(PNM &rgb, PNM &ir, cnv_t *cnv);

PNM detect_dust1(PNM &ir, double thr);
void expand_dust(PNM &mask);
void interp1(PNM &rgb, PNM &mask, cnv_t *cnv);

/* calculate RMS in a square */
double rms(PNM &pnm, int x0, int y0, int rad);

/* calculate correlation between two pictures */
double corr(PNM &ir, PNM &rgb, int x0, int y0, int rad, cnv_t *cnv);

/* Find shift DX,DY between Image and IR channel (todo - scaling!) */
cnv_t ir_shift(PNM &rgb, PNM &ir, int debug);
