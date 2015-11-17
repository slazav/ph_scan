#include "pnm.h"

/**************************************************************************/

typedef struct{ /* 4-parameter conversion */
  double kx,ky,dx,dy;
} cnv_t;

/* Find conversion from IR to RGB image */
cnv_t ir_shift(const PNM &rgb, const PNM &ir, int neg=1);

/* reduce dispersion of the IR channel using RGB image */
int ir_uncorr(const PNM &rgb, PNM &ir, const cnv_t &cnv);

/* tune RGB channel at weak variations of IR channel */
int ir_mult(PNM &rgb, const PNM &ir, const cnv_t &cnv, double thr);

PNM detect_dust1(PNM &ir, double thr);
void expand_dust(PNM &mask);
void interp(PNM &rgb, PNM &mask, const cnv_t &cnv);

/* calculate RMS in a square */
double rms(PNM &pnm, int x0, int y0, int rad);

/* calculate correlation between two pictures */
double corr(PNM &ir, PNM &rgb, int x0, int y0, int rad, cnv_t *cnv);


