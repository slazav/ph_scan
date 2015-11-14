#ifndef PNM_H
#define PNM_H

#include <stdio.h>
#include <stdlib.h>

/* a simple library for working with pnm files */

/********************************************************************/
typedef struct{
  unsigned char *buf;
  int w, h,
      mcol, // max color
      bpp,  // bytes per point
      type; // pnm type (5 or 6)
} pnm_t;

typedef struct{ /* 4-parameter conversion */
  double kx,ky,dx,dy;
} cnv_t;

/********************************************************************/
/* load/save/create/delete pnm (see pnm.c) */
pnm_t * pnm_load(const char *fname);
int     pnm_save(pnm_t *pnm, const char *fname);
pnm_t * pnm_create(int w, int h, int bpp); // bpp-bytes per point (1,2,3,6)
void    pnm_del(pnm_t *pnm);

/********************************************************************/
/* reduce dispersion of the IR channel using RGB image */
int ir_uncorr(pnm_t *rgb, pnm_t *ir, cnv_t *cnv);

/* detect dust (everything below mI-thr)*/
pnm_t * detect_dust1(pnm_t *ir, double thr);
void    expand_dust(pnm_t *dst);

/* interpolate dust */
void interp1(pnm_t *rgb, pnm_t *mask, cnv_t *cnv);

/* Find shift DX,DY between Image and IR channel (todo - scaling!) */
cnv_t ir_shift(pnm_t *rgb, pnm_t *ir, int debug);

/********************************************************************/

/* read functions for 8 and 16 bit data (no internal checks). */
static __inline__ int
pnm_get8r(pnm_t *pnm, int x, int y){ return pnm->buf[pnm->bpp*(y*pnm->w+x)+0]; }
static __inline__ int
pnm_get8g(pnm_t *pnm, int x, int y){ return pnm->buf[pnm->bpp*(y*pnm->w+x)+1]; }
static __inline__ int
pnm_get8b(pnm_t *pnm, int x, int y){ return pnm->buf[pnm->bpp*(y*pnm->w+x)+2]; }

static __inline__ int
pnm_get16r(pnm_t *pnm, int x, int y){
  return pnm->buf[pnm->bpp*(y*pnm->w+x)+0]*256 + pnm->buf[pnm->bpp*(y*pnm->w+x)+1]; }
static __inline__ int
pnm_get16g(pnm_t *pnm, int x, int y){
  return pnm->buf[pnm->bpp*(y*pnm->w+x)+2]*256 + pnm->buf[pnm->bpp*(y*pnm->w+x)+3]; }
static __inline__ int
pnm_get16b(pnm_t *pnm, int x, int y){
  return pnm->buf[pnm->bpp*(y*pnm->w+x)+4]*256 + pnm->buf[pnm->bpp*(y*pnm->w+x)+5]; }

static __inline__ int
pnm_getr(pnm_t *pnm, int x, int y){
  if (pnm->bpp==6 || pnm->bpp==2) return pnm_get16r(pnm,x,y);
  if (pnm->bpp==3 || pnm->bpp==1) return pnm_get8r(pnm,x,y);
  return -1;
}
static __inline__ int
pnm_getg(pnm_t *pnm, int x, int y){
  if (pnm->bpp==6) return pnm_get16g(pnm,x,y);
  if (pnm->bpp==3) return pnm_get8g(pnm,x,y);
  return -1;
}
static __inline__ int
pnm_getb(pnm_t *pnm, int x, int y){
  if (pnm->bpp==6) return pnm_get16b(pnm,x,y);
  if (pnm->bpp==3) return pnm_get8b(pnm,x,y);
  return -1;
}

/* same, but with all needed checks */
#define DEF_SAFE_GET(NAME,BPP) static __inline__ int\
  pnm_get ## NAME ## _safe(pnm_t *pnm, int x, int y){\
    if (x<0 || y<0 || x>=pnm->w || y>=pnm->h || pnm->bpp!=BPP) return -1;\
    return pnm_get ## NAME (pnm, x,y);\
  }
DEF_SAFE_GET(8r, 3)
DEF_SAFE_GET(8g, 3)
DEF_SAFE_GET(8b, 3)
DEF_SAFE_GET(16r, 6)
DEF_SAFE_GET(16g, 6)
DEF_SAFE_GET(16b, 6)


/* write functions for 8 and 16 bit data (no internal checks). */
static __inline__ void
pnm_set8r(pnm_t *pnm, int x, int y, int v){
  pnm->buf[pnm->bpp*(y*pnm->w+x)+0] = v&0xFF; }
static __inline__ void
pnm_set8g(pnm_t *pnm, int x, int y, int v){
  pnm->buf[pnm->bpp*(y*pnm->w+x)+1] = v&0xFF; }
static __inline__ void
pnm_set8b(pnm_t *pnm, int x, int y, int v){
  pnm->buf[pnm->bpp*(y*pnm->w+x)+2] = v&0xFF; }
static __inline__ void
pnm_set8rgb(pnm_t *pnm, int x, int y, int r, int g, int b){
  pnm_set8r(pnm, x, y, r);
  pnm_set8g(pnm, x, y, g);
  pnm_set8b(pnm, x, y, b);
}

static __inline__ void
pnm_set16r(pnm_t *pnm, int x, int y, int v){
  pnm->buf[pnm->bpp*(y*pnm->w+x)+0] = (v>>8)&0xFF;
  pnm->buf[pnm->bpp*(y*pnm->w+x)+1] = v&0xFF; }
static __inline__ void
pnm_set16g(pnm_t *pnm, int x, int y, int v){
  pnm->buf[pnm->bpp*(y*pnm->w+x)+2] = (v>>8)&0xFF;
  pnm->buf[pnm->bpp*(y*pnm->w+x)+3] = v&0xFF; }
static __inline__ void
pnm_set16b(pnm_t *pnm, int x, int y, int v){
  pnm->buf[pnm->bpp*(y*pnm->w+x)+4] = (v>>8)&0xFF;
  pnm->buf[pnm->bpp*(y*pnm->w+x)+5] = v&0xFF; }
static __inline__ void
pnm_set16rgb(pnm_t *pnm, int x, int y, int r, int g, int b){
  pnm_set16r(pnm, x, y, r);
  pnm_set16g(pnm, x, y, g);
  pnm_set16b(pnm, x, y, b);
}

static __inline__ void
pnm_setr(pnm_t *pnm, int x, int y, int v){
  if (pnm->bpp==6 || pnm->bpp==2) pnm_set16r(pnm,x,y,v);
  if (pnm->bpp==3 || pnm->bpp==1) pnm_set8r(pnm,x,y,v);
}
static __inline__ void
pnm_setg(pnm_t *pnm, int x, int y, int v){
  if (pnm->bpp==6) return pnm_set16g(pnm,x,y,v);
  if (pnm->bpp==3) return pnm_set8g(pnm,x,y,v);
}
static __inline__ void
pnm_setb(pnm_t *pnm, int x, int y, int v){
  if (pnm->bpp==6) pnm_set16b(pnm,x,y,v);
  if (pnm->bpp==3) pnm_set8b(pnm,x,y,v);
}
static __inline__ void
pnm_setrgb(pnm_t *pnm, int x, int y, int r, int g, int b){
  if (pnm->bpp==6) pnm_set16rgb(pnm,x,y,r,g,b);
  if (pnm->bpp==3) pnm_set8rgb(pnm,x,y,r,g,b);
}


/* same, but with all needed checks */
#define DEF_SAFE_SET(NAME,BPP) static __inline__ void\
  pnm_set ## NAME ## _safe(pnm_t *pnm, int x, int y, int v){\
    if (x<0 || y<0 || x>=pnm->w || y>=pnm->h || pnm->bpp!=BPP) return;\
    pnm_set ## NAME (pnm, x,y,v);\
  }
#define DEF_SAFE_SETRGB(NAME,BPP) static __inline__ void\
  pnm_set ## NAME ## rgb_safe(pnm_t *pnm, int x, int y, int r, int g, int b){\
    if (x<0 || y<0 || x>=pnm->w || y>=pnm->h || pnm->bpp!=BPP) return;\
    pnm_set ## NAME (pnm, x,y,r,g,b);\
  }
DEF_SAFE_SET(8r, 3)
DEF_SAFE_SET(8g, 3)
DEF_SAFE_SET(8b, 3)
DEF_SAFE_SETRGB(8rgb, 3)
DEF_SAFE_SET(16r, 6)
DEF_SAFE_SET(16g, 6)
DEF_SAFE_SET(16b, 6)
DEF_SAFE_SETRGB(16rgb, 6)

#endif
