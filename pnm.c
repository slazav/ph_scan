#include "pnm.h"
#include <string.h>
#include <math.h>

/**************************************************************************/
/* load a file*/
pnm_t *
pnm_load(const char *fname){
  char c;
  FILE *IN;

  /* allocate the structure */
  pnm_t *pnm = (pnm_t *)malloc(sizeof(pnm_t));
  if (pnm==NULL) return NULL;

  /* open file ant read the header */
  IN=fopen(fname, "r");
  if (IN==NULL) {
    fprintf(stderr, "Can't open file: %s\n", fname);
    free(pnm); return NULL;
  }

  fread(&c,1,1,IN);  /* P */
  fscanf(IN, "%d", &(pnm->type));
  do{ fread(&c,1,1,IN); } while (c!='\n');  /* \n */

  c=fgetc(IN);
  if (c=='#'){
    do{ fread(&c,1,1,IN); } while (c!='\n');}  /* skip comment */
  else ungetc(c, IN);

  fscanf(IN, "%d %d", &(pnm->w), &(pnm->h));
  fscanf(IN, "%d", &(pnm->mcol));
  do{ fread(&c,1,1,IN); } while (c!='\n');  /* \n */

  if      (pnm->type==5 && pnm->mcol==255) pnm->bpp=1; /*  8bit grayscale */
  else if (pnm->type==6 && pnm->mcol==255) pnm->bpp=3; /* 24bit color */
  else if (pnm->type==5 && pnm->mcol>255)  pnm->bpp=2; /* 16bit grayscale */
  else if (pnm->type==6 && pnm->mcol>255)  pnm->bpp=6; /* 48bit color */
  else {
    fprintf(stderr, "Unsupported file type: %s: %d\n", fname, pnm->type);
    fclose(IN); free(pnm); return NULL;
  }

  /* allocate the buffer */
  pnm->buf=(unsigned char *)malloc(pnm->w*pnm->h*pnm->bpp);
  if (pnm->buf==NULL) {
    fprintf(stderr, "Can't allocate memory!\n"); exit(0);
    fclose(IN); free(pnm); return NULL;
  }

  /* read data */
  if (fread(pnm->buf, pnm->bpp, pnm->w*pnm->h, IN)!=pnm->w*pnm->h){
    fprintf(stderr, "Can't read the file: %s\n", fname);
    fclose(IN); free(pnm); return NULL;
  }

  /* close the file */
  fclose(IN);
  return pnm;
}

int
pnm_save(pnm_t *pnm, const char *fname){
  FILE *OUT = fopen(fname, "w");
  if (OUT==NULL){
    fprintf(stderr, "can't open file: %s\n", fname);
    return 1;
  }
  fprintf(OUT, "P%d\n%d %d\n%d\n",pnm->type,pnm->w,pnm->h,pnm->mcol);
    if (fwrite(pnm->buf, pnm->bpp, pnm->w*pnm->h, OUT)!=pnm->w*pnm->h)
      {fprintf(stderr, "Write Error!\n"); return 1;}
    fclose(OUT);
  return 0;
}

pnm_t *
pnm_create(int w, int h, int bpp){
  /* allocate the structure */
  pnm_t *pnm = (pnm_t *)malloc(sizeof(pnm_t));
  if (pnm==NULL) return NULL;
  pnm->w = w;
  pnm->h = h;
  pnm->bpp = bpp;
  pnm->type = bpp<3?5:6;  // 1,2 or 3,6
  pnm->mcol = (bpp%2==1)?255:256*256-1; // 1,3 or 2,6

  /* allocate the buffer */
  pnm->buf=(unsigned char *)malloc(pnm->w*pnm->h*pnm->bpp);
  if (pnm->buf==NULL) {
    fprintf(stderr, "Can't allocate memory!\n"); exit(0);
    free(pnm); return NULL;
  }
  memset(pnm->buf, 0, pnm->w*pnm->h*pnm->bpp);
  return pnm;
}

void
pnm_del(pnm_t *pnm){
  if (pnm->buf) free(pnm->buf);
  if (pnm) free(pnm);
}


/**************************************************************************/

/* reduce dispersion of the IR channel using RGB image */
int
ir_uncorr(pnm_t *rgb, pnm_t *ir, cnv_t *cnv){
  int x,y,xd,yd,n=0;
  /* mean values */
  double mR=0, mG=0, mB=0, mI=0;
  /* correlations */
  double RR=0, RG=0, RB=0,
               GG=0, GB=0,
                     BB=0,
         IR=0, IB=0, IG=0;
  double dR,dG,dB,dI;
  double d0,d1,d2,d3;
  double A,B,C;
  /* Weights of other channels in IR:
       The correction is I1 = I+ a R + b G + c B;
       We want to minimize <(I1 - <I1>)^2>
       d/da=0:  <(I1 - <I1>)(R - <R>)> =
                <(I -<I>)(R-<R>)> + a <(R -<R>)^2>
                   + b <(G-<G>)(R -<R>)> + b <(B-<B>)(R -<R>)> = 0

       corr(I,R) + a corr(R,R) + b corr(G,R) + c corr(B,R) = 0
       corr(I,G) + a corr(R,G) + b corr(G,G) + c corr(B,G) = 0
       corr(I,B) + a corr(R,B) + b corr(G,B) + c corr(B,B) = 0

       | RR GR BR |  a    | IR |
       | RG GG BG |  b  = | IG |
       | RB GB BB |  c    | IB |

       d0 = RR GG BB + RB GR BG + RG GB BR - RB GG BR - RG GR BB - RR GB BG;
       d1 = IR GG BB + IB GR BG + IG GB BR - IB GG BR - IG GR BB - IR GB BG;
       d2 = RR IG BB + RB IR BG + RG IB BR - RB IG BR - RG IR BB - RR IB BG;
       d3 = RR GG IB + RB GR IG + RG GB IR - RB GG IR - RG GR IB - RR GB IG;

  */

  /* calculate average values */
  for (x=0; x<ir->w; x++){
    for (y=0; y<ir->h; y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb->w || yd<0 || yd>=rgb->h) continue;
      mR += (double)pnm_getr(rgb,xd,yd); /* rgb/grey */
      mI += (double)pnm_getr(ir,x,y);
      if (rgb->type==6){ /* rgb */
        mG += (double)pnm_getg(rgb,xd,yd);
        mB += (double)pnm_getb(rgb,xd,yd);
      }
      n++;
    }
  }
  mR/=n; mG/=n; mB/=n; mI/=n;

  /* calculate correlations */
  for (x=0; x<ir->w; x++){
    for (y=0; y<ir->h; y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb->w || yd<0 || yd>=rgb->h) continue;

      dR = (double)pnm_getr(rgb,xd,yd) - mR;
      dI = (double)pnm_getr(ir,x,y) - mI;
      RR += dR*dR;
      IR += dI*dR;
      if (rgb->type==6){ /* rgb */
        dG = (double)pnm_getg(rgb,xd,yd) - mG;
        dB = (double)pnm_getb(rgb,xd,yd) - mB;
        RG += dR*dG;
        RB += dR*dB;
        GG += dG*dG;
        GB += dG*dB;
        BB += dB*dB;
        IG += dI*dG;
        IB += dI*dB;
      }
    }
  }
  /* solve linear system (for greyscale image almost all correlations are 0) */
  d0 = RR*GG*BB + RB*RG*GB + RG*GB*RB - RB*GG*RB - RG*RG*BB - RR*GB*GB;
  d1 = IR*GG*BB + IB*RG*GB + IG*GB*RB - IB*GG*RB - IG*RG*BB - IR*GB*GB;
  d2 = RR*IG*BB + RB*IR*GB + RG*IB*RB - RB*IG*RB - RG*IR*BB - RR*IB*GB;
  d3 = RR*GG*IB + RB*RG*IG + RG*GB*IR - RB*GG*IR - RG*RG*IB - RR*GB*IG;
  A = d1/d0;
  B = d2/d0;
  C = d3/d0;
  //fprintf(stderr, "> %f %f %f\n", A,B,C);

  /* modify IR image */
  for (x=0; x<ir->w; x++){
    for (y=0; y<ir->h; y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb->w || yd<0 || yd>=rgb->h) continue;
      dR = (double)pnm_getr(rgb,xd,yd)-mR; /* rgb/grey */
      dI = (double)pnm_getr(ir,x,y);
      if (rgb->type==6){ /* rgb */
        dG = (double)pnm_getg(rgb,xd,yd)-mG;
        dB = (double)pnm_getb(rgb,xd,yd)-mB;
        pnm_setr(ir, x,y, (int)(dI - A*dR - B*dG - C*dB));
      }
      else{
        pnm_setr(ir, x,y, (int)(dI - A*dR));
      }
    }
  }
  return 0;
}

pnm_t *
detect_dust1(pnm_t *ir, double thr){
  int x,y;
  double mI=0,dI=0;
  pnm_t * ret = pnm_create(ir->w, ir->h, 1);

  /* calculate average value of IR */
  for (x=0; x<ir->w; x++){
    for (y=0; y<ir->h; y++){
      mI += (double)pnm_getr(ir,x,y)/(double)(ir->w*ir->h);
    }
  }

  /* detect dust (everything below mI-thr)*/
  thr = thr * (256*256-1); /* 0..1 -> pixels */

  for (x=0; x<ir->w; x++){
    for (y=0; y<ir->h; y++){
      dI = mI - (double)pnm_getr(ir,x,y);
      if (dI > thr) pnm_setr(ret, x,y, 255);
    }
  }
  return ret;
}

void
expand_dust(pnm_t *mask){
  int x,y;
  for (x=1; x<mask->w-1; x++){
    for (y=1; y<mask->h-1; y++){
      if (pnm_getr(mask,x-1,y-1)>1) continue;
      if (pnm_getr(mask,x-1,y-1)>1 ||
          pnm_getr(mask,x-1,y  )>1 ||
          pnm_getr(mask,x-1,y+1)>1 ||
          pnm_getr(mask,x  ,y-1)>1 ||
          pnm_getr(mask,x  ,y+1)>1 ||
          pnm_getr(mask,x+1,y-1)>1 ||
          pnm_getr(mask,x+1,y  )>1 ||
          pnm_getr(mask,x+1,y+1)>1) pnm_setr(mask,x,y,1);
    }
  }
  for (x=1; x<mask->w-1; x++){
    for (y=1; y<mask->h-1; y++){
       if (pnm_getr(mask,x+1,y+1)==1) pnm_setr(mask,x,y,0xFFFF);
    }
  }
}

void
interp1(pnm_t *rgb, pnm_t *mask, cnv_t *cnv){
  int x,y,x0,y0,xm,ym,xp,yp;
  int rm,rp,gm,gp,bm,bp,dm,dp;
  int intx,inty;

  for (x=0; x<mask->w; x++){
    for (y=0; y<mask->h; y++){
      /* do we need interpolation of this point? */
      if (pnm_getr(mask,x,y)==0) continue;

      /* find a cross: nearest good points in 4 directions*/
      int m1=0,m2=0,m3=0,m4=0,i;
      for (i=x; i>0; i--)       if (pnm_getr(mask,i,y)==0) {m1=i; break;}
      for (i=x; i<mask->w; i++) if (pnm_getr(mask,i,y)==0) {m2=i; break;}
      for (i=y; i>0; i--)       if (pnm_getr(mask,x,i)==0) {m3=i; break;}
      for (i=y; i<mask->h; i++) if (pnm_getr(mask,x,i)==0) {m4=i; break;}

      /* convert coordinates to RGB image */
      x0 = cnv->kx*x + cnv->dx;
      y0 = cnv->ky*y + cnv->dy;
      xm = cnv->kx*m1 + cnv->dx - 1;
      xp = cnv->kx*m2 + cnv->dx + 1;
      ym = cnv->ky*m3 + cnv->dy - 1;
      yp = cnv->ky*m4 + cnv->dy + 1;
      if (x0<0 || x0>=rgb->w || y0<0 || y0>=rgb->h) continue;

      /* find the best interpolation direction */
      intx=inty=1;
      if (xm==x0 || xm<0 || xm>=rgb->w || xp==x0 || xp<0 || xp>=rgb->w) intx=0;
      if (ym==x0 || ym<0 || ym>=rgb->h || yp==y0 || yp<0 || yp>=rgb->h) inty=0;
      if (intx && inty){ /* choose shortest way*/
        if (xp-xm > yp-ym) intx=0; else inty=0;
      }

      /* edge colors and distances */
      if (intx){
        dm = x0-xm; dp = xp-x0;
        rm = pnm_getr(rgb,xm,y0);
        rp = pnm_getr(rgb,xp,y0);
        if (rgb->type==6){ /* rgb */
          gm = pnm_getg(rgb,xm,y0);
          gp = pnm_getg(rgb,xp,y0);
          bm = pnm_getb(rgb,xm,y0);
          bp = pnm_getb(rgb,xp,y0);
          //draw borders
          //pnm_setrgb(rgb, xm,y0, 0xFFFF,0xFFFF,0);
          //pnm_setrgb(rgb, xp,y0, 0xFFFF,0xFFFF,0);
        }
      }
      else if (inty){
        dm = y0-ym; dp = yp-y0;
        rm = pnm_getr(rgb,x0,ym);
        rp = pnm_getr(rgb,x0,yp);
        if (rgb->type==6){ /* rgb */
          gm = pnm_getg(rgb,x0,ym);
          gp = pnm_getg(rgb,x0,yp);
          bm = pnm_getb(rgb,x0,ym);
          bp = pnm_getb(rgb,x0,yp);
          //draw borders
          //pnm_setrgb(rgb, x0,ym, 0xFFFF,0,0xFFFF);
          //pnm_setrgb(rgb, x0,yp, 0xFFFF,0,0xFFFF);
        }
      }
      else continue;

      /* interpolation */
      if (rgb->type==6){
        pnm_setrgb(rgb, x0,y0,
          (rm*dp+rp*dm)/(dm+dp),
          (gm*dp+gp*dm)/(dm+dp),
          (bm*dp+bp*dm)/(dm+dp));
        //fill the dust
        //pnm_setrgb(rgb, x0,y0, 0xFFFF,0xFFFF,0);
      } else {
        pnm_setr(rgb, x0,y0,
          (rm*dp+rp*dm)/(dm+dp));
      }
    }
  }
}

/**************************************************************************/

#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
/* calculate RMS in a square */
double
rms(pnm_t *pnm, int x0, int y0, int rad) {
  int x,y,n=0;
  double mm=0, rr=0;
  for (x=MAX(0,x0-rad); x<MIN(x0+rad,pnm->w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,pnm->h); y++){
      mm += (double)pnm_getr(pnm, x,y); n++;
    }
  }
  mm/=n; n=0;
  for (x=MAX(0,x0-rad); x<MIN(x0+rad,pnm->w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,pnm->h); y++){
      rr += (double)pow(pnm_getr(pnm, x,y)-mm, 2); n++;
    }
  }
  return sqrt(rr/n);
}

/* calculate correlation between two pictures */
double
corr(pnm_t *ir, pnm_t *rgb, int x0, int y0, int rad, cnv_t *cnv) {
  int x,y,xd,yd,n=0;
  double mR=0, mI=0;
  double RR=0, II=0, IR=0;
  double dR,dI;

  for (x=MAX(0,x0-rad); x<MIN(x0+rad,ir->w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,ir->h); y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb->w || yd<0 || yd>=rgb->h) continue;
      mR += (double)pnm_getr(rgb, xd,yd);
      mI += (double)pnm_getr(ir, x,y);
      n++;
    }
  }
  if (n==0) return -1;
  mR/=n; mI/=n;
  for (x=MAX(0,x0-rad); x<MIN(x0+rad,ir->w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,ir->h); y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb->w || yd<0 || yd>=rgb->h) continue;
      dR = (double)pnm_getr(rgb, xd,yd) - mR;
      dI = (double)pnm_getr(ir, x,y) - mI;
      RR += dR*dR;
      IR += dI*dR;
      II += dI*dI;
    }
  }
  return IR/sqrt(RR*II);
}

/***************************************************************************/
/* Find shift DX,DY between Image and IR channel (todo - scaling!) */
cnv_t
ir_shift(pnm_t *rgb, pnm_t *ir, int debug){

  int rad=5;
  int maxsh=10;
  int nsteps=10;
  cnv_t cnv;

  int stepx=ir->w/nsteps, stepy=ir->h/nsteps;
  int i,j,ii,jj, im,jm, dx,dy;
  double v, vm;
  int n=0;
  double DX=0, DY=0;
  int xm[nsteps*nsteps], ym[nsteps*nsteps];

  /* Find squares with large RMS values */
  for (i=rad; i<ir->w-rad; i+=stepx){
    for (j=rad; j<ir->h-rad; j+=stepy){

      /* Find max rms in one square */
      im=i; jm=j; vm=0;
      for (ii=i+rad; ii<MIN(i+stepx,ir->w-rad); ii++){
        for (jj=j+rad; jj<MIN(j+stepy,ir->h-rad); jj++){
          v=rms(ir, ii,jj, rad);
          if ( v > vm){ vm = v; im=ii; jm=jj; }
        }
      }
      //fprintf(stderr, "rms> %f %d %d\n", vm, im, jm);

      /* compare the square with the red channel
         and find the best shift */
      dx=0; dy=0; vm=-1;
      for (ii=-maxsh; ii<=maxsh; ii++){
        for (jj=-maxsh; jj<=maxsh; jj++){
          cnv_t cnv1;
          cnv1.dx=ii; cnv1.dy=jj;
          cnv1.kx = cnv1.ky = 1;
          v = corr(ir, rgb, im, jm, rad, &cnv1);
          if ( v > vm){ vm = v; dx=ii; dy=jj; }
        }
      }

      /* debug -- print squares on the image */
      if (debug) {
        for (ii=im-rad; ii<im+rad; ii++){
          pnm_setrgb(rgb, ii,jm-rad, 0xFFFF,0xFFFF,0xFFFF);
          pnm_setrgb(rgb, ii,jm+rad, 0xFFFF,0xFFFF,0xFFFF);
          pnm_setrgb(rgb, ii+dx,jm-rad+dy, 0xFFFF,0,0);
          pnm_setrgb(rgb, ii+dx,jm+rad+dy, 0xFFFF,0,0);
        }
        for (jj=jm-rad; jj<jm+rad; jj++){
          pnm_setrgb(rgb, im-rad,jj, 0xFFFF,0xFFFF,0xFFFF);
          pnm_setrgb(rgb, im+rad,jj, 0xFFFF,0xFFFF,0xFFFF);
          pnm_setrgb(rgb, im-rad+dx,jj+dy, 0xFFFF,0,0);
          pnm_setrgb(rgb, im+rad+dx,jj+dy, 0xFFFF,0,0);
        }
      }

      DX+=dx;
      DY+=dy;
      n++;
    }
  }

  cnv.kx = 1;
  cnv.ky = 1;
  cnv.dx = DX/n;
  cnv.dy = DY/n;
  fprintf(stderr, "> %f %f %d\n", DX/n, DY/n, n);
  return cnv;
}
