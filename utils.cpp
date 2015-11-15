#include <cstdio>
#include <cmath>
#include "utils.h"

/**************************************************************************/

/* reduce dispersion of the IR channel using RGB image */
int
ir_uncorr(PNM &rgb, PNM &ir, cnv_t *cnv){
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
  for (x=0; x<ir.w; x++){
    for (y=0; y<ir.h; y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb.w || yd<0 || yd>=rgb.h) continue;
      mR += (double)rgb.get(0,xd,yd); /* rgb/grey */
      mI += (double)ir.get(0,x,y);
      if (rgb.is_rgb()){ /* rgb */
        mG += (double)rgb.get(1,xd,yd);
        mB += (double)rgb.get(2,xd,yd);
      }
      n++;
    }
  }
  mR/=n; mG/=n; mB/=n; mI/=n;

  /* calculate correlations */
  for (x=0; x<ir.w; x++){
    for (y=0; y<ir.h; y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb.w || yd<0 || yd>=rgb.h) continue;

      dR = (double)rgb.get(0,xd,yd) - mR;
      dI = (double)ir.get(0,x,y) - mI;
      RR += dR*dR;
      IR += dI*dR;
      if (rgb.is_rgb()){ /* rgb */
        dG = (double)rgb.get(0,xd,yd) - mG;
        dB = (double)rgb.get(0,xd,yd) - mB;
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
  for (x=0; x<ir.w; x++){
    for (y=0; y<ir.h; y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb.w || yd<0 || yd>=rgb.h) continue;
      dR = (double)rgb.get(0,xd,yd)-mR; /* rgb/grey */
      dI = (double)ir.get(0,x,y);
      if (rgb.is_rgb()){ /* rgb */
        dG = (double)rgb.get(1,xd,yd)-mG;
        dB = (double)rgb.get(2,xd,yd)-mB;
        ir.set(0, x,y, (int)(dI - A*dR - B*dG - C*dB));
      }
      else{
        ir.set(0, x,y, (int)(dI - A*dR));
      }
    }
  }
  return 0;
}

PNM
detect_dust1(PNM &ir, double thr){
  int x,y;
  double mI=0,dI=0;
  PNM ret(ir.w, ir.h, 1);

  /* calculate average value of IR */
  for (x=0; x<ir.w; x++){
    for (y=0; y<ir.h; y++){
      mI += (double)ir.get(0,x,y)/(double)(ir.w*ir.h);
    }
  }

  /* detect dust (everything below mI-thr)*/
  thr = thr * (256*256-1); /* 0..1 -> pixels */

  for (x=0; x<ir.w; x++){
    for (y=0; y<ir.h; y++){
      dI = mI - (double)ir.get(0,x,y);
      if (dI > thr) ret.set(0,x,y, 255);
    }
  }
  return ret;
}

void
expand_dust(PNM &mask){
  int x,y;
  for (x=1; x<mask.w-1; x++){
    for (y=1; y<mask.h-1; y++){
      if (mask.get(0,x-1,y-1)>1) continue;
      if (mask.get(0,x-1,y-1)>1 ||
          mask.get(0,x-1,y  )>1 ||
          mask.get(0,x-1,y+1)>1 ||
          mask.get(0,x  ,y-1)>1 ||
          mask.get(0,x  ,y+1)>1 ||
          mask.get(0,x+1,y-1)>1 ||
          mask.get(0,x+1,y  )>1 ||
          mask.get(0,x+1,y+1)>1) mask.set(0,x,y,1);
    }
  }
  for (x=1; x<mask.w-1; x++){
    for (y=1; y<mask.h-1; y++){
       if (mask.get(0,x+1,y+1)==1) mask.set(0,x,y,0xFFFF);
    }
  }
}

void
interp1(PNM &rgb, PNM &mask, cnv_t *cnv){
  int x,y,x0,y0,xm,ym,xp,yp;
  int rm,rp,gm,gp,bm,bp,dm,dp;
  int intx,inty;

  for (x=0; x<mask.w; x++){
    for (y=0; y<mask.h; y++){
      /* do we need interpolation of this point? */
      if (mask.get(0,x,y)==0) continue;

      /* find a cross: nearest good points in 4 directions*/
      int m1=0,m2=0,m3=0,m4=0,i;
      for (i=x; i>0; i--)      if (mask.get(0,i,y)==0) {m1=i; break;}
      for (i=x; i<mask.w; i++) if (mask.get(0,i,y)==0) {m2=i; break;}
      for (i=y; i>0; i--)      if (mask.get(0,x,i)==0) {m3=i; break;}
      for (i=y; i<mask.h; i++) if (mask.get(0,x,i)==0) {m4=i; break;}

      /* convert coordinates to RGB image */
      x0 = cnv->kx*x + cnv->dx;
      y0 = cnv->ky*y + cnv->dy;
      xm = cnv->kx*m1 + cnv->dx - 1;
      xp = cnv->kx*m2 + cnv->dx + 1;
      ym = cnv->ky*m3 + cnv->dy - 1;
      yp = cnv->ky*m4 + cnv->dy + 1;
      if (x0<0 || x0>=rgb.w || y0<0 || y0>=rgb.h) continue;

      /* find the best interpolation direction */
      intx=inty=1;
      if (xm==x0 || xm<0 || xm>=rgb.w || xp==x0 || xp<0 || xp>=rgb.w) intx=0;
      if (ym==x0 || ym<0 || ym>=rgb.h || yp==y0 || yp<0 || yp>=rgb.h) inty=0;
      if (intx && inty){ /* choose shortest way*/
        if (xp-xm > yp-ym) intx=0; else inty=0;
      }

      /* edge colors and distances */
      if (intx){
        dm = x0-xm; dp = xp-x0;
        rm = rgb.get(0,xm,y0);
        rp = rgb.get(0,xp,y0);
        if (rgb.is_rgb()){ /* rgb */
          gm = rgb.get(1,xm,y0);
          gp = rgb.get(1,xp,y0);
          bm = rgb.get(2,xm,y0);
          bp = rgb.get(2,xp,y0);
          //draw borders
          //rgb.setrgb(xm,y0, 0xFFFF,0xFFFF,0);
          //rgb.setrgb(xp,y0, 0xFFFF,0xFFFF,0);
        }
      }
      else if (inty){
        dm = y0-ym; dp = yp-y0;
        rm = rgb.get(0,x0,ym);
        rp = rgb.get(0,x0,yp);
        if (rgb.is_rgb()){ /* rgb */
          gm = rgb.get(1,x0,ym);
          gp = rgb.get(1,x0,yp);
          bm = rgb.get(2,x0,ym);
          bp = rgb.get(2,x0,yp);
          //draw borders
          //rgb.setrgb(x0,ym, 0xFFFF,0,0xFFFF);
          //rgb.setrgb(x0,yp, 0xFFFF,0,0xFFFF);
        }
      }
      else continue;

      /* interpolation */
      if (rgb.is_rgb()){
        rgb.setrgb(x0,y0,
          (rm*dp+rp*dm)/(dm+dp),
          (gm*dp+gp*dm)/(dm+dp),
          (bm*dp+bp*dm)/(dm+dp));
        //fill the dust
        //pnm_setrgb(rgb, x0,y0, 0xFFFF,0xFFFF,0);
      } else {
        rgb.set(0,x0,y0,
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
rms(PNM &pnm, int x0, int y0, int rad) {
  int x,y,n=0;
  double mm=0, rr=0;
  for (x=MAX(0,x0-rad); x<MIN(x0+rad,pnm.w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,pnm.h); y++){
      mm += (double)pnm.get(0,x,y); n++;
    }
  }
  mm/=n; n=0;
  for (x=MAX(0,x0-rad); x<MIN(x0+rad,pnm.w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,pnm.h); y++){
      rr += (double)pow(pnm.get(0,x,y)-mm, 2); n++;
    }
  }
  return sqrt(rr/n);
}

/* calculate correlation between two pictures */
double
corr(PNM &ir, PNM &rgb, int x0, int y0, int rad, cnv_t *cnv) {
  int x,y,xd,yd,n=0;
  double mR=0, mI=0;
  double RR=0, II=0, IR=0;
  double dR,dI;

  for (x=MAX(0,x0-rad); x<MIN(x0+rad,ir.w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,ir.h); y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb.w || yd<0 || yd>=rgb.h) continue;
      mR += (double)rgb.get(0,xd,yd);
      mI += (double)ir.get(0,x,y);
      n++;
    }
  }
  if (n==0) return -1;
  mR/=n; mI/=n;
  for (x=MAX(0,x0-rad); x<MIN(x0+rad,ir.w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,ir.h); y++){
      xd = cnv->kx*x + cnv->dx;
      yd = cnv->ky*y + cnv->dy;
      if (xd<0 || xd>=rgb.w || yd<0 || yd>=rgb.h) continue;
      dR = (double)rgb.get(0,xd,yd) - mR;
      dI = (double)ir.get(0,x,y) - mI;
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
ir_shift(PNM &rgb, PNM &ir, int debug){

  int rad=5;
  int maxsh=10;
  int nsteps=10;
  cnv_t cnv;

  int stepx=ir.w/nsteps, stepy=ir.h/nsteps;
  int i,j,ii,jj, im,jm, dx,dy;
  double v, vm;
  int n=0;
  double DX=0, DY=0;
  int xm[nsteps*nsteps], ym[nsteps*nsteps];

  /* Find squares with large RMS values */
  for (i=rad; i<ir.w-rad; i+=stepx){
    for (j=rad; j<ir.h-rad; j+=stepy){

      /* Find max rms in one square */
      im=i; jm=j; vm=0;
      for (ii=i+rad; ii<MIN(i+stepx,ir.w-rad); ii++){
        for (jj=j+rad; jj<MIN(j+stepy,ir.h-rad); jj++){
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
          rgb.setrgb(ii,jm-rad, 0xFFFF,0xFFFF,0xFFFF);
          rgb.setrgb(ii,jm+rad, 0xFFFF,0xFFFF,0xFFFF);
          rgb.setrgb(ii+dx,jm-rad+dy, 0xFFFF,0,0);
          rgb.setrgb(ii+dx,jm+rad+dy, 0xFFFF,0,0);
        }
        for (jj=jm-rad; jj<jm+rad; jj++){
          rgb.setrgb(im-rad,jj, 0xFFFF,0xFFFF,0xFFFF);
          rgb.setrgb(im+rad,jj, 0xFFFF,0xFFFF,0xFFFF);
          rgb.setrgb(im-rad+dx,jj+dy, 0xFFFF,0,0);
          rgb.setrgb(im+rad+dx,jj+dy, 0xFFFF,0,0);
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
