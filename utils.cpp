#include <cstdio>
#include <cmath>

#include <map>
#include <queue>
#include <set>
#include "utils.h"

using namespace std;

/***************************************************************************/
/* Find conversion from IR to RGB image */
cnv_t
ir_shift(PNM &rgb, PNM &ir, int neg){

  /* build a histogram of the IR channel */
  int hist[ir.get_mcol()+1];
  memset(hist,0, sizeof(hist));
  for (int i=0; i<ir.w; i++){
    for (int j=0; j<ir.h; j++){
      hist[ir.get(0,i,j)]++;
    }
  }

  /* find 0.1% dark level */
  int sum=0, ll;
  for (ll=0;ll<sizeof(hist);ll++){
    if ( (sum+=hist[ll]) >= ir.w*ir.h/1000) break;
  }

  /* put all points below the level into a vector */
  PTS ir_dark;
  for (int i=0; i<ir.w; i++){
    for (int j=0; j<ir.h; j++){
      if (ir.get(0,i,j) < ll){ ir_dark.push_back(PT(i,j));
      ir.set(0,i,j,0);}
    }
  }

  /* */
  int rad=10;
  cnv_t ret;
  int mval = ir_dark.size()*ir.get_mcol();
  int s0m=mval;

  for (int x1=-10; x1<=10; x1++){
    for (int y1=-10; y1<=10; y1++){
//       for (int x2=-10; x2<=10; x2++){
//         for (int y2=-10; y2<=10; y2++){

          cnv_t cnv; /* (0.0) -> (x1,y1), (ir.w,ir.h)->(rgb.w+x2,rgb.h+y2) */
//          cnv.kx=(rgb.w+x2-x1)/(double)ir.w;
//          cnv.ky=(rgb.h+y2-y1)/(double)ir.h;
          cnv.dx=x1; cnv.dy=y1;
          cnv.kx=1;
          cnv.ky=1;

          /* calculate sums */
          long long s0=0;
          for (PTS::const_iterator i=ir_dark.begin(); i!=ir_dark.end(); i++){
            int xd = cnv.kx * i->x + cnv.dx;
            int yd = cnv.ky * i->y + cnv.dy;
            if (xd<0 || xd>=rgb.w || yd<0 || yd>=rgb.h) continue;
            s0 += rgb.get(0,xd,yd);
          }
          if (neg) s0 = mval - s0;
          if (s0<s0m) {s0m=s0; ret=cnv; }
//        }
//      }
    }
  }
//  fprintf(stderr, "> %f %f %f %f\n", ret.kx, ret.ky, ret.dx, ret.dy);
  return ret;
}

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

/**************************************************************************/
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
       if (mask.get(0,x,y)==1) mask.set(0,x,y,0xFFFF);
    }
  }
}

/********************/

/* find a one-color spot around a point p */
set<PT>
get_spot(const PNM &mask, const PT& p, int max){
  set<PT> ret;
  queue<PT> q;

  int h = mask.get(0,p.x,p.y);
  q.push(p); ret.insert(p);

  while (!q.empty()){
    PT p1 = q.front();
    q.pop();
    for (int i=0; i<8; i++){
      PT p2 = p1.adj(i);
      if (p2.x<0 || p2.x>=mask.w || p2.y<0 || p2.y>=mask.h) continue;
      if ((mask.get(0,p2.x,p2.y) == h)&&(ret.insert(p2).second)) q.push(p2);
    }
    if ((max!=0)&&(ret.size()>max)) break;
  }
  return ret;
}

/* find a border of a point set */
set<PT>
border(const set<PT> & pset){
  set<PT> ret;
  set<PT>::const_iterator it;
  for (it = pset.begin(); it != pset.end(); it++){
    for (int i=0; i<8; i++){
      PT p = it->adj(i);
      if (pset.count(p)==0) ret.insert(p);
    }
  }
  return ret;
}

void
interp(PNM &rgb, PNM &mask, cnv_t *cnv){
  int x,y,x0,y0,xm,ym,xp,yp;
  int rm,rp,gm,gp,bm,bp,dm,dp;
  int intx,inty;

  for (x=0; x<mask.w; x++){
    for (y=0; y<mask.h; y++){
      //  do we need interpolation of this point?
      if (mask.get(0,x,y)==0) continue;


      // find the whole interpolation area and its border
      set<PT> pset = get_spot(mask, PT(x,y), 10000);
      set<PT> bord = border(pset);
      set<PT>::iterator pi, bi;

      for (pi = pset.begin(); pi != pset.end(); pi++){
        // convert to rgb coordinates
        int xd = cnv->kx*pi->x + cnv->dx;
        int yd = cnv->ky*pi->y + cnv->dy;
        if (xd<0 || xd>=rgb.w || yd<0 || yd>=rgb.h) continue;

        double SR=0, SG=0, SB=0, S=0;
        for (bi = bord.begin(); bi != bord.end(); bi++){

          // convert to rgb coordinates
          int bxd = cnv->kx*bi->x + cnv->dx;
          int byd = cnv->ky*bi->y + cnv->dy;
          if (bxd<0 || bxd>=rgb.w || byd<0 || byd>=rgb.h) continue;

          double w = 1.0/(double)(pow(bxd - xd,2) + pow(byd - yd,2));

          S  += w;
          SR += w*rgb.get(0,bxd,byd);
          if (rgb.is_rgb()){
            SG += w*rgb.get(1,bxd,byd);
            SB += w*rgb.get(2,bxd,byd);
          }
        }

        rgb.set(0, xd, yd, SR/S);
        mask.set(0, pi->x, pi->y,0);
        if (rgb.is_rgb()){
          rgb.set(1, xd, yd, SG/S);
          rgb.set(2, xd, yd, SB/S);
        }
      }
    }
  }
}


