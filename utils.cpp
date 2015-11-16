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
ir_shift(const PNM &rgb, const PNM &ir, int neg){
  PT p;

  /* build a histogram of the IR channel */
  int hist[ir.get_mcol()+1];
  memset(hist,0, sizeof(hist));
  for (p.x=0; p.x<ir.w; p.x++){
    for (p.y=0; p.y<ir.h; p.y++){
      hist[ir.get(0,p)]++;
    }
  }

  /* find 0.1% dark level */
  int sum=0, ll;
  for (ll=0;ll<sizeof(hist);ll++){
    if ( (sum+=hist[ll]) >= ir.w*ir.h/1000) break;
  }

  /* put all points below the level into a vector */
  PTS ir_dark;
  for (p.x=0; p.x<ir.w; p.x++){
    for (p.y=0; p.y<ir.h; p.y++){
      if (ir.get(0,p) < ll) ir_dark.push_back(p);
    }
  }

  /* */
  int rad=10;
  cnv_t ret;
  int mval = ir.get_mcol();
  int s0m=mval;

  for (int x1=-rad; x1<=rad; x1++){
    for (int y1=-rad; y1<=rad; y1++){
//       for (int x2=-rad; x2<=rad; x2++){
//         for (int y2=-rad; y2<=rad; y2++){

          cnv_t cnv; /* (0.0) -> (x1,y1), (ir.w,ir.h)->(rgb.w+x2,rgb.h+y2) */
          cnv.dx=x1; cnv.dy=y1;

//          cnv.kx=(rgb.w+x2-x1)/(double)ir.w;
//          cnv.ky=(rgb.h+y2-y1)/(double)ir.h;
            cnv.kx=1; cnv.ky=1;

          /* calculate sums */
          long long s0=0;
          int n=0;
          for (PTS::const_iterator i=ir_dark.begin(); i!=ir_dark.end(); i++){
            PT pd(cnv.kx * i->x + cnv.dx,
                  cnv.ky * i->y + cnv.dy);
            if (!rgb.is_in(pd)) continue;
            s0 += rgb.get(0,pd);
            n++;
          }
          s0/=n;
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
ir_uncorr(const PNM &rgb, PNM &ir, const cnv_t &cnv){
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
      PT pd(cnv.kx*x + cnv.dx,
            cnv.ky*y + cnv.dy);
      if (!rgb.is_in(pd)) continue;
      mR += (double)rgb.get(0,pd); /* rgb/grey */
      mI += (double)ir.get(0,PT(x,y));
      if (rgb.is_rgb()){ /* rgb */
        mG += (double)rgb.get(1,pd);
        mB += (double)rgb.get(2,pd);
      }
      n++;
    }
  }
  mR/=n; mG/=n; mB/=n; mI/=n;

  /* calculate correlations */
  for (x=0; x<ir.w; x++){
    for (y=0; y<ir.h; y++){
      PT pd(cnv.kx*x + cnv.dx,
            cnv.ky*y + cnv.dy);
      if (!rgb.is_in(pd)) continue;

      dR = (double)rgb.get(0,pd) - mR;
      dI = (double)ir.get(0,PT(x,y)) - mI;
      RR += dR*dR;
      IR += dI*dR;
      if (rgb.is_rgb()){ /* rgb */
        dG = (double)rgb.get(0,pd) - mG;
        dB = (double)rgb.get(0,pd) - mB;
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
      PT p(x,y);
      PT pd(cnv.kx*x + cnv.dx,
            cnv.ky*y + cnv.dy);
      if (!rgb.is_in(pd)) continue;
      dR = (double)rgb.get(0,pd)-mR; /* rgb/grey */
      dI = (double)ir.get(0,p);
      if (rgb.is_rgb()){ /* rgb */
        dG = (double)rgb.get(1,pd)-mG;
        dB = (double)rgb.get(2,pd)-mB;
        ir.set(0, p, (int)(dI - A*dR - B*dG - C*dB));
      }
      else{
        ir.set(0, p, (int)(dI - A*dR));
      }
    }
  }
  return 0;
}

PNM
detect_dust1(PNM &ir, double thr){
  PT p;
  double mI=0,dI=0;
  PNM ret(ir.w, ir.h, 1);

  /* calculate average value of IR */
  for (p.x=0; p.x<ir.w; p.x++){
    for (p.y=0; p.y<ir.h; p.y++){
      mI += (double)ir.get(0,p)/(double)(ir.w*ir.h);
    }
  }

  /* detect dust (everything below mI-thr)*/
  thr = thr * (256*256-1); /* 0..1 -> pixels */

  for (p.x=0; p.x<ir.w; p.x++){
    for (p.y=0; p.y<ir.h; p.y++){
      dI = mI - (double)ir.get(0,p);
      if (dI > thr) ret.set(0,p, 255);
    }
  }
  return ret;
}

/**************************************************************************/
void
expand_dust(PNM &mask){
  PT p;
  for (p.x=1; p.x<mask.w-1; p.x++){
    for (p.y=1; p.y<mask.h-1; p.y++){
      if (mask.get(0,p)>1) continue;
      for (int i=0;i<8;i++){
        if (mask.get(0,p.adj(0))>1) mask.set(0,p,1);
      }
    }
  }
  for (p.x=1; p.x<mask.w-1; p.x++){
    for (p.y=1; p.y<mask.h-1; p.y++){
       if (mask.get(0,p)==1) mask.set(0,p,0xFF);
    }
  }
}

/********************/

/* find a one-color spot around a point p */
set<PT>
get_spot(const PNM &mask, const PT& p, int max){
  set<PT> ret;
  queue<PT> q;

  int h = mask.get(0,p);
  q.push(p); ret.insert(p);

  while (!q.empty()){
    PT p1 = q.front();
    q.pop();
    for (int i=0; i<8; i++){
      PT p2 = p1.adj(i);
      if (!mask.is_in(p2)) continue;
      if ((mask.get(0,p2) == h)&&(ret.insert(p2).second)) q.push(p2);
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
interp(PNM &rgb, PNM &mask, const cnv_t &cnv){
  PT p;

  for (p.x=0; p.x<mask.w; p.x++){
    for (p.y=0; p.y<mask.h; p.y++){
      //  do we need interpolation of this point?
      if (mask.get(0,p)!=255) continue;


      // find the whole interpolation area and its border
      set<PT> pset = get_spot(mask, p, 1000);
      set<PT> bord = border(pset);
      set<PT>::iterator pi, bi;

      for (pi = pset.begin(); pi != pset.end(); pi++){
        // convert to rgb coordinates
        PT pd(cnv.kx*pi->x + cnv.dx,
              cnv.ky*pi->y + cnv.dy);
        if (!rgb.is_in(pd)) continue;

        double SR=0, SG=0, SB=0, S=0;
        for (bi = bord.begin(); bi != bord.end(); bi++){

          // convert to rgb coordinates
          PT bpd(cnv.kx*bi->x + cnv.dx,
                 cnv.ky*bi->y + cnv.dy);
          if (!rgb.is_in(bpd)) continue;
          double w = 1.0/(double)(pow(bpd.x-pd.x,2) + pow(bpd.y-pd.y,2));
          w = w*w;

          S  += w;
          SR += w*rgb.get(0,bpd);
          if (rgb.is_rgb()){
            SG += w*rgb.get(1,bpd);
            SB += w*rgb.get(2,bpd);
          }
        }
        mask.set(0,*pi,254);

        /* interpolate only darker points! */
        //if (rgb.get(0, pd) < SR/S) rgb.set(0, pd, SR/S);
        //if (rgb.is_rgb()){
        //  if (rgb.get(1, pd) < SG/S) rgb.set(1, pd, SG/S);
        //  if (rgb.get(2, pd) < SB/S) rgb.set(2, pd, SB/S);
        //}

        /* interpolate only darker points! */
        rgb.set(0, pd, SR/S);
        if (rgb.is_rgb()){
          rgb.set(1, pd, SG/S);
          rgb.set(2, pd, SB/S);
        }

      }
    }
  }
}


