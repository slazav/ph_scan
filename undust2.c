#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define LIMIT(x,min,max) MIN(MAX(x, min), max)

/* some global values and useful macro to work with data arrays */
int w,h;
int cw   = 6; // bytes per point
int cwir = 2; // bytes per point in IR channel
unsigned char *buf, *irbuf;

#define GET_R(x,y)   (buf[cw*((y)*w+(x))+0]*256 + buf[cw*((y)*w+(x))+1])
#define GET_G(x,y)   (buf[cw*((y)*w+(x))+2]*256 + buf[cw*((y)*w+(x))+3])
#define GET_B(x,y)   (buf[cw*((y)*w+(x))+4]*256 + buf[cw*((y)*w+(x))+5])
#define GET_BW(x,y)  (buf[cw*((y)*w+(x))+0]*256 + buf[cw*((y)*w+(x))+1])

#define SET_RGB(x,y,r,g,b) (\
  buf[cw*((y)*w+(x))+0] = ((r)>>8)&0xFF, buf[cw*((y)*w+(x))+1]=((r)&0xFF),\
  buf[cw*((y)*w+(x))+2] = ((g)>>8)&0xFF, buf[cw*((y)*w+(x))+3]=((g)&0xFF),\
  buf[cw*((y)*w+(x))+4] = ((b)>>8)&0xFF, buf[cw*((y)*w+(x))+5]=((b)&0xFF) )
#define SET_BW(x,y,v) (\
  buf[cw*((y)*w+(x))+0] = ((v)>>8)&0xFF,\
  buf[cw*((y)*w+(x))+1] = ((v)&0xFF) )

#define GET_IR(x,y)   (irbuf[cwir*((y)*w+(x))]*256 + irbuf[cwir*((y)*w+(x))+1])
#define SET_IR(x,y,v) (\
  irbuf[2*((y)*w+(x))+0] = ((v)>>8)&0xFF,\
  irbuf[2*((y)*w+(x))+1] = ((v)&0xFF) )


/* calculate RMS of a square in the IR channel */
double
rms_ir(int x0, int y0, int rad) {
  int x,y,n=0;
  double mm=0, rr=0;
  for (x=MAX(0,x0-rad); x<MIN(x0+rad,w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,h); y++){
      mm += (double)GET_IR(x,y); n++;
    }
  }
  mm/=n; n=0;
  for (x=MAX(0,x0-rad); x<MIN(x0+rad,w); x++){
    for (y=MAX(0,y0-rad); y<MIN(y0+rad,h); y++){
      rr += (double)pow(GET_IR(x,y)- mm, 2); n++;
    }
  }
  return sqrt(rr/n);
}

/* calculate correlation between IR and Red channels */
double
corr(int x0, int y0, double dx, double dy, int rad) {
  int x,y,xd,yd,n=0;
  double mR=0, mG=0, mB=0, mI=0;
  double RR=0, II=0, IR=0;
  double dR,dI;

  for (x=x0; x<MIN(x0+rad,w); x++){
    for (y=y0; y<MIN(y0+rad,h); y++){
      xd=x+dx; yd=y+dy;
      if (xd<0 || xd>=w || yd<0 || yd>=h) continue;
      mR += (double)GET_R(xd,yd);
      mI += (double)GET_IR(x,y);
      n++;
    }
  }
  if (n==0) return -1;
  mR/=n; mI/=n;
  for (x=x0; x<x0+rad && x<w; x++){
    for (y=y0; y<y0+rad && y<h; y++){
      xd=x+dx; yd=y+dy;
      if (xd<0 || xd>=w || yd<0 || yd>=h) continue;
      dR = GET_R(xd,yd) - mR;
      dI = GET_IR(x,y) - mI;
      RR += (double)dR*dR;
      IR += (double)dI*dR;
      II += (double)dI*dI;
    }
  }
  return IR/sqrt(RR*II);
}

/* calculate correlation between IR and Red channels */
double
corr2(int x0, int y0, double dx, double dy, int rad) {
  int x,y,xd,yd,n=0;
  double mR=0, mG=0, mB=0, mI=0;
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

  for (x=x0; x<MIN(x0+rad,w); x++){
    for (y=y0; y<MIN(y0+rad,h); y++){
      xd=x+dx; yd=y+dy;
      if (xd<0 || xd>=w || yd<0 || yd>=h) continue;
      mR += (double)GET_R(xd,yd);
      mG += (double)GET_G(xd,yd);
      mB += (double)GET_B(xd,yd);
      mI += (double)GET_IR(x,y);
      n++;
    }
  }
  if (n==0) return -1;
  mR/=n; mG/=n; mB/=n; mI/=n;
  for (x=x0; x<x0+rad && x<w; x++){
    for (y=y0; y<y0+rad && y<h; y++){
      xd=x+dx; yd=y+dy;
      if (xd<0 || xd>=w || yd<0 || yd>=h) continue;
      dR = GET_R(xd,yd) - mR;
      dG = GET_G(xd,yd) - mG;
      dB = GET_B(xd,yd) - mB;
      dI = GET_IR(x,y) - mI;
      RR += (double)dR*dR;
      RG += (double)dR*dG;
      RB += (double)dR*dB;
      GG += (double)dG*dG;
      GB += (double)dG*dB;
      BB += (double)dB*dB;
      IR += (double)dI*dR;
      IG += (double)dI*dG;
      IB += (double)dI*dB;
    }
  }
  d0 = RR*GG*BB + RB*RG*GB + RG*GB*RB - RB*GG*RB - RG*RG*BB - RR*GB*GB;
  d1 = IR*GG*BB + IB*RG*GB + IG*GB*RB - IB*GG*RB - IG*RG*BB - IR*GB*GB;
  d2 = RR*IG*BB + RB*IR*GB + RG*IB*RB - RB*IG*RB - RG*IR*BB - RR*IB*GB;
  d3 = RR*GG*IB + RB*RG*IG + RG*GB*IR - RB*GG*IR - RG*RG*IB - RR*GB*IG;
  A = d1/d0;
  B = d2/d0;
  C = d3/d0;

  for (x=x0; x<x0+rad && x<w; x++){
    for (y=y0; y<y0+rad && y<h; y++){
      xd=x+dx; yd=y+dy;
      if (xd<0 || xd>=w || yd<0 || yd>=h) continue;
      dR = GET_R(xd,yd);
      dG = GET_G(xd,yd);
      dB = GET_B(xd,yd);
      dI = GET_IR(x,y);
      SET_IR(x,y, (int)(dI - A*dR - B*dG - C*dB));
    }
  }

  fprintf(stderr, "> %f %f %f\n", A,B,C);
  return 0;
}



int
main(int argc, char *argv[]){

  unsigned char *stat;
  unsigned char c;
  int      i,j;
  double DX=0, DY=0;
  int    irrad=5;
  int    irthr=round(0.055*256*256);
  const char *outd = NULL;

  /***************************************************************************/
  /* 1. Process command-line options */
  void usage(void){
    fprintf(stderr, "usage: %s OPTIONS in_file ir_file > out_file\n", argv[0]);
    fprintf(stderr, "           -T # -- infrared threshold 0 .. 1.0 -- default 0.06\n");
    fprintf(stderr, "           -R # -- interpolation radius, px    -- default 30\n");
    fprintf(stderr, "           -D name -- write dust picture in the file\n");
    exit(0);
  }
  while (1){
    i=getopt(argc, argv, "T:R:D:");
    if (i==-1) break;
    switch(i){
      case 'T': irthr=round(256*256*atof(optarg)); break;
      case 'R': irrad=atoi(optarg); break;
      case 'D': outd=optarg; break;
      default:  usage();
    }
  }
  if (argc-optind!=2) usage();

  /***************************************************************************/
  /* 2. Read files */
  {
    int      type, irtype;
    unsigned mp, irw,irh,irmp;
    FILE     *IN=NULL, *IR=NULL;

    /* Open files and process headers */
    IN=fopen(argv[optind],   "r");
    IR=fopen(argv[optind+1], "r");
    if (IN==NULL) {fprintf(stderr, "Can't open file %s\n", argv[optind]); exit(0);}
    if (IR==NULL) {fprintf(stderr, "Can't open file %s\n", argv[optind+1]); exit(0);}

    /* pnm type */
    fread(&c,1,1,IN);  /* P */
    fscanf(IN, "%d", &type);
    do{ fread(&c,1,1,IN); } while (c!='\n');  /* \n */

    fread(&c,1,1,IR);  /* P */
    fscanf(IR, "%d", &irtype);
    do{ fread(&c,1,1,IR); } while (c!='\n');  /* \n */

    if (type==5) cw=2; /* grayscale */
    else if (type==6) cw=6; /* color */
    else {fprintf(stderr, "Bad file type: %d\n", type); exit(0);}

    if (irtype!=5) /* grayscale */
      {fprintf(stderr, "Bad IR file type: %d\n", irtype); exit(0);}

    /* pnm comment */
    c=fgetc(IN);
    if (c=='#'){ do{ fread(&c,1,1,IN); } while (c!='\n');}  /* skip comment */
    else ungetc(c, IN);

    c=fgetc(IR);
    if (c=='#'){ do{ fread(&c,1,1,IR); } while (c!='\n');}  /* skip comment */
    else ungetc(c, IR);

    /* pnm size and maxcolor */
    fscanf(IN, "%d %d", &w, &h);
    fscanf(IN, "%d", &mp);
    do{ fread(&c,1,1,IN); } while (c!='\n');  /* \n */

    fscanf(IR, "%d %d", &irw, &irh);
    fscanf(IR, "%d", &irmp);
    do{ fread(&c,1,1,IR); } while (c!='\n');  /* \n */

    if (mp<256) {fprintf(stderr, "Bad number of colors: %d\n", mp); exit(0);}
    if (irw!=w || irh!=h)
        {fprintf(stderr, "IR picture size differs from that of the image: %d x %d vs %d x %d\n", irw,irh,w,h); exit(0);}

    /* allocate buffers and load data */

    buf=(unsigned char *)malloc(h*w*cw);
    if (buf==NULL) {fprintf(stderr, "Can't allocate memory!\n"); exit(0);}

    irbuf=(unsigned char *)malloc(h*w*2);
    if (irbuf==NULL) {fprintf(stderr, "Can't allocate memory!\n"); exit(0);}

    stat=(unsigned char *)malloc(h*w);
    if (stat==NULL) {fprintf(stderr, "Can't allocate memory!\n"); exit(0);}

    if (fread(buf,  cw, w*h, IN)!=w*h) {fprintf(stderr, "can't read the file\n"); exit(0);}
    if (fread(irbuf, 2, w*h, IR)!=w*h) {fprintf(stderr, "can't read IR channel\n"); exit(0);}
    if (IN) fclose(IN);
    if (IR) fclose(IR);
  }

  /***************************************************************************/
  /* 3. Find shift DX,DY between Image and IR channel (todo - scaling!) */
  {
    int rad=5;
    int maxsh=20;
    int nsteps=10;

    int stepx=w/nsteps, stepy=h/nsteps;
    int ii,jj, im,jm, dx,dy;
    double v, vm;
    int n=0;
    DX=0; DY=0;
    int xm[nsteps*nsteps], ym[nsteps*nsteps];

    /* Find squares with large RMS values */
    for (i=0; i<w; i+=stepx){
      for (j=0; j<h; j+=stepy){

        /* Find max rms in one square */
        im=i; jm=j; vm=0;
        for (ii=i+rad; ii<MIN(i+stepx,w-rad); ii++){
          for (jj=j+rad; jj<MIN(j+stepy,h-rad); jj++){
            v=rms_ir(ii,jj, rad);
            if ( v > vm){ vm = v; im=ii; jm=jj; }
          }
        }
        //fprintf(stderr, "rms> %f %d %d\n", vm, im, jm);

        /* compare the square with the red channel
           and find the best shift */
        dx=0; dy=0; vm=-1;
        for (ii=-maxsh; ii<=maxsh; ii++){
          for (jj=-maxsh; jj<=maxsh; jj++){
            v = corr(im, jm, ii,jj, rad);
            if ( v > vm){ vm = v; dx=ii; dy=jj; }
          }
        }
        //fprintf(stderr, "corr> %f %d %d\n", vm, dx, dy);
        DX+=dx;
        DX+=dy;
        n++;

        /* test -- print squares on the image */
        if(0) {
          for (ii=im; ii<im+rad; ii++){
            if (ii+dx<0 || ii+dx>=w) continue;
            if (jm+dy<0 || jm+dy>=h) continue;
            SET_RGB(ii,jm,     0xFFFF,0xFFFF,0xFFFF);
            SET_RGB(ii,jm+rad, 0xFFFF,0xFFFF,0xFFFF);
            SET_RGB(dx+ii,dy+jm,     0xFFFF,0xFFFF,0);
            SET_RGB(dx+ii,dy+jm+rad, 0xFFFF,0xFFFF,0);
          }
          for (jj=jm; jj<jm+rad; jj++){
            if (im+dx<0 || im+dx>=w) continue;
            if (jj+dy<0 || jj+dy>=h) continue;
            SET_RGB(im,jj,     0xFFFF,0xFFFF,0xFFFF);
            SET_RGB(im+rad,jj, 0xFFFF,0xFFFF,0xFFFF);
            SET_RGB(dx+im,dy+jj,     0xFFFF,0xFFFF,0);
            SET_RGB(dx+im+rad,dy+jj, 0xFFFF,0xFFFF,0);
          }
        }
      }
    }
    DX/=n;
    DY/=n;
    fprintf(stderr, "shift> %f %f\n", DX, DY);
  }

  /***************************************************************************/
  /* Remove correlation between IR and other channels
  /* subtract corr(IR,Red)/sigma(Red) * Red*/

  corr2(0, 0, round(DX), round(DY), w+h);

  /***************************************************************************/
  /* interpolation */
  {
    int k, m1,m2,m3,m4,v0,d1,d2,xd,yd;
    int r1,r2,g1,g2,b1,b2;
    int intx,inty;
    for (i=0; i<w; i++){
      for (j=0; j<h; j++){
        /* find left, right, top and bottom maxima */
        m1=m2=m3=m4=0;
        v0=GET_IR(i,j);
        for (k=0; k<irrad; k++){ /* note 1 px additional space for m1..m4*/
          if (i-k-1>=0 && GET_IR(i-k,j)>GET_IR(i-m1,j)) m1=k;
          if (i+k+1<w  && GET_IR(i+k,j)>GET_IR(i+m2,j)) m2=k;
          if (j-k-1>=0 && GET_IR(i,j-k)>GET_IR(i,j-m3)) m3=k;
          if (j+k+1<h  && GET_IR(i,j+k)>GET_IR(i,j+m4)) m4=k;
        }
        /* check do we need interpolation */
        intx=0; inty=0;
        if (GET_IR(i-m1,j)-v0 > irthr && GET_IR(i+m2,j)-v0 > irthr) intx=1;
        if (GET_IR(i,j-m3)-v0 > irthr && GET_IR(i,j+m4)-v0 > irthr) inty=1;

        if (intx==0 && inty==0) continue;
//        if (intx && inty){ /* shortest way */
//          if (m1+m2 > m3+m4) intx=0; else inty=0;
//        }

        //m1++; m2++;m3++;m4++;
        if (intx){
          r1 = GET_R(i-m1,j); r2 = GET_R(i+m2,j);
          g1 = GET_G(i-m1,j); g2 = GET_G(i+m2,j);
          b1 = GET_B(i-m1,j); b2 = GET_B(i+m2,j);
          d1 = m1; d2 = m2;
          for (k=-m1+1; k<m2; k++){
            d1 = m1+k; d2 = m2-k;
            xd=i+k+DX; yd=j+DY;
            if (xd<0 || xd>=w || yd<0 || yd>=h || d1+d2<1) continue;
            SET_RGB(xd,yd, (r1*d2+r2*d1)/(d1+d2),
                           (g1*d2+g2*d1)/(d1+d2),
                           (b1*d2+b2*d1)/(d1+d2));
             //SET_RGB(xd,yd, 0,0xFFFF,0xFFFF);
          }
        }
        if (inty){
          r1 = GET_R(i,j-m3); r2 = GET_R(i,j+m4);
          g1 = GET_G(i,j-m3); g2 = GET_G(i,j+m4);
          b1 = GET_B(i,j-m3); b2 = GET_B(i,j+m4);
          d1 = m3; d2 = m4;
          for (k=-m3+1; k<m4; k++){
            d1 = m3+k; d2 = m4-k;
            xd=i+DX; yd=j+k+DY;
            if (xd<0 || xd>=w || yd<0 || yd>=h || d1+d2<1) continue;
              SET_RGB(xd,yd, (r1*d2+r2*d1)/(d1+d2),
                           (g1*d2+g2*d1)/(d1+d2),
                           (b1*d2+b2*d1)/(d1+d2));
            //SET_RGB(xd,yd, 0xFFFF,0xFFFF,0);
          }
        }
        if (inty || intx){
//          xd=i+DX; yd=j+DY;
//          if (xd>=0 && xd<w && yd>=0 && yd<h){
//            SET_RGB(xd,yd, (r1*d2+r2*d1)/(d1+d2),
//                           (g1*d2+g2*d1)/(d1+d2),
//                           (b1*d2+b2*d1)/(d1+d2));
//          }
          SET_IR(i,j, 0x0);
        }
      }
    }
  }


//      xd=x+dx; yd=y+dy;
//      if (xd<0 || xd>=w || yd<0 || yd>=h) continue;


  /***************************************************************************/
  /* Save the data*/

  printf("P%d\n%d %d\n%d\n",cw==2?5:6,w,h,256*256-1);
  if (fwrite(buf, cw, w*h, stdout)!=w*h) {fprintf(stderr, "Write Error!\n"); exit(0);}

  if (outd){
    FILE *OUT = fopen(outd, "w");
    if (OUT==NULL){
      fprintf(stderr, "can't open file: %s\n", outd);
      exit(0);
    }
    fprintf(OUT, "P%d\n%d %d\n%d\n",5,w,h,256*256-1);
    if (fwrite(irbuf, 2, w*h, OUT)!=w*h) {fprintf(stderr, "Write Error!\n"); exit(0);}
    fclose(OUT);
  }

}
//