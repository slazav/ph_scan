#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>


#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define LIMIT(x,min,max) MIN(MAX(x, min), max)


int get16(unsigned char *buf){
  return( buf[0]*256 + buf[1] );
}
int set16(unsigned char *buf, int v){
  buf[0] = (v>>8)&0xFF;
  buf[1] = (v&0xFF);
}


int
main(int argc, char *argv[]){

  unsigned char *buf, *irbuf, *stat;
  unsigned char c;
  int      i,j;
  int      cw=6, type, irtype;
  unsigned w,h,mp, irw,irh,irmp;
  FILE     *IN=NULL, *IR=NULL;
  int irrad=30;
  int irthr=round(0.06*256*256);

  void usage(void){
    fprintf(stderr, "usage: %s OPTIONS in_file ir_file > out_file\n", argv[0]);
    fprintf(stderr, "           -T name -- infrared threshold 0 .. 1.0 -- default 0.06\n");
    fprintf(stderr, "           -R name -- interpolation radius, px         -- default 30\n");
    exit(0);
  }

  while (1){
    i=getopt(argc, argv, "T:R:");
    if (i==-1) break;
    switch(i){
      case 'T': irthr=round(256*256*atof(optarg)); break;
      case 'R': irrad=atoi(optarg); break;
      default:  usage();
    }
  }

  if (argc-optind!=2) usage();

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

  /* allocate buffers */

  buf=(unsigned char *)malloc(w*cw);
  if (buf==NULL) {fprintf(stderr, "Can't allocate memory!\n"); exit(0);}

  irbuf=(unsigned char *)malloc(w*2);
  if (irbuf==NULL) {fprintf(stderr, "Can't allocate memory!\n"); exit(0);}

  stat=(unsigned char *)malloc(w);
  if (stat==NULL) {fprintf(stderr, "Can't allocate memory!\n"); exit(0);}


  printf("P%d\n%d %d\n%d\n",type,w,h,mp);

  /* infrared correction */
  for (i=0;i<h;i++){
    int k,n,s,sn;

    /* read data */
    if (fread(buf, cw, w, IN)!=w) {fprintf(stderr, "can't read the file\n"); exit(0);}
    if (fread(irbuf, 2, w, IR)!=w) {fprintf(stderr, "can't read IR channel\n"); exit(0);}

    /* find points for interpolation */
    for (j=1; j<w; j++){
      /* calculate s, mean value of previous irrad non-zero points */
      k=j; s=0; sn=0;
      for (n=0;n<irrad; n++){
        k--; if (k<0) break;
        if (stat[k]==0) continue;
        s+=get16(irbuf+2*k); sn++;
      }
      s = sn? s/sn : 0;

      /* put the point to 0 if it differs from s more then irthr */
      if (s-get16(irbuf+2*j) > irthr) stat[j]=0;
      else stat[j]=1;
    }

    /* extend the interpolation range by n points */
    for (n=0;n<2;n++){
      for (j=1; j<w; j++){    if (stat[j]==0) stat[j-1]=0; }
      for (j=w-2; j>=0; j--){ if (stat[j]==0) stat[j+1]=0; }
    }

    /* interpolation */
    for (j=1; j<w; j++){
      int n1m=j, n1p=j, n2m=j, n2p=j;

      /* skip points which do not need interpolation*/
      if (stat[j]) continue;

      /* find boundaries of interpolation area, n2m,n2p*/
      /* and indices of outer points, n1m,n1p to get colors */
      for (k=j; k>0 && k>j-irrad; k--){
        if (stat[k]) { n1m=k; break; }
        else n2m=k;
      }
      for (k=j; k<w && k<j+irrad; k++){
        if (stat[k]){ n1p=k; break; }
        else n2p=k;
      }
      if (n1m==j) n1m=n1p;
      if (n1p==j) n1p=n1m;

      /* do the interpolation */
      if (cw==6){
        int r,g,b, dm,dp;
        int rm = get16(buf+cw*n1m);
        int rp = get16(buf+cw*n1p);
        int gm = get16(buf+cw*n1m+2);
        int gp = get16(buf+cw*n1p+2);
        int bm = get16(buf+cw*n1m+4);
        int bp = get16(buf+cw*n1p+4);
        for (k=n2m; k<=n2p; k++){
          dm=k-n2m;
          dp=n2p-k;
          r = (rm*dp + rp*dm)/(dm+dp+1);
          g = (gm*dp + gp*dm)/(dm+dp+1);
          b = (bm*dp + bp*dm)/(dm+dp+1);

          stat[k]=1;
          set16(buf+cw*k+0, r);
          set16(buf+cw*k+2, g);
          set16(buf+cw*k+4, b);
        }
      }
      if (cw==2){
        int col,dm,dp;
        int cm = get16(buf+cw*n1p);
        int cp = get16(buf+cw*n1p);
        for (k=n2m; k<=n2p; k++){
          dm=k-n2m+1;
          dp=n2p-k+1;
          col = (cm*dp + cp*dm)/(dm+dp);
          stat[k]=1;
          set16(buf+cw*k, col);
        }
      }
    }

    if (fwrite(buf, cw, w, stdout)!=w) {fprintf(stderr, "Write Error!\n"); exit(0);}
  }
  if (IN) fclose(IN);
  if (IR) fclose(IR);
  exit(0);
}


