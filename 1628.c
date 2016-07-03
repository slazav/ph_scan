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


int
main(int argc, char *argv[]){

  unsigned char *buf, *irbuf;
  unsigned char c;
  unsigned w,h,mp;
  int      i,j,n;
  int      dh=100;
  int      neg=1;
  int      brd=200;
  double   lth = 1, hth = 0; // thresholds, per cent
  int      aval[3]={130,140,150};
  int      cw=6, type;
  int      l=100,d=100;
  unsigned max[3], min[3];
  double   avr[3];
  double A[3], B[3], C[3];
  long int H[3][0x10000]; //  histogram

  FILE     *IN=NULL;
  long     pos;

  void usage(void){
    fprintf(stderr, "usage: %s OPTIONS in_file > out_file\n", argv[0]);
    fprintf(stderr, "options:   -b #  -- border (in pixels)  -- default 200\n");
    fprintf(stderr, "           -R #  -- A-value for R 1-254 -- default 64\n");
    fprintf(stderr, "           -G #  -- A-value for G 1-254 -- default 64\n");
    fprintf(stderr, "           -B #  -- A-value for B 1-254 -- default 64\n");
    fprintf(stderr, "           -L #  -- low threshold, %% -- default 1.0\n");
    fprintf(stderr, "           -H #  -- high threshold, %% -- default 0.0\n");
    fprintf(stderr, "           -n    -- input is negative (default)\n");
    fprintf(stderr, "           -p    -- input is positive\n");
    exit(0);
  }

  while (1){
    i=getopt(argc, argv, "R:G:B:b:L:H:pn");
    if (i==-1) break;
    switch(i){
      case 'R': 
        aval[0]=atoi(optarg); 
        if ((aval[0]<1)||(aval[0]>254)) usage();
        break;
      case 'G': 
        aval[1]=atoi(optarg); 
        if ((aval[1]<1)||(aval[1]>254)) usage();
        break;
      case 'B': 
        aval[2]=atoi(optarg); 
        if ((aval[2]<1)||(aval[2]>254)) usage();
        break;
      case 'L': lth=atof(optarg); if ((lth<0)||(lth>=100)) usage(); break;
      case 'H': hth=atof(optarg); if ((hth<0)||(hth>=100)) usage(); break;
      case 'b': brd=atoi(optarg); break;
      case 'n': neg=1; break;
      case 'p': neg=0; break;
      default:  usage();
    }
  }

  if (argc-optind!=1) usage();

  /************************************************/
  /* read file information */

  /* Process PNM header*/
  IN=fopen(argv[optind], "r");
  if (IN==NULL) {fprintf(stderr, "Can't open file %s\n", argv[optind]); exit(0);}

  fread(&c,1,1,IN);  /* P */
  fscanf(IN, "%d", &type);
  do{ fread(&c,1,1,IN); } while (c!='\n');  /* \n */

  c=fgetc(IN);
  if (c=='#'){ do{ fread(&c,1,1,IN); } while (c!='\n');}  /* skip comment */
  else ungetc(c, IN);

  fscanf(IN, "%d %d", &w, &h);
  fscanf(IN, "%d", &mp);
  do{ fread(&c,1,1,IN); } while (c!='\n');  /* \n */

  if (type==5) cw=2; /* grayscale */
  else if (type==6) cw=6; /* color */
  else {fprintf(stderr, "Bad file type: %d\n", type); exit(0);}

  if (mp!=0xFFFF) {fprintf(stderr, "Bad number of colors: %d\n", mp); exit(0);}
  buf=(unsigned char *)malloc(w*cw);
  if (buf==NULL) {fprintf(stderr, "Can't allocate memory!\n"); exit(0);}
  pos = ftell(IN);


  /************************************************/
  /* build histogram */

  memset(H, 0, sizeof(H));
  n = 0;
  if (h-2*brd-1<=0 || w-2*brd-1<=0){
    fprintf(stderr, "Too large border!\n"); exit(0); }

  for (i=0;i<h;i++){
    if (fread(buf, cw, w, IN)!=w) {fprintf(stderr, "Read Error!\n"); exit(0);}
    if (i<brd) continue;
    if (i>=h-brd) break;
    for (j=brd;j<w-brd;j++){
      if (cw==2){
        int col=get16(buf+cw*j);
        H[0][col&0xFFFF]++;
      }
      else {
        int r=get16(buf+cw*j);
        int g=get16(buf+cw*j+2);
        int b=get16(buf+cw*j+4);
        H[0][r&0xFFFF]++;
        H[1][g&0xFFFF]++;
        H[2][b&0xFFFF]++;
      }
      n++;
    }
  }

  /************************************************/
  /* Calculate MIN/MAX/AVRG for color values */

  {
    for (i=0;i<3;i++){max[i]=mp-1; min[i]=0; avr[i]=0;}
    long int sr1=0, sg1=0, sb1=0;
    long int sr2=0, sg2=0, sb2=0;
    lth*=n/100.0;
    hth*=n/100.0;
    for (i=0; i<0x10000; i++){
      j = mp-i+1;
      if (cw==2){
        sr1+=H[0][i];
        sr2+=H[0][j];
        if (sr1<lth+1) min[0]=i;
        if (sr2<hth+1) max[0]=j;
        avr[0]+=(double)(H[0][i])*i;
      }
      else {
        sr1+=H[0][i];
        sr2+=H[0][j];
        sg1+=H[1][i];
        sg2+=H[1][j];
        sb1+=H[2][i];
        sb2+=H[2][j];
        if (sr1<lth+1) min[0]=i;
        if (sr2<hth+1) max[0]=j;
        if (sg1<lth+1) min[1]=i;
        if (sg2<hth+1) max[1]=j;
        if (sb1<lth+1) min[2]=i;
        if (sb2<hth+1) max[2]=j;
        avr[0]+=(double)H[0][i]*i;
        avr[1]+=(double)H[1][i]*i;
        avr[2]+=(double)H[2][i]*i;
      }
    }
    avr[0]/=n;
    avr[1]/=n;
    avr[2]/=n;
  }

//  fprintf(stderr, "%d %d %d  %d %d %d  %d %d %d\n",
//          min[0], (unsigned)avr[0], max[0],
//          min[1], (unsigned)avr[1], max[1],
//          min[2], (unsigned)avr[2], max[2]);


  /************************************************/
  /* return back and do the color correction */
  fseek(IN, pos, SEEK_SET);

  for (i=0; i<cw/2;i++){
    B[i]=(255.0*max[i]/(max[i]-min[i]) - 1.0*aval[i]*avr[i]/(avr[i]-min[i])) /
         (255.0/(max[i]-min[i]) - 1.0*aval[i]/(avr[i]-min[i]));
    A[i]=255.0*(max[i]-B[i])/(max[i]-min[i]);
    C[i]=A[i]*(B[i]-min[i]);
  }

  printf("P%d\n%d %d\n255\n",type,w,h);

  for (i=0;i<h;i++){
    if (fread(buf, cw, w, IN)!=w) {fprintf(stderr, "can't read the file\n"); exit(0);}

    for (j=0;j<w;j++){
      if (cw==2){
        int col=get16(buf+cw*j);
        col=A[0]-C[0]/(B[0]-col);
        if (neg!=0) { col=255-col; }
        buf[cw/2*j]=LIMIT(col,0,255);
      }
      else {
        int r,g,b;
        r=get16(buf+cw*j);
        g=get16(buf+cw*j+2);
        b=get16(buf+cw*j+4);
        r=A[0]-C[0]/(B[0]-r);
        g=A[1]-C[1]/(B[1]-g);
        b=A[2]-C[2]/(B[2]-b);

        if (neg!=0) { r=255-r; g=255-g; b=255-b; }

        buf[cw/2*j]=LIMIT(r,0,255);
        buf[cw/2*j+1]=LIMIT(g,0,255);
        buf[cw/2*j+2]=LIMIT(b,0,255);
      }
    }
    if (fwrite(buf, cw/2, w, stdout)!=w) {fprintf(stderr, "Write Error!\n"); exit(0);}
  }
  if (IN) fclose(IN);
  exit(0);
}


