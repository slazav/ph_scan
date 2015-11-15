#include "pnm.h"
#include <cstdio>

/**************************************************************************/
/* load a file*/
PNM::PNM(const char *fname){
  int type, mcol;
  char c;
  FILE *IN;

  /* open file ant read the header */
  IN=fopen(fname, "r");
  if (IN==NULL) {
    fprintf(stderr, "Can't open file: %s\n", fname);
    create(0,0,1); return;
  }

  fread(&c,1,1,IN);  /* P */
  fscanf(IN, "%d", &type);
  do{ fread(&c,1,1,IN); } while (c!='\n');  /* \n */

  c=fgetc(IN);
  if (c=='#'){
    do{ fread(&c,1,1,IN); } while (c!='\n');}  /* skip comment */
  else ungetc(c, IN);

  fscanf(IN, "%d %d", &w, &h);
  fscanf(IN, "%d", &mcol);
  do{ fread(&c,1,1,IN); } while (c!='\n');  /* \n */

  if      (type==5 && mcol==255) bpp=1; /*  8bit grayscale */
  else if (type==6 && mcol==255) bpp=3; /* 24bit color */
  else if (type==5 && mcol>255)  bpp=2; /* 16bit grayscale */
  else if (type==6 && mcol>255)  bpp=6; /* 48bit color */
  else {
    fprintf(stderr, "Unsupported file type: %s: %d\n", fname, type);
    create(0,0,1); return;
  }

  create(w,h,bpp);

  /* read data */
  if (fread(buf, bpp, w*h, IN)!=w*h){
    fprintf(stderr, "Can't read the file: %s\n", fname);
    fclose(IN);  create(0,0,1); return;
  }

  /* close the file */
  fclose(IN);
}

int
PNM::get_type() const{
  if (bpp==1 || bpp==2) return 5; /* grayscale */
  if (bpp==3 || bpp==6) return 6; /* RGB       */
  fprintf(stderr, "Unsupported bpp: %d\n", bpp);
  return 0;
}
int
PNM::get_mcol() const{
  if (bpp==1 || bpp==3) return (2<<7)-1;  /* 8 bit  */
  if (bpp==2 || bpp==6) return (2<<15)-1; /* 16 bit */
  fprintf(stderr, "Unsupported bpp: %d\n", bpp);
  return 0;
}

int
PNM::save(const char *fname) const{
  FILE *OUT = fopen(fname, "w");
  if (OUT==NULL){
    fprintf(stderr, "can't open file: %s\n", fname);
    return 1;
  }
  fprintf(OUT, "P%d\n%d %d\n%d\n",get_type(),w,h,get_mcol());
    if (fwrite(buf, bpp, w*h, OUT)!=w*h)
      {fprintf(stderr, "Write Error!\n"); return 1;}
    fclose(OUT);
  return 0;
}

int
PNM::calc_mmm(int &min, int &mean, int &max, int ch) const{
  int x,y,v;
  long long sum=0;
  min=get_mcol(); max=0;
  for (x=0; x<w; x++){
    for (y=0; y<h; y++){
      v=get(ch,x,y);
      sum += v;
      if (min>v) min=v;
      if (max<v) max=v;
    }
  }
  mean = sum/w/h;
}
