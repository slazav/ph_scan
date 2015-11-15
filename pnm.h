#ifndef PNM_H
#define PNM_H

#include <cstring>
#include <cassert>

/* a simple library for working with pnm files */
/********************************************************************/

class PNM {
  public:
    int w, h, bpp;  // width, height, bytes per point (1,2,3 or 6)
    unsigned char *buf; // raw data (w*h*bpp bytes)

  public:

  /* Load from a file */
  PNM(const char *fname);

  /* Create new */
  PNM(int w, int h, int bpp) {
    create(w,h,bpp);
    memset(buf, 0, w*h*bpp);
  }

  /* Copy constructor + assignment */
  PNM(const PNM & other) { copy(other); }
  PNM & operator= (const PNM & other){
    if (&other == this) return *this;
    destroy();
    copy(other);
    return *this;
  }
  /* Destuctor */
  ~PNM() {destroy();}

  int get_type() const;
  int get_mcol() const;
  int is_rgb() const {return get_type()==6;}
  int save(const char *fname) const;

  int calc_mmm(int &min, int &mean, int &max, int ch) const;

  /********************************************************************/
  /* get/set functions */
  /* read functions for 8 and 16 bit data (no internal checks). */
  int get8(int ch, int x, int y) const{
    return buf[bpp*(y*w+x)+ch]; }
  int get16(int ch, int x, int y) const{
    return buf[bpp*(y*w+x)+2*ch]*256
         + buf[bpp*(y*w+x)+2*ch+1]; }

  int get(int ch, int x, int y) const{
    if (bpp==6 || bpp==2)
      return buf[bpp*(y*w+x)+2*ch]*256
           + buf[bpp*(y*w+x)+2*ch+1];
    if (bpp==3 || bpp==1)
      return buf[bpp*(y*w+x)+ch];
    return -1;
  }

  /* write functions for 8 and 16 bit data (no internal checks). */
  void set8(int ch, int x, int y, int v){
    buf[bpp*(y*w+x)+ch] = v&0xFF; }
  void set16(int ch, int x, int y, int v){
    buf[bpp*(y*w+x)+2*ch]   = (v>>8)&0xFF;
    buf[bpp*(y*w+x)+2*ch+1] = v&0xFF; }
  void set(int ch, int x, int y, int v){
    if (bpp==6 || bpp==2){
      buf[bpp*(y*w+x)+2*ch]   = (v>>8)&0xFF;
      buf[bpp*(y*w+x)+2*ch+1] = v&0xFF; }
    if (bpp==3 || bpp==1){
      buf[bpp*(y*w+x)+ch] = v&0xFF; }
  }

  void set8rgb(int x, int y, int r, int g, int b){
    buf[bpp*(y*w+x)+0] = r&0xFF;
    buf[bpp*(y*w+x)+1] = g&0xFF;
    buf[bpp*(y*w+x)+2] = b&0xFF;
  }
  void set16rgb(int x, int y, int r, int g, int b){
    buf[bpp*(y*w+x)+0] = (r>>8)&0xFF;
    buf[bpp*(y*w+x)+1] =  r&0xFF;
    buf[bpp*(y*w+x)+2] = (g>>8)&0xFF;
    buf[bpp*(y*w+x)+3] =  g&0xFF;
    buf[bpp*(y*w+x)+4] = (b>>8)&0xFF;
    buf[bpp*(y*w+x)+5] =  b&0xFF;
  }
  void setrgb(int x, int y, int r, int g, int b){
    if (bpp==6){
      buf[bpp*(y*w+x)+0] = (r>>8)&0xFF;
      buf[bpp*(y*w+x)+1] =  r&0xFF;
      buf[bpp*(y*w+x)+2] = (g>>8)&0xFF;
      buf[bpp*(y*w+x)+3] =  g&0xFF;
      buf[bpp*(y*w+x)+4] = (b>>8)&0xFF;
      buf[bpp*(y*w+x)+5] =  b&0xFF;
      return;
    }
    if (bpp==3){
      buf[bpp*(y*w+x)+0] = r&0xFF;
      buf[bpp*(y*w+x)+1] = g&0xFF;
      buf[bpp*(y*w+x)+2] = b&0xFF;
      return;
    }
  }

  /********************************************************************/
  private:
    int *refcounter;
    /// create
    void create(int _w, int _h, int _bpp){
      w=_w; h=_h; bpp=_bpp;
      assert(w>=0 && h>=0);
      assert(bpp==1 || bpp==2 || bpp==3 || bpp==6);
      buf = new unsigned char[w*h*bpp];
      assert(buf);
      refcounter   = new int;
      *refcounter  = 1;
    }
    /// copy
    void copy(const PNM & other){
      w=other.w; h=other.h; bpp = other.bpp;
      buf = other.buf;
      refcounter = other.refcounter;
      (*refcounter)++;
      assert(*refcounter>0);
    }
    /// destroy image
    void destroy(void){
      (*refcounter)--;
      if (*refcounter<=0){
        delete[] buf;
        delete refcounter;
      }
    }
};

#endif
