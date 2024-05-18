#ifndef PNM_H
#define PNM_H

#include <cstring>
#include <cassert>
#include <vector>

/********************************************************************/
/* a 2D point (simplified mapsoft/iPoint class) */
struct PT{
  int x,y;
  PT(int _x, int _y): x(_x), y(_y) { }
  PT(): x(0), y(0) { }

  bool operator< (const PT & other) const {
    return (x<other.x) || ((x==other.x)&&(y<other.y));
  }
  bool operator== (const PT & other) const {
    return (x==other.x)&&(y==other.y);
  }

  PT adj(const int dir) const;
  int is_adj(const PT & p) const;
};
typedef std::vector<PT> PTS;


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

  void calc_mmm(int &min, int &mean, int &max, int ch) const;
  int is_in(const PT &p) const {
    return p.x>=0 && p.x<w && p.y>=0 && p.y<h; }

  /********************************************************************/
  /* get function (no internal checks). */
  int get(int ch, const PT &p) const{
    if (bpp==6 || bpp==2)
      return buf[bpp*(p.y*w+p.x)+2*ch]*256
           + buf[bpp*(p.y*w+p.x)+2*ch+1];
    if (bpp==3 || bpp==1)
      return buf[bpp*(p.y*w+p.x)+ch];
    return -1;
  }

  /* set function (no internal checks). */
  void set(int ch, const PT &p, int v){
    if (bpp==6 || bpp==2){
      buf[bpp*(p.y*w+p.x)+2*ch]   = (v>>8)&0xFF;
      buf[bpp*(p.y*w+p.x)+2*ch+1] = v&0xFF; }
    if (bpp==3 || bpp==1){
      buf[bpp*(p.y*w+p.x)+ch] = v&0xFF; }
  }

  /* setrgb function (no internal checks). */
  void setrgb(const PT &p, int r, int g, int b){
    if (bpp==6){
      buf[bpp*(p.y*w+p.x)+0] = (r>>8)&0xFF;
      buf[bpp*(p.y*w+p.x)+1] =  r&0xFF;
      buf[bpp*(p.y*w+p.x)+2] = (g>>8)&0xFF;
      buf[bpp*(p.y*w+p.x)+3] =  g&0xFF;
      buf[bpp*(p.y*w+p.x)+4] = (b>>8)&0xFF;
      buf[bpp*(p.y*w+p.x)+5] =  b&0xFF;
      return;
    }
    if (bpp==3){
      buf[bpp*(p.y*w+p.x)+0] = r&0xFF;
      buf[bpp*(p.y*w+p.x)+1] = g&0xFF;
      buf[bpp*(p.y*w+p.x)+2] = b&0xFF;
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
