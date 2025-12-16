#ifndef Image_h
#define Image_h

#include "Color.h"
#include <math.h>
#include <string>

struct Pixel{
         unsigned char r, g, b;
};

class Image {
public:
	Image(const std::string& name);
	Image(int xres, int yres);
	~Image();
	void set(int x, int y, const Color& c) {
	        Pixel p;
                if(c.r() < 0.0f)
                   p.r = 0;
                else
                   p.r = unsigned char(c.r() >= 1.0f? 255 : c.r() * 255.0f);
                if(c.g() < 0.0f)
                   p.g = 0;
                else
                   p.g = unsigned char(c.g() >= 1.0f? 255 : c.g() * 255.0f);
                if(c.b() < 0.0f)
                   p.b = 0;
                else
                   p.b = unsigned char(c.b() >= 1.0f? 255 : c.b() * 255.0f);
		data[y][x] = p;
	}
	void write(const std::string& filename) const;
	double aspect_ratio() const {
		return double(xres)/double(yres);
	}
	int getXresolution() {
		return xres;
	}
	int getYresolution() {
		return yres;
	}
	Color interpolate(double x, double y) const {
		x *= xres; y *= yres;
		int ix = int(floor(x))%xres;
		if(ix<0)
			ix += xres;
		int ix1 = (ix+1)%xres;
		int iy = int(floor(y))%yres;
		if(iy<0)
			iy += yres;
		int iy1 = (iy+1)%yres;
		double fx = x-ix;
		double fy = y-iy;

		Color c00 = Color(data[iy][ix].r,   data[iy][ix].g,   data[iy][ix].b);
		Color c01 = Color(data[iy][ix1].r,  data[iy][ix1].g,  data[iy][ix1].b);
		Color c10 = Color(data[iy1][ix].r,  data[iy1][ix].g,  data[iy1][ix].b);
		Color c11 = Color(data[iy1][ix1].r, data[iy1][ix1].g, data[iy1][ix1].b);
		Color c = c00*(1-fx)*(1-fy) + c01*fx*(1-fy) + c10*(1-fx)*fy + c11*fx*fy;
		return c*(1.0f/255);
	}

	Pixel get(int x, int y)
	{
		Pixel p = data[y][x];
		return p;
	}

protected:
	Pixel** data;
	int xres, yres;
	Image(const Image&);
	Image& operator=(const Image&);
};

#endif
