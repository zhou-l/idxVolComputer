#ifndef COLOR_H_
#define COLOR_H_

class Color
{
public:
	double r() const{ return _r;}
	double g() const{ return _g;}
	double b() const{ return _b;}

	void set_r(double r){ _r = r; }
	void set_g(double g){ _g = g; }
	void set_b(double b){ _b = b; }

	Color(void):_r(0.0f),_g(0.0f),_b(0.0f){};
	~Color(void){};
	Color(const Color& cz):_r(cz._r),_g(cz._g),_b(cz._b){}

	Color( double r, double g, double b) :_r(r),_g(g), _b(b){}

	Color& operator = (const Color& v)
	{
		_r = v._r;
		_g = v._g;
		_b = v._b;
		return *this;
	}

	// scalar product
	Color operator *(double s) const
	{
		return Color( _r * s, _g * s, _b * s);
	}
	// color product
	Color operator * (Color c) const
	{
		return Color( _r * c._r, _g * c._g, _b * c._b );
	}
	// color plus
	Color operator +(const Color& c) const
	{
		return Color( _r + c._r, _g + c._g, _b + c._b);
	}
	// color minus
	Color operator -(const Color& c) const
	{
		return Color( _r - c._r, _g - c._g, _b - c._b);
	}

	friend inline double Max( const Color& c) 
	{
		double t = (c._r > c._g)?c._r:c._g;
		return (t > c._b)?t:c._b; 
	}

	friend inline double Min( const Color& c) 
	{
		double t = (c._r < c._g)?c._r:c._g;
		return (t < c._b)?t:c._b; 
	}


	// scalar left product
	friend inline Color operator * ( double s, const Color& c)
	{
		return Color( s * c._r, s * c._g, s * c._b );
	}

private:
	double _r,_g,_b;
};

#endif