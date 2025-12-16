
/*
 * pbrt source code Copyright(c) 1998-2004 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_GEOMETRY_H
#define PBRT_GEOMETRY_H
// geometry.h*
#include <float.h>
#include "global.h"
// Geometry Declarations
class Vector {
public:
	// Vector Public Methods
	Vector(float _x=0, float _y=0, float _z=0)
		: x(_x), y(_y), z(_z) {
	}
	explicit Vector(const Point &p);
	Vector operator+(const Vector &v) const {
		return Vector(x + v.x, y + v.y, z + v.z);
	}
	
	Vector& operator+=(const Vector &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	Vector operator-(const Vector &v) const {
		return Vector(x - v.x, y - v.y, z - v.z);
	}
	
	Vector& operator-=(const Vector &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	bool operator==(const Vector &v) const {
		return x == v.x && y == v.y && z == v.z;
	}
	Vector operator*(float f) const {
		return Vector(f*x, f*y, f*z);
	}
	
	Vector &operator*=(float f) {
		x *= f; y *= f; z *= f;
		return *this;
	}
	Vector operator/(float f) const {
		assert(f!=0);
		float inv = 1.f / f;
		return Vector(x * inv, y * inv, z * inv);
	}
	
	Vector &operator/=(float f) {
		assert(f!=0);
		float inv = 1.f / f;
		x *= inv; y *= inv; z *= inv;
		return *this;
	}
	Vector operator-() const {
		return Vector(-x, -y, -z);
	}
	float operator[](int i) const {
		assert(i >= 0 && i <= 2);
		return (&x)[i];
	}
	
	float &operator[](int i) {
		assert(i >= 0 && i <= 2);
		return (&x)[i];
	}
	float LengthSquared() const { return x*x + y*y + z*z; }
	float Length() const { return sqrtf(LengthSquared()); }
	explicit Vector(const Normal &n);
	// Vector Public Data
	float x, y, z;
};
class Point {
public:
	// Point Methods
	Point(float _x=0, float _y=0, float _z=0)
		: x(_x), y(_y), z(_z) {
	}
	Point operator+(const Vector &v) const {
		return Point(x + v.x, y + v.y, z + v.z);
	}
	
	Point &operator+=(const Vector &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	Vector operator-(const Point &p) const {
		return Vector(x - p.x, y - p.y, z - p.z);
	}
	
	Point operator-(const Vector &v) const {
		return Point(x - v.x, y - v.y, z - v.z);
	}
	
	Point &operator-=(const Vector &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	Point &operator+=(const Point &p) {
		x += p.x; y += p.y; z += p.z;
		return *this;
	}
	Point operator+(const Point &p) {
		return Point(x + p.x, y + p.y, z + p.z);
	}
	Point operator* (float f) const {
		return Point(f*x, f*y, f*z);
	}
	Point &operator*=(float f) {
		x *= f; y *= f; z *= f;
		return *this;
	}
	Point operator/ (float f) const {
		float inv = 1.f/f;
		return Point(inv*x, inv*y, inv*z);
	}
	Point &operator/=(float f) {
		float inv = 1.f/f;
		x *= inv; y *= inv; z *= inv;
		return *this;
	}
	float operator[](int i) const { return (&x)[i]; }
	float &operator[](int i) { return (&x)[i]; }
	// Point Public Data
	float x,y,z;
};


class Segment{
public:
	Point P0, P1;
};

class  Normal {
public:
	// Normal Methods
	Normal(float _x=0, float _y=0, float _z=0)
		: x(_x), y(_y), z(_z) {}
	Normal operator-() const {
		return Normal(-x, -y, -z);
	}
	Normal operator+ (const Normal &v) const {
		return Normal(x + v.x, y + v.y, z + v.z);
	}
	
	Normal& operator+=(const Normal &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	Normal operator- (const Normal &v) const {
		return Normal(x - v.x, y - v.y, z - v.z);
	}
	
	Normal& operator-=(const Normal &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	Normal operator* (float f) const {
		return Normal(f*x, f*y, f*z);
	}
	
	Normal &operator*=(float f) {
		x *= f; y *= f; z *= f;
		return *this;
	}
	Normal operator/ (float f) const {
		float inv = 1.f/f;
		return Normal(x * inv, y * inv, z * inv);
	}
	
	Normal &operator/=(float f) {
		float inv = 1.f/f;
		x *= inv; y *= inv; z *= inv;
		return *this;
	}
	float LengthSquared() const { return x*x + y*y + z*z; }
	float Length() const        { return sqrtf(LengthSquared()); }
	
	explicit Normal(const Vector &v)
	  : x(v.x), y(v.y), z(v.z) {}
	float operator[](int i) const { return (&x)[i]; }
	float &operator[](int i) { return (&x)[i]; }
	// Normal Public Data
	float x,y,z;
};
class  Ray {
public:
	// Ray Public Methods
	Ray(): mint(RAY_EPSILON), maxt(INFINITY), time(0.f) {}
	Ray(const Point &origin, const Vector &direction,
		float start = RAY_EPSILON, float end = INFINITY, float t = 0.f)
		: o(origin), d(direction), mint(start), maxt(end), time(t) { }
	Point operator()(float t) const { return o + d * t; }
	// Ray Public Data
	Point o;
	Vector d;
	mutable float mint, maxt;
	float time;
};
class RayDifferential : public Ray {
public:
	// RayDifferential Methods
	RayDifferential() { hasDifferentials = false; }
	RayDifferential(const Point &org, const Vector &dir)
			: Ray(org, dir) {
		hasDifferentials = false;
	}
	explicit RayDifferential(const Ray &ray) : Ray(ray) {
		hasDifferentials = false;
	}
	// RayDifferential Public Data
	bool hasDifferentials;
	Ray rx, ry;
};
class BBox {
public:
	// BBox Public Methods
	BBox() {
		pMin = Point( INFINITY,  INFINITY,  INFINITY);
		pMax = Point(-INFINITY, -INFINITY, -INFINITY);
	}
	BBox(const Point &p) : pMin(p), pMax(p) { }
	BBox(const Point &p1, const Point &p2) {
		pMin = Point(min(p1.x, p2.x),
					 min(p1.y, p2.y),
					 min(p1.z, p2.z));
		pMax = Point(max(p1.x, p2.x),
					 max(p1.y, p2.y),
					 max(p1.z, p2.z));
	}
	friend inline ostream &
		operator<<(ostream &os, const BBox &b);
	friend  BBox Union(const BBox &b, const Point &p);
	friend  BBox Union(const BBox &b, const BBox &b2);
	bool Overlaps(const BBox &b) const {
		bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
		bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
		bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
		return (x && y && z);
	}
	bool Inside(const Point &pt) const {
		return (pt.x >= pMin.x && pt.x <= pMax.x &&
	            pt.y >= pMin.y && pt.y <= pMax.y &&
	            pt.z >= pMin.z && pt.z <= pMax.z);
	}
	void Expand(float delta) {
		pMin -= Vector(delta, delta, delta);
		pMax += Vector(delta, delta, delta);
	}
	float Volume() const {
		Vector d = pMax - pMin;
		return d.x * d.y * d.z;
	}
	int MaximumExtent() const {
		Vector diag = pMax - pMin;
		if (diag.x > diag.y && diag.x > diag.z)
			return 0;
		else if (diag.y > diag.z)
			return 1;
		else
			return 2;
	}
	void BoundingSphere(Point *c, float *rad) const;
	bool IntersectP(const Ray &ray,
	                float *hitt0 = NULL,
					float *hitt1 = NULL) const;
	// BBox Public Data
	Point pMin, pMax;
};
// Geometry Inline Functions
inline Vector::Vector(const Point &p)
	: x(p.x), y(p.y), z(p.z) {
}
inline ostream &operator<<(ostream &os, const Vector &v) {
	os << v.x << ", " << v.y << ", " << v.z;
	return os;
}
inline Vector operator*(float f, const Vector &v) {
	return v*f;
}
inline float Dot(const Vector &v1, const Vector &v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
inline float AbsDot(const Vector &v1, const Vector &v2) {
	return fabsf(Dot(v1, v2));
}
inline Vector Cross(const Vector &v1, const Vector &v2) {
	return Vector((v1.y * v2.z) - (v1.z * v2.y),
                  (v1.z * v2.x) - (v1.x * v2.z),
                  (v1.x * v2.y) - (v1.y * v2.x));
}
inline Vector Cross(const Vector &v1, const Normal &v2) {
	return Vector((v1.y * v2.z) - (v1.z * v2.y),
                  (v1.z * v2.x) - (v1.x * v2.z),
                  (v1.x * v2.y) - (v1.y * v2.x));
}
inline Vector Cross(const Normal &v1, const Vector &v2) {
	return Vector((v1.y * v2.z) - (v1.z * v2.y),
                  (v1.z * v2.x) - (v1.x * v2.z),
                  (v1.x * v2.y) - (v1.y * v2.x));
}
inline Vector Normalize(const Vector &v) {
	return v / v.Length();
}
inline void CoordinateSystem(const Vector &v1, Vector *v2, Vector *v3) {
	if (fabsf(v1.x) > fabsf(v1.y)) {
		float invLen = 1.f / sqrtf(v1.x*v1.x + v1.z*v1.z);
		*v2 = Vector(-v1.z * invLen, 0.f, v1.x * invLen);
	}
	else {
		float invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
		*v2 = Vector(0.f, v1.z * invLen, -v1.y * invLen);
	}
	*v3 = Cross(v1, *v2);
}
inline float Distance(const Point &p1, const Point &p2) {
	return (p1 - p2).Length();
}
inline float DistanceSquared(const Point &p1, const Point &p2) {
	return (p1 - p2).LengthSquared();
}
inline ostream &operator<<(ostream &os, const Point &v) {
	os << v.x << ", " << v.y << ", " << v.z;
	return os;
}
inline Point operator*(float f, const Point &p) {
	return p*f;
}
inline Normal operator*(float f, const Normal &n) {
	return Normal(f*n.x, f*n.y, f*n.z);
}
inline Normal Normalize(const Normal &n) {
	return n / n.Length();
}
inline Vector::Vector(const Normal &n)
  : x(n.x), y(n.y), z(n.z) { }
inline float Dot(const Normal &n1, const Vector &v2) {
	return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}
inline float Dot(const Vector &v1, const Normal &n2) {
	return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}
inline float Dot(const Normal &n1, const Normal &n2) {
	return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}
inline float AbsDot(const Normal &n1, const Vector &v2) {
	return fabsf(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}
inline float AbsDot(const Vector &v1, const Normal &n2) {
	return fabsf(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}
inline float AbsDot(const Normal &n1, const Normal &n2) {
	return fabsf(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}
inline ostream &operator<<(ostream &os, const Normal &v) {
	os << v.x << ", " << v.y << ", " << v.z;
	return os;
}
inline ostream &operator<<(ostream &os, Ray &r) {
	os << "org: " << r.o << "dir: " << r.d << " range [" <<
		r.mint << "," << r.maxt << "] time = " <<
		r.time;
	return os;
}
inline ostream &operator<<(ostream &os, const BBox &b) {
	os << b.pMin << " -> " << b.pMax;
	return os;
}
inline Vector SphericalDirection(float sintheta,
                              float costheta, float phi) {
	return Vector(sintheta * cosf(phi),
	              sintheta * sinf(phi),
				  costheta);
}
inline Vector SphericalDirection(float sintheta,
                                   float costheta,
								   float phi,
								   const Vector &x,
								   const Vector &y,
								   const Vector &z) {
	return sintheta * cosf(phi) * x +
		sintheta * sinf(phi) * y + costheta * z;
}
inline float SphericalTheta(const Vector &v) {
	return acosf(Clamp(v.z, -1.f, 1.f));
}
inline float SphericalPhi(const Vector &v) {
	float p = atan2f(v.y, v.x);
	return (p < 0.f) ? p + 2.f*M_PI : p;
}

#define perp(u,v)  ((u).x * (v).y - (u).y * (v).x)  // perp product  (2D)
// inSegment(): determine if a point is inside a segment
//    Input:  a point P, and a collinear segment S
//    Return: 1 = P is inside S
//            0 = P is  not inside S
inline int
inSegment(Point P, Segment S)
{
	if (S.P0.x != S.P1.x) {    // S is not  vertical
		if (S.P0.x <= P.x && P.x <= S.P1.x)
			return 1;
		if (S.P0.x >= P.x && P.x >= S.P1.x)
			return 1;
	}
	else {    // S is vertical, so test y  coordinate
		if (S.P0.y <= P.y && P.y <= S.P1.y)
			return 1;
		if (S.P0.y >= P.y && P.y >= S.P1.y)
			return 1;
	}
	return 0;
}

// intersect2D_2Segments(): find the 2D intersection of 2 finite segments
//    Input:  two finite segments S1 and S2
//    Output: *I0 = intersect point (when it exists)
//            *I1 =  endpoint of intersect segment [I0,I1] (when it exists)
//    Return: 0=disjoint (no intersect)
//            1=intersect  in unique point I0
//            2=overlap  in segment from I0 to I1
//            3=intersect, but outside segments S1 & S2
inline int intersect2D_2Segments(Segment S1, Segment S2, Point& I0, Point& I1)
{
	Vector    u = S1.P1 - S1.P0;
	Vector    v = S2.P1 - S2.P0;
	Vector    w = S1.P0 - S2.P0;

	//float     D = perp(u, v);
	//float     SMALL_NUM = 1e-8f;

	double D = (double(u.x) * double(v.y) - double(u.y) * double(v.x));
	double SMALL_NUM = 1.0e-8;
	// test if  they are parallel (includes either being a point)
	if (fabs(D) < SMALL_NUM) {           // S1 and S2 are parallel
		if (perp(u, w) != 0 || perp(v, w) != 0)  {
			return 0;                    // they are NOT collinear
		}
		// they are collinear or degenerate
		// check if they are degenerate  points
		float du = Dot(u, u);
		float dv = Dot(v, v);
		if (du == 0 && dv == 0) {            // both segments are points
			if (S1.P0.x != S2.P0.x && S1.P0.y != S2.P0.y)         // they are distinct  points
				return 0;
			I0 = S1.P0;                 // they are the same point
			return 1;
			//return 3;
		}
		if (du == 0) {                     // S1 is a single point
			if (inSegment(S1.P0, S2) == 0)  // but is not in S2
				return 0;
			I0 = S1.P0;
			return 1;
			//return 3;
		}
		if (dv == 0) {                     // S2 a single point
			if (inSegment(S2.P0, S1) == 0)  // but is not in S1
				return 0;
			I0 = S2.P0;
			return 1;
			//return 3;
		}
		// they are collinear segments - get  overlap (or not)
		float t0, t1;                    // endpoints of S1 in eqn for S2
		Vector w2 = S1.P1 - S2.P0;
		if (v.x != 0) {
			t0 = w.x / v.x;
			t1 = w2.x / v.x;
		}
		else {
			t0 = w.y / v.y;
			t1 = w2.y / v.y;
		}
		if (t0 > t1) {                   // must have t0 smaller than t1
			float t = t0; t0 = t1; t1 = t;    // swap if not
		}
		if (t0 > 1 || t1 < 0) {
			return 0;      // NO overlap
		}
		t0 = t0<0 ? 0 : t0;               // clip to min 0
		t1 = t1>1 ? 1 : t1;               // clip to max 1
		if (t0 == t1) {                  // intersect is a point
			I0 = S2.P0 + t0 * v;
			return 1;
			//return 3;
		}

		// they overlap in a valid subsegment
		I0 = S2.P0 + t0 * v;
		I1 = S2.P0 + t1 * v;
		return 2;
	}
	// the segments are skew and may intersect in a point
	// get the intersect parameter for S1
	float     sI = perp(v, w) / D;
	I0 = S1.P0 + sI * u;                // compute S1 intersect point
	if (sI < 0 || sI > 1)                // no intersect with S1
		return 3; // intersection outside the segments

	// get the intersect parameter for S2
	float     tI = perp(u, w) / D;
	if (tI < 0 || tI > 1)                // no intersect with S2
		return 3; // intersection outside the segments

	return 1;  // return for degenerated case
}


inline bool lineSegRectHit(FLOATVECTOR2 pStart, FLOATVECTOR2 pEnd, FLOATVECTOR2 bottom_left, FLOATVECTOR2 top_right)
{
	Segment L;
	L.P0 = Point(pStart.x, pStart.y);
	L.P1 = Point(pEnd.x, pEnd.y);

	Point rect[4];
	rect[0] = Point(bottom_left.x, bottom_left.y);
	rect[1] = Point(bottom_left.x, top_right.y);
	rect[2] = Point(top_right.x, top_right.y);
	rect[3] = Point(top_right.x, bottom_left.y);

	Point s0, s1;
	Segment rectL[4]; // assign line segments of the rectangle clockwise
	for (int i = 0; i < 4; i++)
	{
		rectL[i].P0 = rect[i%4];
		rectL[i].P1 = rect[(i + 1) % 4];
		int intersect = intersect2D_2Segments(L, rectL[i], s0, s1);
		if (intersect == 1) // if the test segment intersects with ANY of the rectangle edges, return true
			return true;
	}
	return false; // no intersections!
}

inline bool lineSegPolyHit(FLOATVECTOR2 pStart, FLOATVECTOR2 pEnd, const std::vector<FLOATVECTOR2>& poly)
{
	Segment L;
	L.P0 = Point(pStart.x, pStart.y);
	L.P1 = Point(pEnd.x, pEnd.y);

	size_t num_verts = poly.size();
	std::vector<Point> ppoly(num_verts);
	for (size_t i = 0; i < num_verts; i++)
		ppoly[i] = Point(poly[i].x, poly[i].y);

	Point s0, s1;
	Segment polyL; // assign line segments of the rectangle clockwise
	for (int i = 0; i < poly.size(); i++)
	{
		polyL.P0 = ppoly[i % num_verts];
		polyL.P1 = ppoly[(i + 1) % num_verts];
		int intersect = intersect2D_2Segments(L, polyL, s0, s1);
		if (intersect == 1) // if the test segment intersects with ANY of the rectangle edges, return true
			return true;
	}
	return false; // no intersections!
}

#endif // PBRT_GEOMETRY_H
