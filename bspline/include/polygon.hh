
#ifndef	_POLYGON_HH_
#define	_POLYGON_HH_

#include	<debug.hh>
#include	<stdtypes.h>
#include	<array.hh>

class Polygon {

friend	ostream &operator << ( ostream &out, const Polygon &polygon );
friend	istream &operator >> ( istream &in, Polygon &polygon );

    public:
	struct Point {
		Point( double x = 0, double y = 0 ) : x(x), y(y) { }
		double	x;
		double	y;
	};

	typedef	Array<Point>	Points;

    public:
	Polygon( uint size );
	Polygon( double tx = 0, double ty = 0, double bx = 1, double by = 1,
		uint size = 0 );
	Polygon( const Polygon &p );

	Points	&points();
	const Points &points() const;

	bool	is_inside( double x, double y ) const;
	bool	is_outside( double x, double y ) const;

    private:
	bool	intersect( double a1, double a2, double b1, double b2,
			  double c1, double c2, double d1, double d2 ) const;
	bool	is_on( double x1, double x2, 
		      double a1, double a2, double b1, double b2 ) const;

	Points	_points;
	double	tx, ty;
	double	bx, by;
};

inline Polygon::Polygon( uint size ) 
: _points(size), tx(0), ty(0), bx(1), by(1)
{
}

inline Polygon::Polygon(double tx, double ty, double bx, double by, uint size)
: _points(size), tx(tx), ty(ty), bx(bx), by(by)
{
}

inline Polygon::Polygon( const Polygon &p )
: _points(p._points), tx(p.tx), ty(p.ty), bx(p.bx), by(p.by)
{
}

inline Polygon::Points &Polygon::points()
{
	return _points;
}

inline const Polygon::Points &Polygon::points() const
{
	return _points;
}

inline ostream &operator << ( ostream &out, const Polygon::Point &point )
{
	out << point.x << ' ' << point.y;
	return out;
}

inline istream &operator >> ( istream &in, Polygon::Point &point )
{
	in >> point.x;
	in >> point.y;
	return in;
}

inline ostream &operator << ( ostream &out, const Polygon &polygon )
{
	out << polygon.tx << ' ' << polygon.ty << ' ';
	out << polygon.bx << ' ' << polygon.by << endl;
	out << polygon.points();
	return out;
}

inline istream &operator >> ( istream &in, Polygon &polygon )
{
	in >> polygon.tx >> polygon.ty;
	in >> polygon.bx >> polygon.by;
	in >> polygon.points();
	return in;
}

#endif // _POLYGON_HH_
