
#ifndef	_BSPLINE_HH_
#define	_BSPLINE_HH_

#include	<stdtypes.h>
#include	<array.hh>
#include	<list.hh>
#include	<vector.hh>
#include	<angles.hh>
#include	<map.hh>
#include	<matrix.hh>

typedef Array<double>	KnotVector;
typedef	Vector<double>	Point;
typedef	Array<Point>	Points;

class BSplineCurve;
class BSplineSurface;

class Base {

friend	ostream &operator << ( ostream &o, const Base &b );
friend	istream &operator >> ( istream &i, Base &b );

    public:
	Base() : kv(0), n(0) { }
	Base( uint n, const KnotVector &kv );
	Base( uint n, uint k, double t0 = 0.0, double t1 = 1.0 );

	double	operator () ( uint i, double t ) const;
	double	derivate( uint i, double t, uint k, bool &error ) const;

	double	min_t() const;	// minimum parameter value to evaluate
	double	max_t() const;	// maximum parameter value to evaluate
	uint	K() const;	// number of base functions

    private:
	double	compute( uint n, uint i, double t ) const;
	double	compute_derivative( uint i, uint n, double t, uint k ) const;

	void	update_m( const KnotVector &kv );

	typedef	Map<double,uint> Multis;

	KnotVector	kv;	// knot parameter values
	Multis		m;	// multiplicity of knots
	uint		n;	// degree of base function
};

class BSplineCurve {

    public:
	BSplineCurve() : base(0), points(0) { }
	BSplineCurve( const BSplineCurve & b );
virtual	~BSplineCurve();

	void	set_param( uint n, const KnotVector &kv );
	void	set_param( uint n, uint k, double t0 = 0.0, double t1 = 1.0 );
	void	set_points( const Points &p );

	Points	&get_points() const { return *points; }

	Point	operator () ( double t ) const;
	Point	derivate( double t, uint k, bool &error ) const;

	double	min_t() const { return base->min_t(); }
	double	max_t() const { return base->max_t(); }
	uint	K() const { return base->K(); }

    private:
	Base	*base;
	Points	*points;
};

class BSplineSurface {

friend	ostream &operator << ( ostream &o, const BSplineSurface &b );
friend	istream &operator >> ( istream &i, BSplineSurface &b );

    public:
	BSplineSurface() : _ubase(0), _vbase(0), points2(0) { }
	BSplineSurface( const BSplineSurface &b );
virtual	~BSplineSurface();

	void	uset_param( uint n, const KnotVector &kv );
	void	uset_param( uint n, uint k, double t0 = 0.0, double t1 = 1.0 );
	void	vset_param( uint n, const KnotVector &kv );
	void	vset_param( uint n, uint k, double t0 = 0.0, double t1 = 1.0 );
	void	set_points( const Matrix<Point> &p2 );

	Matrix<Point> &get_points() const { return *points2; }

	Point 	operator () ( double u, double v ) const;
	Point	derivate( double u, double v, 
			 uint uk, uint vk, bool &error ) const;

	// first base quantities
	bool	first_q( double u, double v, 
		      double &E, double &F, double &G ) const;

	// second base quantities
	bool	second_q( double u, double v, 
			 double &L, double &M, double &N ) const;

	Point	normal( double u, double v, bool &error ) const;
	double	curvature( double u, double v, double du, double dv,
		  bool &error ) const;
	double	gaussian( double u, double v, bool &error ) const;
	double	minkowski( double u, double v, bool &error ) const;

	const Base &ubase() const { return *_ubase; }
	double 	umin_t() const { return _ubase->min_t(); }
	double 	umax_t() const { return _ubase->max_t(); }
	const Base &vbase() const { return *_vbase; }
	double 	vmin_t() const { return _vbase->min_t(); }
	double 	vmax_t() const { return _vbase->max_t(); }
	uint 	uK() const { return _ubase->K(); }
	uint 	vK() const { return _vbase->K(); }

    private:
	Base	*_ubase;
	Base	*_vbase;
	Matrix<Point> *points2;
};

ostream &operator << ( ostream &o, const BSplineSurface &b );
istream &operator >> ( istream &i, BSplineSurface &b );

#endif // _BSPLINE_HH_
