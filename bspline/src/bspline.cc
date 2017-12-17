
#include	"bspline.hh"
#include	<funcs.hh>

//
//	Base class
//

//----------------------------------------------------------------------
//	Base::Base( uint n, const KnotVector &kv )
//----------------------------------------------------------------------

Base::Base( uint n, const KnotVector &kv ) : kv(kv), n(n)
{
	TEST_EXPR( (uint) kv.size() >= 2*n+2 );
	update_m( kv );
}

//----------------------------------------------------------------------
//	Base::Base( uint n, uint k, double t0, double t1 )
//----------------------------------------------------------------------

Base::Base( uint n, uint k, double t0, double t1 ) : kv(k), n(n)
{
	TEST_EXPR( k >= 2*n+2 );
	TEST_EXPR( t1 > t0 );

	for( uint i = 0; i < k; i++ )
		kv[i] = t0 + (t1-t0) * i / (k-1);
	update_m( kv );
}

//----------------------------------------------------------------------
//	double Base::operator () ( uint i, double t ) const
//----------------------------------------------------------------------

double Base::operator () ( uint i, double t ) const
{
	TEST_EXPR( t >= kv.front() && t <= kv.back() );
	TEST_EXPR( i <= kv.size() - 1 - (n+1) );

	return compute( i, n, t );
}

//----------------------------------------------------------------------
//	double Base::derivate( uint i, double t, uint k, bool &error ) const
//----------------------------------------------------------------------

double Base::derivate( uint i, double t, uint k, bool &error ) const
{
	TEST_EXPR( t >= kv.front() && t <= kv.back() );
	TEST_EXPR( i <= kv.size() - 1 - (n+1) );

	double rv;
	Multis::iterator multi = m.find( t );
	error = false;

	if( multi != m.end() && k != 0 && k + *multi > n ) {
		error = true;
		rv = 0;
	} else if( k <= n ) {
		rv = compute_derivative( i, n, t, k );
	} else {
		rv = 0;
	}

	return rv;
}

//----------------------------------------------------------------------
//	double Base::min_t() const
//----------------------------------------------------------------------

double Base::min_t() const
{
	return kv[n];
}

//----------------------------------------------------------------------
//	double Base::max_t() const
//----------------------------------------------------------------------

double Base::max_t() const
{
	return kv[kv.size()-1-n];
}

//----------------------------------------------------------------------
//	double Base::compute( uint i, uint n, double t )
//----------------------------------------------------------------------

double Base::compute( uint i, uint n, double t ) const
{
	double	value;

	if( t < kv[i] || t > kv[i+n+1] ) {
		value = 0;
	} else if( n == 0 ) {
		if( t >= kv[i] &&  (t < kv[i+1] || 
		   (equal(t,kv[i+1]) && equal(t,kv.back()))) ) {
			value = 1;
		} else {
			value = 0;
		}
	} else {
		value = 0;
		if( !equal( kv[i+n], kv[i] ) ) {
			value += (t-kv[i]) / (kv[i+n]-kv[i]) * 
				compute(i,n-1,t);
		}
		if( !equal( kv[i+n+1], kv[i+1] ) ) {
			value += (kv[i+n+1]-t) / (kv[i+n+1]-kv[i+1]) *
				compute(i+1,n-1,t);
		}
	}

	return value;
}

//----------------------------------------------------------------------
//	double Base::compute_derivative( uint i, uint n, 
//					double t, uint k ) const
//----------------------------------------------------------------------

double Base::compute_derivative( uint i, uint n, double t, uint k ) const
{
	double rv;

	if( t < kv[i] || t > kv[i+n+1] ) {
		rv = 0;
	} else if( k == 0 ) {
		rv = compute( i, n, t );
	} else {
		rv = 0;
		if( !equal( kv[i+n], kv[i] ) ) {
			rv = compute_derivative( i, n-1, t, k-1 ) /
				(kv[i+n] - kv[i]);
		}
		if( !equal( kv[i+n+1], kv[i+1] ) ) {
			rv -= compute_derivative( i+1, n-1, t, k-1 ) /
				(kv[i+n+1] - kv[i+1]);
		}
		rv *= n;
	}

	return rv;
}

//----------------------------------------------------------------------
//	void Base::update_m( const KnotVector &kv )
//----------------------------------------------------------------------

void Base::update_m( const KnotVector &kv )
{
	empty( m );

	for( KnotVector::iterator knot = kv.begin(); knot != kv.end(); ) {
		uint	multi = 1;
		double	value = *knot++;
		
		while( knot != kv.end() && equal(value, *knot) ) {
			multi++;
			knot++;
		}
		m[ value ] = multi;
	}
}

//----------------------------------------------------------------------
//	uint Base::K() const
//----------------------------------------------------------------------

uint Base::K() const
{
	return kv.size()-1-n;
}

//
//	BSplineCurve class
//

//----------------------------------------------------------------------
//	BSplineCurve::BSplineCurve( const BSplineCurve &b )
//----------------------------------------------------------------------

BSplineCurve::BSplineCurve( const BSplineCurve &b ) : base(0), points(0)
{
	if( b.base ) base = new Base( *b.base );
	if( b.points ) points = new Points( *b.points );
}

//----------------------------------------------------------------------
//	BSplineCurve::~BSplineCurve()
//----------------------------------------------------------------------

BSplineCurve::~BSplineCurve()
{
	delete base;
	delete points;
}

//----------------------------------------------------------------------
//	void BSplineCurve::set_param( uint n, const KnotVector &kv )
//----------------------------------------------------------------------

void BSplineCurve::set_param( uint n, const KnotVector &kv )
{
	delete base;
	base = new Base( n, kv );
}

//----------------------------------------------------------------------
//	void BSplineCurve::set_param( uint n, uint k, double t0, double t1 )
//----------------------------------------------------------------------

void BSplineCurve::set_param( uint n, uint k, double t0, double t1 )
{
	delete base;
	base = new Base( n, k, t0, t1 );
}

//----------------------------------------------------------------------
//	void BSplineCurve::set_points( const Points &p )
//----------------------------------------------------------------------

void BSplineCurve::set_points( const Points &p )
{
	delete points;
	points = new Points( p );
}

//----------------------------------------------------------------------
//	Point BSplineCurve::operator () ( double t ) const
//----------------------------------------------------------------------

Point BSplineCurve::operator () ( double t ) const
{
	TEST_EXPR( base != NULL && points != NULL );
	TEST_EXPR( K() == (uint) points->size() );
	TEST_EXPR( min_t() <= t && t <= max_t() );

	Point	result = 0;
	int	i = 0;
	for( Points::iterator p = points->begin(); p != points->end(); p++ ) {
		result += (*p) * (*base)(i++, t);
	}
	return result;
}

//----------------------------------------------------------------------
//	Point BSplineCurve::derivate( double t, uint k, bool &error ) const
//----------------------------------------------------------------------

Point BSplineCurve::derivate( double t, uint k, bool &error ) const
{
	TEST_EXPR( base != NULL && points != NULL );
	TEST_EXPR( K() == (uint) points->size() );
	TEST_EXPR( min_t() <= t && t <= max_t() );

	Point rv = 0;
	uint i = 0;
	error = false;

	for( Points::iterator p = points->begin(); 
	    p != points->end() && !error; p++ ) {
		rv += (*p) * base->derivate(i++, t, k, error);
	}

	return rv;
}

//
//	BSplineSurface class
//

//----------------------------------------------------------------------
//	BSplineSurface::BSplineSurface( const BSplineSurface &b )
//----------------------------------------------------------------------

BSplineSurface::BSplineSurface( const BSplineSurface &b ) 
: _ubase(0), _vbase(0), points2(0)
{
	if( b._ubase ) _ubase = new Base( *b._ubase );
	if( b._vbase ) _vbase = new Base( *b._vbase );
	if( b.points2 ) points2 = new Matrix<Point>( *b.points2 );
}

//----------------------------------------------------------------------
//	BSplineSurface::~BSplineSurface()
//----------------------------------------------------------------------

BSplineSurface::~BSplineSurface()
{
	delete _ubase;
	delete _vbase;
	delete points2;
}

//----------------------------------------------------------------------
//	void BSplineSurface::uset_param( uint n, const KnotVector &kv )
//----------------------------------------------------------------------

void BSplineSurface::uset_param( uint n, const KnotVector &kv )
{
	delete _ubase;
	_ubase = new Base(n, kv);
}

//----------------------------------------------------------------------
//	void BSplineSurface::uset_param( uint n, uint k, double t0, double t1 )
//----------------------------------------------------------------------

void BSplineSurface::uset_param( uint n, uint k, double t0, double t1 )
{
	delete _ubase;
	_vbase = new Base(n, k, t0, t1);
}

//----------------------------------------------------------------------
//	void BSplineSurface::vset_param( uint n, const KnotVector &kv )
//----------------------------------------------------------------------

void BSplineSurface::vset_param( uint n, const KnotVector &kv )
{
	delete _vbase;
	_vbase = new Base(n, kv);
}

//----------------------------------------------------------------------
//	void BSplineSurface::vset_param( uint n, uint k, double t0, double t1 )
//----------------------------------------------------------------------

void BSplineSurface::vset_param( uint n, uint k, double t0, double t1 )
{
	delete _vbase;
	_vbase = new Base(n, k, t0, t1);
}

//----------------------------------------------------------------------
//	void BSplineSurface::set_points( const Matrix<Point> &p2 )
//----------------------------------------------------------------------

void BSplineSurface::set_points( const Matrix<Point> &p2 )
{
	delete points2;
	points2 = new Matrix<Point>( p2 );
}

//----------------------------------------------------------------------
//	Point BSplineSurface::operator () ( double u, double v ) const
//----------------------------------------------------------------------

Point BSplineSurface::operator () ( double u, double v ) const
{
	TEST_EXPR( _ubase != NULL && _vbase != NULL && points2 != NULL );
	TEST_EXPR( uK() == points2->n() && vK() == points2->m() );

	uint unum = points2->n();
	uint vnum = points2->m();

	Point res_uv = 0;
	for( uint i = 0; i < unum; i++ ) {
		Point	res_v = 0;
		for( uint j = 0; j < vnum; j++ ) {
			res_v += (*points2)(i, j) * (*_vbase)(j, v);
		}
		res_uv += res_v * (*_ubase)(i, u);
	}

	return res_uv;
}

//----------------------------------------------------------------------
//	Point BSplineSurface::derivate( double u, double v, 
//				       uint uk, uint vk, bool &error ) const
//----------------------------------------------------------------------

Point BSplineSurface::derivate( double u, double v, 
			       uint uk, uint vk, bool &error ) const
{
	TEST_EXPR( _ubase != NULL && _vbase != NULL && points2 != NULL );
	TEST_EXPR( uK() == points2->n() && vK() == points2->m() );

	uint unum = points2->n();
	uint vnum = points2->m();
	error = false;

	Point res_uv = 0;
	Point res_v;
	for( uint i = 0; i < unum && !error; i++ ) {
		res_v = 0;
		for( uint j = 0; j < vnum && !error; j++ ) {
			res_v += (*points2)(i, j) * 
				_vbase->derivate(j, v, vk, error);
		}
		if( !error ) {
			res_uv += res_v * _ubase->derivate(i, u, uk, error);
		}
	}

	return res_uv;
}

//----------------------------------------------------------------------
//	bool BSplineSurface::first_q( double u, double v, 
//				     double &E, double &F, double &G ) const
//----------------------------------------------------------------------

bool BSplineSurface::first_q( double u, double v, 
			     double &E, double &F, double &G ) const
{
	bool error = false;
	Point r_u, r_v;
	r_u = derivate( u, v, 1, 0, error );
	if( !error ) r_v = derivate( u, v, 0, 1, error );
	if( !error ) {
		E = r_u * r_u;
		F = r_u * r_v;
		G = r_v * r_v;
	}

	return error;
}

//----------------------------------------------------------------------
//	bool BSplineSurface::second_q( double u, double v, 
//				      double &L, double &M, double &N ) const
//----------------------------------------------------------------------

bool BSplineSurface::second_q( double u, double v, 
			      double &L, double &M, double &N ) const
{
	bool error = false;
	Point r_uu, r_uv, r_vv, n;
	r_uu = derivate( u, v, 2, 0, error );
	if( !error ) r_uv = derivate( u, v, 1, 1, error );
	if( !error ) r_vv = derivate( u, v, 0, 2, error );
	if( !error ) n = normal( u, v, error );

	if( !error ) {
		L = r_uu * n;
		M = r_uv * n;
		N = r_vv * n;
	}

	return error;
}

//----------------------------------------------------------------------
//	Point BSplineSurface::normal( double u, double v, bool &error ) const
//----------------------------------------------------------------------

Point BSplineSurface::normal( double u, double v, bool &error ) const
{
	Point r_u, r_v;
	r_u = derivate( u, v, 1, 0, error );
	if( !error ) r_v = derivate( u, v, 0, 1, error );

	return error ? Point(0) : ::normal( vector_mul(r_u, r_v) );
}

//----------------------------------------------------------------------
//	double BSplineSurface::curvature( double u, double v, 
//				double du, double dv,  bool &error ) const
//----------------------------------------------------------------------

double BSplineSurface::curvature( double u, double v, double du, double dv,
				 bool &error ) const
{
	double E, F, G, L, M, N;
	double rv = 0;

	error = first_q( u, v, E, F, G );
	error = error || second_q( u, v, L, M, N );
	if( !error ) {
		double denominator = E*du*du + 2*F*du*dv + G*dv*dv;
		if( !(error = equal(denominator, 0)) ) {
			rv = (L*du*du + 2*M*du*dv + N*dv*dv) / denominator;
		}
	}

	return rv;
}

//----------------------------------------------------------------------
//	double BSplineSurface::gaussian( double u, double v, 
//					bool &error ) const
//----------------------------------------------------------------------

double BSplineSurface::gaussian( double u, double v, bool &error ) const
{
	double E, F, G, L, M, N;
	double rv = 0;
	error = first_q( u, v, E, F, G );
	error = error || second_q( u, v, L, M, N );
	if( !error ) {
		double denominator = E*G - F*F;
		if( !(error = equal(denominator, 0)) ) {
			rv = (L*N - M*M) / denominator;
		}
	}

	return rv;
}

//----------------------------------------------------------------------
//	double BSplineSurface::minkowski( double u, double v, 
//					bool &error ) const
//----------------------------------------------------------------------

double BSplineSurface::minkowski( double u, double v, bool &error ) const
{
	double E, F, G, N, M, L;
	double rv = 0;

	error = first_q( u, v, E, F, G );
	error = error || second_q( u, v, L, M, N );
	if( !error ) {
		double denominator = E*G - F*F;
		if( !(error = equal(denominator, 0)) ) {
			rv = (E*N - 2*F*M + G*L) / denominator;
		}
	}

	return rv;
}

//
//	Streams
//

//----------------------------------------------------------------------
//	ostream &operator << ( ostream &o, const Point &b )
//----------------------------------------------------------------------

ostream &operator << ( ostream &o, const Point &b )
{
	o << b[0] << ' ';
	o << b[1] << ' ';
	o << b[2] << ' ';
	return o;
}

//----------------------------------------------------------------------
//	istream &operator >> ( istream &i, Point &b )
//----------------------------------------------------------------------

istream &operator >> ( istream &i, Point &b )
{
	i >> b[0];
	i >> b[1];
	i >> b[2];
	return i;
}

//----------------------------------------------------------------------
//	ostream &operator << ( ostream &o, const Base &b )
//----------------------------------------------------------------------

ostream &operator << ( ostream &o, const Base &b )
{
	o << b.n << endl;
	o << b.kv;
	return o;
}

//----------------------------------------------------------------------
//	istream &operator >> ( istream &i, Base &b )
//----------------------------------------------------------------------

istream &operator >> ( istream &i, Base &b )
{
	i >> b.n;
	i >> b.kv;
	b.update_m( b.kv );
	return i;
}

//----------------------------------------------------------------------
//	ostream &operator << ( ostream &o, const BSplineSurface &b )
//----------------------------------------------------------------------

ostream &operator << ( ostream &o, const BSplineSurface &b )
{
	TEST_EXPR( b._ubase != NULL && b._vbase != NULL && b.points2 != NULL );

	o << *b._ubase;
	o << *b._vbase;
	o << *b.points2;
	o << endl;
	return o;
}

//----------------------------------------------------------------------
//	istream &operator >> ( istream &i, BSplineSurface &b )
//----------------------------------------------------------------------

istream &operator >> ( istream &i, BSplineSurface &b )
{
	delete b._ubase;
	delete b._vbase;
	delete b.points2;

	b._ubase = new Base;
	b._vbase = new Base;
	b.points2 = new Matrix<Point>;

	i >> *b._ubase;
	i >> *b._vbase;
	i >> *b.points2;
	return i;
}

//----------------------------------------------------------------------
/*
int main()
{
	KnotVector ukv(8);
	ukv[0] = -3;
	ukv[1] = -2;
	ukv[2] = -1;
	ukv[3] = 0;
	ukv[4] = 1;
	ukv[5] = 2;
	ukv[6] = 3;
	ukv[7] = 4;
	KnotVector vkv = ukv;

	Matrix<Point> points(4,4);

	points(0,0) = Point(-1, 0, -1);
	points(0,1) = Point(-1/3., 1, -1);
	points(0,2) = Point(1/3., 1, -1);
	points(0,3) = Point(1, 0, -1);

	points(1,0) = Point(-1, 1/3., -1/3.);
	points(1,1) = Point(-1/3., 1, -1/3.);
	points(1,2) = Point(1/3., 1, -1/3.);
	points(1,3) = Point(1, 1/3., -1/3.);

	points(2,0) = Point(-1, 1/3., 1/3.);
	points(2,1) = Point(-1/3., 1, 1/3.);
	points(2,2) = Point(1/3., 1, 1/3.);
	points(2,3) = Point(1, 1/3., 1/3.);

	points(3,0) = Point(-1, 0, 1);
	points(3,1) = Point(-1/3., 1, 1);
	points(3,2) = Point(1/3., 1, 1);
	points(3,3) = Point(1, 0, 1);

	cout << points << endl;

	BSplineSurface surface;
	surface.uset_param( 3, ukv );
	surface.vset_param( 3, vkv );
	surface.set_points( points );

	cout << surface( 0, 0 ) << endl; 
	cout << surface( 0.5, 0 ) << endl;

	bool error;
	cout << surface.normal( 0, 0, error ) << endl;
	cout << error << endl;
	cout << surface.curvature( 0, 0, 1, 1, error ) << endl;
	cout << error << endl;
	cout << surface.gaussian( 0, 0, error ) << endl;
	cout << error << endl;
	cout << surface.minkowski( 0, 0, error ) << endl;
	cout << error << endl;

	KnotVector kv(8);
	kv[0] = 0;
	kv[1] = 1;
	kv[2] = 2;
	kv[3] = 3;
	kv[4] = 4;
	kv[5] = 5;
	kv[6] = 6;
	kv[7] = 7;

	Base	base( 3, kv );

	cout << "------------\n";
	cout << base.derivate( 0, 2, 2, error ) << endl;
	cout << error << endl;
	cout << base.derivate( 0, 2.5, 3, error ) << endl;
	cout << error << endl;

	return 0;
}
*/
