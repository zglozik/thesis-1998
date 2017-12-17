
#include	"functional.hh"

//----------------------------------------------------------------------
//	compute_integrate( Matrix<double> &A, const Base &base,
//				uint d1, uint d2 ) const
//----------------------------------------------------------------------

void Functional::compute_integrate( Matrix<double> &A, const Base &base,
				   uint d1, uint d2 ) const
{
	const double	step = (base.max_t()-base.min_t()) / (base.K() * 10);
	bool		error;
	double		v1, v2;

	for( uint i = 0; i < base.K(); i++ ) {
		for( uint j = 0; j < base.K(); j++ ) {
			double pv1 = base.derivate(i,base.min_t(),d1,error);
			double pv2 = base.derivate(j,base.min_t(),d2,error);
	
			A( i, j ) = 0;
			for( double u = base.min_t() + step; u <= base.max_t();
			    u += step ) {
				v1 = base.derivate( i, u, d1, error );
				v2 = base.derivate( j, u, d2, error );
				A( i, j ) += (v1*v2 + pv1*pv2) * step / 2;
				pv1 = v1;
				pv2 = v2;
			}
		}
	}
}

//----------------------------------------------------------------------
//	integrate( const BSplineSurface &surface ) const
//----------------------------------------------------------------------

double Functional::integrate( const BSplineSurface &surface ) const
{
	double ustep = (surface.umax_t()-surface.umin_t()) / (surface.uK()*10);
	double vstep = (surface.vmax_t()-surface.vmin_t()) / (surface.vK()*10);
	double rv = 0;

	for( double u = surface.umin_t() + ustep/2.0; u < surface.umax_t();
	    u += ustep ) {
		for( double v = surface.vmin_t() + vstep/2.0;
		    v < surface.vmax_t(); v += vstep ) {
			rv += value(surface, u, v) * ustep * vstep;
		}
	}
	return rv;

}

//----------------------------------------------------------------------
//	class Q1, try to minimize the gradients of the surface
//----------------------------------------------------------------------

void Q1::init( const Base &ubase, const Base &vbase )
{
	_U11.set_size( ubase.K(), ubase.K() );
	_U00.set_size( ubase.K(), ubase.K() );
	_V11.set_size( vbase.K(), vbase.K() );
	_V00.set_size( vbase.K(), vbase.K() );

	compute_integrate( _U11, ubase, 1, 1 );
	compute_integrate( _U00, ubase, 0, 0 );
	compute_integrate( _V11, vbase, 1, 1 );
	compute_integrate( _V00, vbase, 0, 0 );
}

double Q1::derivate( uint k, uint l, uint i, uint j ) const
{
	return 2 * _U11(i,k) * _V00(j,l) + _U00(i,k) * _V11(j,l);
}

double Q1::value( const BSplineSurface &surface, double u, double v ) const
{
	bool	error;
	Point	du, dv;

	du = surface.derivate( u, v, 1, 0, error );
	if( !error ) {
		dv = surface.derivate( u, v, 0, 1, error );
	}
	return error ? 0 : (du*du + dv*dv);
}

//----------------------------------------------------------------------
//	class Q2, try to minimize the acceleration of the surface
//----------------------------------------------------------------------

void Q2::init( const Base &ubase, const Base &vbase )
{
	_U22.set_size( ubase.K(), ubase.K() );
	_U20.set_size( ubase.K(), ubase.K() );
	_U02.set_size( ubase.K(), ubase.K() );
	_U00.set_size( ubase.K(), ubase.K() );
	_V00.set_size( vbase.K(), vbase.K() );
	_V02.set_size( vbase.K(), vbase.K() );
	_V20.set_size( vbase.K(), vbase.K() );
	_V22.set_size( vbase.K(), vbase.K() );

	compute_integrate( _U22, ubase, 2, 2 );
	compute_integrate( _U20, ubase, 2, 0 );
	compute_integrate( _U02, ubase, 0, 2 );
	compute_integrate( _U00, ubase, 0, 0 );
	compute_integrate( _V00, vbase, 0, 0 );
	compute_integrate( _V02, vbase, 0, 2 );
	compute_integrate( _V20, vbase, 2, 0 );
	compute_integrate( _V22, vbase, 2, 2 );
}

double Q2::derivate( uint k, uint l, uint i, uint j ) const
{
	return 2 * _U22(i,k) * _V00(j,l) + _U20(i,k) * _V02(j,l) +
		_U02(i,k) * _V20(j,l) + _U00(i,k) * _V22(j,l);
}

double Q2::value( const BSplineSurface &surface, double u, double v ) const
{
	bool	error;
	Point	du, dv;
	double	rv;

	du = surface.derivate( u, v, 2, 0, error );
	if( !error ) {
		dv = surface.derivate( u, v, 0, 2, error );
	}
	if( error ) {
		rv = 0;
	} else {
		Point tmp = du + dv;
		rv = tmp*tmp;
	}
	return rv;
}

//----------------------------------------------------------------------
//	class Q3, try to minimize the distances between control points
//----------------------------------------------------------------------

void Q3::init( const Base &ubase, const Base &vbase )
{
	_uK = ubase.K();
	_vK = vbase.K();
}

double Q3::derivate( uint k, uint l, uint i, uint j ) const
{
	int	rv;

	if( k == i && l == j ) {
		rv = 0;
		if( k+1 < _uK ) rv += 4;
		if( k > 0 )     rv += 4;
		if( l+1 < _vK ) rv += 4;
		if( l > 0 )     rv += 4;
	} else if( (k == i+1 || k+1 == i) && j == l ||
		k == i && (l == j+1 || l+1 == j) ) {
		rv = -4;
	} else {
		rv = 0;
	}
	
	return rv;
}

double Q3::value( const BSplineSurface &surface, double u, double v ) const
{
	return 0;
}
