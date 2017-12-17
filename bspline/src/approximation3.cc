
//
//	member functions for approximation version 3
//

#include	"approximation.hh"
#include	"householder.hh"
#include	<vtkMath.h>

//----------------------------------------------------------------------
//	approximate3( uint uK, uint vK, uint un, uint vn,
//		    BSplineSurface &surface,
//		    const Exclusion *excl = 0 );
//----------------------------------------------------------------------

double Approximation::approximate3( uint uK, uint vK, uint un, uint vn,
				   BSplineSurface &surface, comp_fun act_fun,
				   const Exclusion *excl = 0 )
{
	double square_e;
	uint steps = 0;

	cout << "Step " << ++steps << ":\n";
	cout << "least squares method..."; cout.flush();
	square_e = approximate2( uK, vK, un, vn, surface, true, excl );
	cout << "\naverage of squares of differences: " << square_e << endl;
	
	cout << "parameter correction..."; cout.flush();
	correct_parameters( surface, excl, act_fun );
	cout << endl;
	while( steps < 6 ) {
		cout << "Step " << ++steps << ":\n";
		cout << "least squares method..."; cout.flush();
		square_e = approximate2(uK, vK, un, vn, surface, false, excl);
		cout << "\naverage of squares of differences: ";
		cout << square_e <<endl;

		cout << "parameter correction..."; cout.flush();
		correct_parameters( surface, excl, act_fun );
		cout << endl;
	}
	return count_diff( surface, excl );
}

//----------------------------------------------------------------------
//	compute_length( const BSplineSurface &surface,
//			uint i, uint j, double &lu, double &lv ) const
//----------------------------------------------------------------------

void Approximation::compute_length( const BSplineSurface &surface,
				   uint i, uint j, 
				   double &lu, double &lv ) const
{
	double u, v;
	Point cur, prev;
	u = parameters(i,j).u;
	lv = 0;
	prev = surface(u,0);
	for( uint vj = 1; vj < measured.m(); vj++ ) {
		cur = surface(u,vj);
		lv += abs( cur - prev );
		prev = cur;
	}
	v = parameters(i,j).v;
	lu = 0;
	prev = surface(0,v);
	for( uint ui = 1; ui < measured.n(); ui++ ) {
		cur = surface(ui,v);
		lu += abs( cur - prev );
		prev = cur;
	}
}

//----------------------------------------------------------------------
//	correct_parameters( const BSplineSurface &surface,
//			    const Exclusion *excl )
//----------------------------------------------------------------------

void Approximation::correct_parameters( const BSplineSurface &surface,
				       const Exclusion *excl, 
				       comp_fun act_fun )
{
	for( uint i = 0; i < measured.n(); i++ ) {
		for( uint j = 0; j < measured.m(); j++ ) {
			if( !excl || excl->is_inside(i,j) ) {
				(this->*act_fun)( surface, i, j );
			}
		}
	}
}

//----------------------------------------------------------------------
//	correct_parameter( const BSplineSurface &surface, uint i, uint j )
//----------------------------------------------------------------------

void Approximation::correct_parameter( const BSplineSurface &surface,
				      uint i, uint j )
{
	bool	error;
	double	lu, lv;
	Point	r_u, r_v;
	Par	p = parameters(i,j);

	compute_length( surface, i, j, lu, lv );
	r_u = surface.derivate( p.u, p.v, 1, 0, error );
	if( !error ) r_v = surface.derivate( p.u, p.v, 0, 1, error );
	if( error ) return;

	r_u = normal(r_u);
	r_v = normal(r_v);
	Point X = surface(p.u,p.v);
	Point M = measured(i,j);
	Point diff = M - X;
	double dus = diff * r_u;
	double dvs = diff * r_v;
	
	parameters(i,j).u = min( max(p.u + dus/lu * (measured.n()-1),
				     surface.umin_t()),
				surface.umax_t() );
	parameters(i,j).v = min( max(p.v + dvs/lv * (measured.m()-1),
				     surface.vmin_t()),
				surface.vmax_t() );
}

//----------------------------------------------------------------------
//	correct_parameter2( const BSplineSurface &surface, uint i, uint j )
//----------------------------------------------------------------------

void Approximation::correct_parameter2( const BSplineSurface &surface,
				       uint i, uint j )
{
	bool	error;
	Point	r_u, r_v;
	Par	p = parameters(i,j);

	r_u = surface.derivate( p.u, p.v, 1, 0, error );
	if( !error ) r_v = surface.derivate( p.u, p.v, 0, 1, error );
	if( error ) return;

	Point Ps = measured( i, j );
	Point X = surface( p.u, p.v );
	Point Ds = Ps - X;
	double F = (r_u*r_u)*(r_v*r_v) - (r_u*r_v)*(r_u*r_v);

	if( !equal(F, 0) ) {
		double dus = ((r_v*r_v)*(Ds*r_u) - (r_u*r_v)*(Ds*r_v)) / F;
		double dvs = ((r_u*r_u)*(Ds*r_v) - (r_u*r_v)*(Ds*r_u)) / F;

		parameters(i,j).u = min( max(p.u + dus, surface.umin_t()),
					surface.umax_t() );
		parameters(i,j).v = min( max(p.v + dvs, surface.vmin_t()),
					surface.vmax_t() );
	}
}

//----------------------------------------------------------------------
//	correct_parameter3( const BSplineSurface &surface, uint i, uint j )
//----------------------------------------------------------------------

void Approximation::correct_parameter3( const BSplineSurface &surface,
				       uint i, uint j )
{
	bool	error;
	Point	r_u, r_v;
	Point	r_uu, r_vv;
	Par	p = parameters(i,j);

	r_u = surface.derivate( p.u, p.v, 1, 0, error );
	if( !error ) r_v = surface.derivate( p.u, p.v, 0, 1, error );
	if( !error ) r_uu = surface.derivate( p.u, p.v, 2, 0, error );
	if( !error ) r_vv = surface.derivate( p.u, p.v, 0, 2, error );
	if( error ) return;
	

	Point Ps = measured( i, j );
	Point X = surface( p.u, p.v );
	Point Ds = Ps - X;

	double dus = -1*(Ds*r_u) / (Ds*r_uu - r_u*r_u);
	double dvs = -1*(Ds*r_v) / (Ds*r_vv - r_v*r_v);

	parameters(i,j).u = min( max(p.u + dus, surface.umin_t()),
				surface.umax_t() );
	parameters(i,j).v = min( max(p.v + dvs, surface.vmin_t()),
				surface.vmax_t() );
}

