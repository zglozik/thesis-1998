
//
//	member functions for approximation version 1
//

#include	"approximation.hh"
#include	"householder.hh"
#include	<vtkMath.h>

//----------------------------------------------------------------------
//	approximate1( uK, vK, un, vn, surface, default_parameters, excl )
//----------------------------------------------------------------------

double Approximation::approximate1( uint uK, uint vK, uint un, uint vn,
				   BSplineSurface &surface,
				   bool default_parameters,
				   const Exclusion *excl )
{
	if( default_parameters ) init_parameters( measured.n(), measured.m() );

	KnotVector uk( uK+un+1 );
	KnotVector vk( vK+vn+1 );

	fill_knotvector( uk, un, 0, measured.n()-1 );
	fill_knotvector( vk, vn, 0, measured.m()-1 );

	surface.uset_param( un, uk );
	surface.vset_param( vn, vk );
	
	uint n = excl ? excl->size() : measured.n() * measured.m();
	uint m = uK * vK;

	Matrix<double> A( n, m+1 );
	Matrix<double> x(m);
	Matrix<Point> controls(uK, vK);

	fill_matrix1( A, surface.ubase(), surface.vbase(), excl );
	for( uint index = 0; index < 3; index++ ) {
		fill_b1( A, index, excl );
		solve_problem( A, x );
		fill_controls( controls, x, uK, vK, index );
	}
	surface.set_points( controls );

	return count_diff( surface, excl );
}

//----------------------------------------------------------------------
//	fill_matrix1( Matrix<double> &A, const Base &ubase, const Base &vbase,
//		const Exclusion *excl ) const
//----------------------------------------------------------------------

void Approximation::fill_matrix1( Matrix<double> &A, 
				 const Base &ubase, const Base &vbase,
				 const Exclusion *excl ) const
{
	uint n = measured.n() * measured.m();
	uint m = ubase.K() * vbase.K();
	uint k = 0;

	for( uint i = 0; i < n; i++ ) {
		uint u0 = i / measured.m();
		uint v0 = i % measured.m();
		if( !excl || excl->is_inside(u0,v0) ) {
			Par p = parameters(u0,v0);
			for( uint j = 0; j < m; j++ ) {
				uint ui = j / vbase.K();
				uint vi = j % vbase.K();
				A(k,j) = ubase(ui,p.u) * vbase(vi,p.v);
			}
			k++;
		}
	}
}

//----------------------------------------------------------------------
//	fill_b1( Matrix<double> &A, uint index, const Exclusion *excl ) const
//----------------------------------------------------------------------

void Approximation::fill_b1( Matrix<double> &A, uint index,
			   const Exclusion *excl ) const
{
	uint k = 0;
	uint l = A.m()-1;

	for( uint i = 0; i < measured.n(); i++ ) {
		for( uint j = 0; j < measured.m(); j++ ) {
			if( !excl || excl->is_inside(i,j) ) {
				A(k++,l) = measured(i,j)[index];
			}
		}
	}
}

