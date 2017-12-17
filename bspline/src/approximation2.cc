
//
//	member functions for approximation version 2
//

#include	"approximation.hh"
#include	"householder.hh"
#include	<vtkMath.h>

//----------------------------------------------------------------------
//	approximate2( uK, vK, un, vn, surface, default_parameters, excl ) 
//----------------------------------------------------------------------

double Approximation::approximate2( uint uK, uint vK, uint un, uint vn,
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
	
	uint n = uK * vK;

	Matrix<double> A( n, n+1 );
	Matrix<double> x( n );
	Matrix<Point> controls( uK, vK );
	Matrix<double> base_product( measured.n()*measured.m(), uK*vK );

	fill_base_product( base_product, surface.ubase(), surface.vbase() );
	fill_matrix2( A, base_product, surface.ubase(), surface.vbase(),
		     excl );
	for( uint index = 0; index < 3; index++ ) {
		fill_b2( A, base_product, index, 
			surface.ubase(), surface.vbase(), excl );
		solve_problem( A, x );
		fill_controls( controls, x, uK, vK, index );
	}
	surface.set_points( controls );

	return count_diff( surface, excl );
}

//----------------------------------------------------------------------
//	fill_matrix2( Matrix<double> &A, const Matrix<double> &base_product,
//		const Base &ubase, const Base &vbase,
//		const Exclusion *excl ) const
//----------------------------------------------------------------------

void Approximation::fill_matrix2( Matrix<double> &A, 
				 const Matrix<double> &base_product,
				 const Base &ubase, const Base &vbase,
				 const Exclusion *excl ) const
{
	for( uint row = 0; row < A.n(); row++ ) {
		uint i = row / vbase.K();
		uint j = row % vbase.K();
		for( uint col = 0; col < A.m()-1; col++ ) {
			uint m = col / vbase.K();
			uint n = col % vbase.K();
			A(row,col)=count_matrix_sum(base_product,ubase,vbase,
						    i,j,m,n,excl);
		}
	}
}

//----------------------------------------------------------------------
//	fill_b2( Matrix<double> &A, const Matrix<double> &bp, uint index,
//		const Base &ubase, const Base &vbase,
//		const Exclusion *excl ) const
//----------------------------------------------------------------------

void Approximation::fill_b2( Matrix<double> &A, const Matrix<double> &bp,
			    uint index, const Base &ubase, const Base &vbase,
			    const Exclusion *excl ) const
{
	uint num_points = measured.n()*measured.m();
	for( uint row = 0; row < A.n(); row++ ) {
		uint i = row / vbase.K();
		uint j = row % vbase.K();
		double sum = 0;
		for( uint points = 0; points < num_points; points++ ) {
			uint k = points / measured.m();
			uint l = points % measured.m();
			if( !excl || excl->is_inside(k,l) ) {
				sum += measured(k,l)[index] *
					base_product(bp,ubase,vbase,k,l,i,j);
			}
		}
		A(row,A.m()-1) = sum;
	}
}
