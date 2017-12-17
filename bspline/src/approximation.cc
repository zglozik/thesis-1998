
//
//	common member functions for Approximation class
//

#include	"approximation.hh"
#include	"householder.hh"
#include	<vtkMath.h>

//----------------------------------------------------------------------
//	init_parameters( uint un, uint vn )
//----------------------------------------------------------------------

void Approximation::init_parameters( uint un, uint vn )
{
	parameters.set_size( un, vn );
	for( uint i = 0; i < un; i++ ) {
		for( uint j = 0; j < vn; j++ ) {
			parameters(i,j).u = i;
			parameters(i,j).v = j;
		}
	}
}

//----------------------------------------------------------------------
//	fill_knotvector( KnotVector &kv, double min, double max ) const
//----------------------------------------------------------------------

void Approximation::fill_knotvector( KnotVector &kv, uint n,
				    double min, double max ) const
{
	TEST_EXPR( (uint) kv.size() >= 2*n+2 );

	uint i = 0;
	for( ; i < n+1; kv[i++] = min );
	for( ; i < kv.size()-n-1; i++ ) {
		kv[i] = min + (i-n) * (max-min) / (kv.size()-2*n-1);
	}
	for( ; i < (uint) kv.size(); kv[i++] = max );
}

//----------------------------------------------------------------------
//	fill_controls( Matrix<Point> &controls, Matrix<double> &x,
//		  uint uK, uint vK, uint index ) const
//----------------------------------------------------------------------

void Approximation::fill_controls( Matrix<Point> &controls, Matrix<double> &x,
				  uint uK, uint vK, uint index ) const
{
	for( uint i = 0; i < uK; i++ ) {
		for( uint j = 0; j < vK; j++ ) {
			controls(i,j)[index] = x(i*vK+j);
		}
	}
}

//----------------------------------------------------------------------
//	double Approximation::count_diff( const BSplineSurface &surface,
//					 const Exclusion *excl ) const
//----------------------------------------------------------------------

double Approximation::count_diff( const BSplineSurface &surface,
				 const Exclusion *excl ) const
{
	double		diff = 0;
	Point		d;
	const Par	*p;

	for( uint i = 0; i < measured.n(); i++ ) {
		for( uint j = 0; j < measured.m(); j++ ) {
			if( !excl || excl->is_inside(i,j) ) {
				p = &parameters( i, j );
				d = surface(p->u,p->v) - measured(i,j);
				diff += d * d;
			}
		}
	}
	uint num_points = excl ? excl->size() : measured.n()*measured.m();
	return diff / num_points;
}

//----------------------------------------------------------------------
//	count_matrix_sum( const Matrix<double> &bp,
//			const Base &ubase, const Base &vbase,
//			uint i, uint j, uint m, uint n,
//			const Exclusion *excl ) const
//----------------------------------------------------------------------

double Approximation::count_matrix_sum( const Matrix<double> &bp,
				       const Base &ubase, const Base &vbase,
				       uint i, uint j, uint m, uint n,
				       const Exclusion *excl ) const
{
	double sum = 0;
	for( uint k = 0; k < measured.n(); k++ ) {
		for( uint l = 0; l < measured.m(); l++ ) {
			if( !excl || excl->is_inside(k,l) ) {
				sum += base_product(bp,ubase,vbase,k,l,m,n)*
					base_product(bp,ubase,vbase,k,l,i,j);
			}
		}
	}
	return sum;
}

//----------------------------------------------------------------------
//	fill_base_product( Matrix<double> &A, const Base &ubase,
//			 const Base &vbase ) const
//----------------------------------------------------------------------

void Approximation::fill_base_product( Matrix<double> &A, const Base &ubase,
				      const Base &vbase ) const
{
	Par	p;
	for( uint i = 0; i < A.n(); i++ ) {
		uint u0 = i / measured.m();
		uint v0 = i % measured.m();
		p = parameters(u0,v0);
		for( uint j = 0; j < A.m(); j++ ) {
			uint ui = j / vbase.K();
			uint vi = j % vbase.K();
			A(i,j) = ubase(ui,p.u) * vbase(vi,p.v);
		}
	}
}

//----------------------------------------------------------------------
//	solve_problem( const Matrix<double> &A, Matrix<double> &x ) const
//----------------------------------------------------------------------

bool Approximation::solve_problem( const Matrix<double> &A, 
				  Matrix<double> &x ) const
{
	HouseHolder	hh;
	Matrix<double>	B = A;

	hh.solve( B, x );
	return true;
}
