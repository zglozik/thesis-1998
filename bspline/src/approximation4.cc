
#include	"approximation.hh"
#include	"householder.hh"
#include	<vtkMath.h>

//----------------------------------------------------------------------
//	approximate4( uK, vK, un, vn, surface, functional, max_err,
//			default_parameters, excl ) 
//----------------------------------------------------------------------

double Approximation::approximate4( uint uK, uint vK, uint un, uint vn,
				   BSplineSurface &surface,
				   Functional &functional,
				   double max_err,
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
	
	Matrix<double> base_product( measured.n()*measured.m(), uK*vK );

	fill_base_product( base_product, surface.ubase(), surface.vbase() );
	functional.init( surface.ubase(), surface.vbase() );

	double	x0, x1, x2, fx2;
	double	fx0, fx1;
	uint	step = 0;
/*
	x0 = 1;
	while( x0 >= 0 ) {
		fx0 = solve_one_step( step++, surface, base_product,
				     functional, x0, excl ) - max_err;
		x0 -= 0.1;
	}
*/

	x0 = 0.5;
	fx0 = solve_one_step( ++step, surface, base_product, functional,
			     x0, excl ) - max_err;
	x1 = fx0 < 0 ? 0 : 1;
	fx1 = solve_one_step( ++step, surface, base_product, functional,
			     x1, excl ) - max_err;
	while( ++step <= 15 ) {
//		x2 = -fx1 * (x1-x0) / (fx1-fx0) + x1;
		x2 = (x0 + x1) / 2.0;
		fx2 = solve_one_step( step, surface, base_product,
				     functional, x2, excl ) - max_err;
		if( fx1*fx2 < 0 ) {
			x0 = x2;
			fx0 = fx2;
		} else {
			x1 = x2;
			fx1 = fx2;
		}
	}
	return fx2 + max_err;
}

//----------------------------------------------------------------------
//	solve_one_step( uint step, BSplineSurface &surface,
//		     const Functional &functional,
//		     double u, const Exclusion *excl )
//----------------------------------------------------------------------

double Approximation::solve_one_step( uint step, BSplineSurface &surface,
				     const Matrix<double> &base_product,
				     const Functional &functional,
				     double u, const Exclusion *excl )
{
	const uint n = surface.uK() * surface.vK();
	Matrix<double> A( n, n+1 );
	Matrix<double> x( n );
	Matrix<Point> controls( surface.uK(), surface.vK() );

	cout << "Step " << step << ". (" << 1-u << "): ...";
	cout.flush();

	fill_matrix4( A, base_product, functional, u,
		     surface.ubase(), surface.vbase(), excl );
	for( uint index = 0; index < 3; index++ ) {
		fill_b4( A, base_product, index, u,
			surface.ubase(), surface.vbase(), excl );
		solve_problem( A, x );
		fill_controls( controls, x, surface.uK(), surface.vK(),
			      index );
	}
	surface.set_points( controls );

	double err = count_diff( surface, excl );
	cout << "\naverage of squares of differences: " << err << endl;

	return err;
}

//----------------------------------------------------------------------
//	fill_matrix4( Matrix<double> &A, const Matrix<double> &base_product,
//		const Functional &functional, double u,
//		const Base &ubase, const Base &vbase,
//		const Exclusion *excl ) const
//----------------------------------------------------------------------

void Approximation::fill_matrix4( Matrix<double> &A, 
				 const Matrix<double> &base_product,
				 const Functional &functional, double u,
				 const Base &ubase, const Base &vbase,
				 const Exclusion *excl ) const
{
	for( uint row = 0; row < A.n(); row++ ) {
		uint i = row / vbase.K();
		uint j = row % vbase.K();
		for( uint col = 0; col < A.m()-1; col++ ) {
			uint m = col / vbase.K();
			uint n = col % vbase.K();
			A(row,col) = u * 2 * count_matrix_sum( base_product,
				ubase, vbase, i, j, m, n, excl );
			A(row,col) += (1-u) * functional.derivate(i,j,m,n);
		}
	}
}

//----------------------------------------------------------------------
//	fill_b4( Matrix<double> &A, const Matrix<double> &bp, uint index,
//		double u, const Base &ubase, const Base &vbase,
//		const Exclusion *excl ) const
//----------------------------------------------------------------------

void Approximation::fill_b4( Matrix<double> &A, const Matrix<double> &bp,
			    uint index, double u,
			    const Base &ubase, const Base &vbase,
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
				sum += u * 2 * measured(k,l)[index] *
					base_product(bp,ubase,vbase,k,l,i,j);
			}
		}
		A(row,A.m()-1) = sum;
	}
}

