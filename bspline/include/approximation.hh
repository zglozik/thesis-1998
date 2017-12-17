
#ifndef	_APPROXIMATION_HH_
#define	_APPROXIMATION_HH_

#include	"bspline.hh"
#include	"matrix.hh"
#include	"exclusion.hh"
#include	"functional.hh"

class Approximation {

	typedef	void (Approximation::*comp_fun)( const BSplineSurface &, 
						uint, uint );

    public:
	Approximation( const Matrix<Point> &measured ) : measured(measured) { }

	double approximate1( uint uK, uint vK, uint un, uint vn,
			    BSplineSurface &surface,
			    bool default_parameters = true,
			    const Exclusion *excl = 0 );
	double approximate2( uint uK, uint vK, uint un, uint vn,
			    BSplineSurface &surface,
			    bool default_parameters = true,
			    const Exclusion *excl = 0 );
	double approximate3( uint uK, uint vK, uint un, uint vn,
			    BSplineSurface &surface,
			    comp_fun act_fun, const Exclusion *excl = 0);
	double approximate4( uint uK, uint vK, uint un, uint vn,
			    BSplineSurface &surface, Functional &functional,
			    double max_err,
			    bool default_parameters = true,
			    const Exclusion *excl = 0 );
	void correct_parameter( const BSplineSurface &surface,
			       uint i, uint j );
	void correct_parameter2( const BSplineSurface &surface,
				uint i, uint j );
	void correct_parameter3( const BSplineSurface &surface,
				uint i, uint j );

    private:
	struct Par {
		Par( double u = 0, double v = 0 ) : u(u), v(v) { }
		double	u;
		double	v;
	};

	// common member functions

	void init_parameters( uint un, uint vn );
	bool solve_problem( const Matrix<double> &A, Matrix<double> &x ) const;
	void fill_knotvector( KnotVector &kv, uint n, 
			     double min, double max ) const;
	double count_diff( const BSplineSurface &surface,
			  const Exclusion *excl ) const;
	void fill_controls( Matrix<Point> &controls, Matrix<double> &x,
			   uint uK, uint vK, uint index ) const;
	void fill_base_product( Matrix<double> &A,
			       const Base &ubase, const Base &vbase ) const;
	inline double base_product( const Matrix<double> &A,
				   const Base &ubase, const Base &vbase,
				   uint u0, uint v0, uint i, uint j ) const;
	double count_matrix_sum( const Matrix<double> &bp, 
				const Base &ubase, const Base &vbase,
				uint i, uint j, uint m, uint n,
				const Exclusion *excl = 0 ) const;

	// member functions for version 1

	void fill_matrix1( Matrix<double> &A, 
			  const Base &ubase, const Base &vbase,
			  const Exclusion *excl = 0 ) const ;
	void fill_b1( Matrix<double> &A, uint index, 
		     const Exclusion *excl = 0 ) const;

	// member functions for version 2

	void fill_matrix2( Matrix<double> &A, const Matrix<double> &bp,
			  const Base &ubase, const Base &vbase,
			  const Exclusion *excl = 0 ) const ;
	void fill_b2( Matrix<double> &A, const Matrix<double> &bp, uint index, 
		     const Base &ubase, const Base &vbase,
		     const Exclusion *excl = 0 ) const;

	// member functions for version 3

	void compute_length( const BSplineSurface &surface,
			    uint i, uint j, double &lu, double &lv ) const;
	void correct_parameters( const BSplineSurface &surface,
				const Exclusion *excl, comp_fun act_fun );

	// member functions for version 4

	double solve_one_step( uint step, BSplineSurface &surface,
			      const Matrix<double> &base_product,
			      const Functional &functional,
			      double u, const Exclusion *excl );
	void fill_matrix4( Matrix<double> &A, 
			  const Matrix<double> &base_product,
			  const Functional &functional, double u,
			  const Base &ubase, const Base &vbase,
			  const Exclusion *excl = 0 ) const;
	void fill_b4( Matrix<double> &A, const Matrix<double> &bp, uint index, 
		     double u, const Base &ubase, const Base &vbase,
		     const Exclusion *excl = 0 ) const;

	const Matrix<Point>	&measured;
	Matrix<Par>		parameters;
};

inline double Approximation::base_product( const Matrix<double> &A,
				const Base &ubase, const Base &vbase,
				uint u0, uint v0, uint i, uint j ) const
{
	return A( u0 * measured.m() + v0, i * vbase.K() + j );
}

#endif // _APPROXIMATION_HH_
