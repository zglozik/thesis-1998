
#ifndef	_FUNCTIONAL_HH_
#define	_FUNCTIONAL_HH_

#include	"matrix.hh"
#include	"bspline.hh"

//----------------------------------------------------------------------
//	generic functional class
//----------------------------------------------------------------------

class Functional {

    public:
	virtual ~Functional() { };

	virtual void	init( const Base &ubase, const Base &vbase ) = 0;
	virtual double	derivate( uint k, uint l, uint i, uint j ) const = 0;
	virtual double	integrate( const BSplineSurface &surface ) const;

    protected:
	void compute_integrate( Matrix<double> &A, const Base &base,
			       uint d1, uint d2 ) const;
	virtual	double	value( const BSplineSurface &surface, 
			      double u, double v ) const = 0;

};

class Q1 : public Functional {

    public:
	virtual	void	init( const Base &ubase, const Base &vbase );
	virtual double	derivate( uint k, uint l, uint i, uint j ) const;

    private:
	virtual	double	value( const BSplineSurface &surface, 
			      double u, double v ) const;

	Matrix<double>	_U11;
	Matrix<double>	_U00;
	Matrix<double>	_V11;
	Matrix<double>	_V00;

};

class Q2 : public Functional {

    public:
	virtual	void	init( const Base &ubase, const Base &vbase );
	virtual double	derivate( uint k, uint l, uint i, uint j ) const;

    private:
	virtual	double	value( const BSplineSurface &surface, 
			      double u, double v ) const;

	Matrix<double>	_U22;
	Matrix<double>	_U20;
	Matrix<double>	_U02;
	Matrix<double>	_U00;

	Matrix<double>	_V22;
	Matrix<double>	_V20;
	Matrix<double>	_V02;
	Matrix<double>	_V00;

};

class Q3 : public Functional {

    public:
	virtual	void	init( const Base &ubase, const Base &vbase );
	virtual double	derivate( uint k, uint l, uint i, uint j ) const;

    private:
	virtual	double	value( const BSplineSurface &surface, 
			      double u, double v ) const;

	uint	_uK, _vK;	// number of control points in two directions
};

#endif
