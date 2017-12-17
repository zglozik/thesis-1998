
#ifndef	_HOUSEHOLDER_HH_
#define	_HOUSEHOLDER_HH_

#include	<matrix.hh>

class HouseHolder {

    public:
	HouseHolder( double E = 1e-16 ) : E(E) { }

	void	solve( Matrix<double> &A, Matrix<double> &x );

    private:
	double	E;
};

#endif // _HOUSEHOLDER_HH_
