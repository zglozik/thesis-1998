
#include	"householder.hh"
#include	<math.h>

//----------------------------------------------------------------------
//	solve( Matrix<double> &A, Matrix<double> &x )
//----------------------------------------------------------------------

void HouseHolder::solve( Matrix<double> &A, Matrix<double> &x )
{
	TEST_EXPR( A.n() >= A.m()-1 && x.n() == A.m()-1 );

	// we should compute R, that QR = A
	double	r, s, v;

	uint n = A.m()-1;
	uint m = A.n();
	for( uint j = 0, k = 0; j < n+1 && k < n; k++ ) {
		r = 0;
		for( uint i = j; i < m; i++ ) r += A(i,k) * A(i,k);
		if( r > E ) {
			x(k) = 1;
			r = sqrt(r);
			if( A(j,k) < 0 ) r *= -1;
			A(j,k) += r;
			s = 1 / (r*A(j,k));
			for( uint l = k+1; l < n+1; l++ ) {
				v = 0;
				for( uint i = j; i < m; i++ ) 
					v += A(i,k) * A(i,l);
				v *= s;
				for( uint i = j; i < m; i++ )
					A(i,l) -= v * A(i,k);
			}
			A(j,k) = -r;
			j++;
		} else {
			x(k) = 0;
		}
	}

	// we are solving Rx = r

	uint k = 0;
	for( uint i = 0; i < n; i++ ) k += (uint) x(i);

	for( int i = n-1; i >= 0; i-- ) {
		if( x(i) != 0 ) {
			k--;
			r = A(k,n);
			for( uint j = i+1; j < n; j++ ) r -= A(k,j) * x(j);
			x(i) = r / A(k,i);
		}
	}
}

//----------------------------------------------------------------------
/*
int main()
{
	HouseHolder hh;

	Matrix<double>	A;

	cin >> A;
	Matrix<double> x(A.m()-1);
	Matrix<double> B = Col<double>(A,0,A.m()-2);

	hh.solve( A, x );
	cout << x << endl;
	cout << B * x << endl;

	return 0;
}
*/
