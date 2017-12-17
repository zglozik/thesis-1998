
#ifndef	_MATRIX_HH_
#define	_MATRIX_HH_

#include	<iostream.h>
#include	<debug.hh>
#include	<stdtypes.h>

template< class ELEM >
class Matrix;

template< class ELEM >
Matrix<ELEM> T( const Matrix<ELEM> &a );

template< class ELEM >
class Matrix {

friend	Matrix T<>( const Matrix<ELEM> &a );
friend	Matrix operator + <>( const Matrix &a, const Matrix &b );
friend	Matrix operator * <>( const Matrix &a, const Matrix &b );
friend	Matrix operator * <>( double a, const Matrix &b );

friend	ostream &operator << <>( ostream &o, const Matrix &m );
friend	istream &operator >> <>( istream &input, Matrix &m );

    public:
	Matrix();
	Matrix( uint n, uint m = 1 );
	Matrix( const Matrix &b );
virtual	~Matrix();

	Matrix &operator = ( const Matrix &b );

	void set_size( uint n, uint m );
virtual	ELEM &operator () ( uint i, uint j = 0 );
virtual	const ELEM &operator () ( uint i, uint j = 0 ) const;
	uint n() const { return _n; }
	uint m() const { return _m; }

	operator ELEM * () const { return data; }

    private:
	ELEM	*data;
	uint	_n, _m;
};

template< class ELEM >
class Row : public Matrix<ELEM> {

    public:
	Row( const Matrix<ELEM> &A, uint i, uint j );

};

template< class ELEM >
class Col : public Matrix<ELEM> {

    public:
	Col( const Matrix<ELEM> &A, uint i, uint j );

};

template< class ELEM >
Matrix<ELEM> T( const Matrix<ELEM> &a );
template< class ELEM >
Matrix<ELEM> operator + ( const Matrix<ELEM> &a, const Matrix<ELEM> &b );
template< class ELEM >
Matrix<ELEM> operator * ( const Matrix<ELEM> &a, const Matrix<ELEM> &b );
template< class ELEM >
Matrix<ELEM> operator * ( double a, const Matrix<ELEM> &b );

template< class ELEM >
inline Matrix<ELEM> operator * ( const Matrix<ELEM> &a, double b )
{
	return b * a;
}

template< class ELEM >
inline Matrix<ELEM> operator / ( const Matrix<ELEM> &a, double b )
{
	return (1/b) * a;
}

template< class ELEM >
extern ostream &operator << ( ostream &o, const Matrix<ELEM> &m );
template< class ELEM >
extern istream &operator >> ( istream &input, Matrix<ELEM> &m );

// inlcude implementation
#include	"matrix.cc"

#endif // _MATRIX_HH_
