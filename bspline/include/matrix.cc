
template< class ELEM >
Matrix<ELEM>::Matrix() : data(0), _n(0), _m(0)
{
}

template< class ELEM >
Matrix<ELEM>::Matrix( uint n, uint m ) : data(0), _n(0), _m(0)
{
	set_size( n, m );
}

template< class ELEM >
Matrix<ELEM>::Matrix( const Matrix<ELEM> &b ) : data(0), _n(0), _m(0)
{
	*this = b;
}

template< class ELEM >
Matrix<ELEM>::~Matrix()
{
	delete [] data;
}

template< class ELEM >
Matrix<ELEM> &Matrix<ELEM>::operator = ( const Matrix<ELEM> &b )
{
	if( this != &b ) {
		set_size( b._n, b._m );
		for( uint i = 0; i < _n; i++ ) {
			for( uint j = 0; j < _m; j++ ) {
				(*this)(i,j) = b(i,j);
			}
		}
	}
	return *this;
}

template< class ELEM >
void Matrix<ELEM>::set_size( uint nn, uint nm )
{
	if( _n != nn || _m != nm ) {
		delete [] data;
		data = new ELEM [nn*nm];
		_n = nn;
		_m = nm;
	}
}

template< class ELEM >
inline ELEM &Matrix<ELEM>::operator () ( uint i, uint j )
{
	TEST_EXPR( data != NULL && i < _n && j < _m );

	return data[i*_m + j];
}

template< class ELEM >
inline const ELEM &Matrix<ELEM>::operator () ( uint i, uint j ) const
{
	TEST_EXPR( data != NULL && i < _n && j < _m );

	return data[i*_m + j];
}

template< class ELEM >
Matrix<ELEM> T( const Matrix<ELEM> &a )
{
	TEST_EXPR( a.data != NULL );

	Matrix<ELEM> b(a.m(), a.n());
	for( uint i = 0; i < a.n(); i++ ) {
		for( uint j = 0; j < a.m(); j++ ) {
			b(j,i) = a(i,j);
		}
	}
	return b;
}

template< class ELEM >
Matrix<ELEM> operator + ( const Matrix<ELEM> &a, const Matrix<ELEM> &b )
{
	TEST_EXPR( a.data != NULL && b.data != NULL );
	TEST_EXPR( a._n == b._n && a._m == b._m );

	Matrix<ELEM> c( a._n, a._m );
	for( uint i = 0; i < c._n; i++ ) {
		for( uint j = 0; j < c._m; j++ ) {
			c(i,j) = a(i,j) + b(i,j);
		}
	}
	return c;
}

template< class ELEM >
Matrix<ELEM> operator * ( const Matrix<ELEM> &a, const Matrix<ELEM> &b )
{
	TEST_EXPR( a.data != NULL && b.data != NULL );
	TEST_EXPR( a._m == b._n );

	Matrix<ELEM> c( a._n, b._m );
	for( uint i = 0; i < a._n; i++ ) {
		for( uint j = 0; j < b._m; j++ ) {
			ELEM value = 0;
			for( uint k = 0; k < a._m; k++ ) {
				value += a(i,k) * b(k,j);
			}
			c(i,j) = value;
		}
	}
	return c;
}

template< class ELEM >
Matrix<ELEM> operator * ( double a, const Matrix<ELEM> &b )
{
	TEST_EXPR( b.data != NULL );

	Matrix<ELEM> c( b._n, b._m );
	for( uint i = 0; i < c._n; i++ ) {
		for( uint j = 0; j < c._m; j++ ) {
			c(i,j) = a * b(i,j);
		}
	}
	return c;
}

template< class ELEM >
ostream &operator << ( ostream &o, const Matrix<ELEM> &m )
{
	TEST_EXPR( m.data != NULL );

	o << m._n << ' ' << m._m << endl;

	int size = m._n * m._m;
	for( ELEM *p = m.data; p-m.data < size; p++ ) {
		o << *p;
		if( (p-m.data) % m._m == 0 ) {
			o << endl;
		} else {
			o << '\t';
		}
	}
	return o;
}

template< class ELEM >
istream &operator >> ( istream &input, Matrix<ELEM> &m )
{
	uint h, w;
	input >> h >> w;
	m.set_size( h, w );

	int size = m._n * m._m;
	for( ELEM *p = m.data; p-m.data < size; p++ ) {
		input >> *p;
	}
	return input;
}

template<class ELEM>
Row<ELEM>::Row( const Matrix<ELEM> &A, uint i, uint j ) : Matrix(j-i+1,A.m())
{
	for( uint k = i; k <= j; k++ ) {
		for( uint l = 0; l < A.m(); l++ ) {
			(*this)(k-i,l) = A(k,l);
		}
	}
}

template<class ELEM>
Col<ELEM>::Col( const Matrix<ELEM> &A, uint i, uint j ) 
: Matrix<ELEM>(A.n(),j-i+1)
{
	for( uint k = 0; k < A.n(); k++ ) {
		for( uint l = i; l <= j; l++ ) {
			(*this)(k,l-i) = A(k,l);
		}
	}
}

//----------------------------------------------------------------------
/*
int main()
{
	Matrix<double> a;
	Matrix<double> b(2,3);
	Matrix<double> c(3,2);

	cin >> a >> b >> c;
	cout << "b*c =\n" << b*c << endl;
	cout << "3*a =\n" << 3*a << endl;
	cout << "b+b =\n" << b+b << endl;
	cout << "a = b+b\n" << (a = b+b) << endl;
	Matrix<double> d = a;
	cout << "d(a)\n" << d << endl;

	return 0;
}

*/
