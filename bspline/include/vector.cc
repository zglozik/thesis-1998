
#include	<iostream.h>

template< class ELEM >
ostream &operator << ( ostream &s, const Vector<ELEM> &v )
{
	s << v.x << '\t' << v.y << '\t' << v.z;
	return s;
}

template< class ELEM >
istream &operator >> ( istream &s, Vector<ELEM> &v )
{
	s >> v.x;
	s >> v.y;
	s >> v.z;
	return s;
}

template< class ELEM >
bool operator < ( const Vector<ELEM> &a, const Vector<ELEM> &b )
{
	return abs(a) < abs(b) && abs(a) != abs(b);
}

template< class ELEM >
bool operator == ( const Vector<ELEM> &a, const Vector<ELEM> &b )
{
	return equal(a[0], b[0]) && equal(a[1], b[1]) && equal(a[2], b[2]);
}

template< class ELEM >
Vector<ELEM> operator +( const Vector<ELEM> &a, const Vector<ELEM> &b )
{
	return Vector<ELEM>( a.x + b.x, a.y + b.y, a.z + b.z );
}

template< class ELEM >
Vector<ELEM> operator +( const Vector<ELEM> &a )
{
	return a;
}

template< class ELEM >
Vector<ELEM> operator -( const Vector<ELEM> &a, const Vector<ELEM> &b )
{
	return Vector<ELEM>( a.x - b.x, a.y - b.y, a.z - b.z );
}

template< class ELEM >
Vector<ELEM> operator -( const Vector<ELEM> &a )
{
	return Vector<ELEM>( -a.x, -a.y, -a.z );
}

template< class ELEM >
ELEM operator *( const Vector<ELEM> &a, const Vector<ELEM> &b )
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

template< class ELEM >
Vector<ELEM> operator*( const ELEM &a, const Vector<ELEM> &b )
{
	return Vector<ELEM>( a*b.x, a*b.y, a*b.z );
}

template< class ELEM >
inline Vector<ELEM> operator*( const Vector<ELEM> &a, const ELEM &b )
{
	return b * a;
}

template< class ELEM >
inline Vector<ELEM> operator/( const Vector<ELEM> &a, const ELEM &b )
{
	return Vector<ELEM>( a.x/b, a.y/b, a.z/b );
}

template< class ELEM >
Vector<ELEM> vector_mul( const Vector<ELEM> &a, const Vector<ELEM> &b )
{
	return Vector<ELEM>( a.y*b.z - b.y*a.z,
			    a.z*b.x - b.z*a.x,
			    a.x*b.y - b.x*a.y );
}

inline double abs( const Vector<double> &v )
{
	return sqrt( v*v );
}

template< class ELEM >
double cos_alpha( const Vector<ELEM> &a, const Vector<ELEM> &b )
{
	TEST_EXPR( abs(a) != 0.0 && abs(b) != 0.0 );
	double ca = a * b / (abs(a) * abs(b));
	return min( max(ca,-1.0), 1.0);
}

template< class ELEM >
double alpha( const Vector<ELEM> &a, const Vector<ELEM> &b )
{
	double ca = cos_alpha(a, b);
	return equal(ca, 1) ? 0 : rad2deg(acos(ca));
}

template< class ELEM >
Vector<ELEM> normal( const Vector<ELEM> &v )
{
	return v == Vector<ELEM>(0) ? Vector<ELEM>(0) : (1/abs(v)) * v;
}

//----------------------------------------------------------------------
/*
int main()
{
	Vector<double>	i( 1, 0, 0 ), j( 0, 1, 0), k( 0, 0, 1);

	cout << i + j << endl;
	cout << (i+j+k) * (j + k) << endl;
	cout << vector_mul( i, j ) << endl;
	cout << vector_mul( j, k ) << endl;
	cout << vector_mul( j, i ) << endl;
	cout << abs(j + k) << endl;
	cout << -j << endl;
	cout << vector_mul(i,j) * k << endl;
	cout << cos_alpha( i, i ) << endl;
	cout << alpha( i, j ) << endl;

	return 0;
}
*/
