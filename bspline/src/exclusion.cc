
#include	"exclusion.hh"

//----------------------------------------------------------------------
//	Exclusion::Exclusion( const Polygons &polygons, uint n, uint m )
//----------------------------------------------------------------------

Exclusion::Exclusion( const Polygons &polygons, uint n, uint m )
{
	init( polygons, n, m );
}

//----------------------------------------------------------------------
//	void Exclusion::init( const Polygons &polygons, uint n, uint m )
//----------------------------------------------------------------------

void Exclusion::init( const Polygons &polygons, uint n, uint m )
{
	this->n = n;
	this->m = m;
	points.resize( 0, n*m );
	for( uint i = 0; i < n; i++ ) {
		for( uint j = 0; j < m; j++ ) {
			if( !is_inside(polygons, 
			   (double) i/(n-1), (double) j/(m-1))) {
				points.put( i*m+j );
			}
		}
	}
}

//----------------------------------------------------------------------
//	uint Exclusion::size() const
//----------------------------------------------------------------------

uint Exclusion::size() const
{
	return points.size();
}

//----------------------------------------------------------------------
//	bool Exclusion::is_inside( uint i, uint j ) const
//----------------------------------------------------------------------

bool Exclusion::is_inside( uint i, uint j ) const
{
	TEST_EXPR( i < n && j < m );

	return points.in( i*m+j );
}

//----------------------------------------------------------------------
//	bool Exclusion::is_outside( uint i, uint j ) const
//----------------------------------------------------------------------

bool Exclusion::is_outside( uint i, uint j ) const
{
	return !is_inside( i, j );
}

//----------------------------------------------------------------------
//	bool Exclusion::is_inside( const Polygons &polygons, 
//				  double u, double v ) const
//----------------------------------------------------------------------

bool Exclusion::is_inside( const Polygons &polygons, double u, double v ) const
{
	bool inside = false;
	for( Polygons::iterator p = polygons.begin();
	    p != polygons.end() && !inside; p++ ) {
		inside = (*p).is_inside( u, v );
	}

	return inside;
}

//----------------------------------------------------------------------
//	istream &operator >> ( istream &in, Exclusion &e )
//----------------------------------------------------------------------

istream &operator >> ( istream &in, Exclusion &e )
{
	uint n, m;
	Polygons polygons;

	in >> n >> m;
	in >> polygons;
	e.init( polygons, n, m );

	return in;
}

//----------------------------------------------------------------------
/*
int main()
{
	Exclusion excl;

	cin >> excl;
	cout << "0, 0:\t" << excl.is_inside(0,0) << endl;
	cout << "15, 15:\t" << excl.is_inside(15,15) << endl;
	cout << "80, 20:\t" << excl.is_inside(80,20) << endl;
	cout << "10, 10:\t" << excl.is_inside(10,10) << endl;
	cout << "11, 11:\t" << excl.is_inside(11,11) << endl;
	cout << "90, 60:\t" << excl.is_inside(90,60) << endl;
	cout << "21, 100:\t" << excl.is_inside(21,100) << endl;

	return 0;
}
*/
