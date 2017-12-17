
#include	<iostream>
#include	<fstream>
#include	"bspline.hh"

static double count_diff( Matrix<Point> &p1, Matrix<Point> &p2 )
{
	double	sum = 0;
	for( uint i = 0; i < p1.n(); ++i ) {
		for( uint j = 0; j < p1.m(); ++j ) {
			sum += abs( p1(i,j) - p2(i,j) );
		}
	}

	return sum / (p1.n()*p1.m());
}

static bool open_file( ifstream &str, const char *name )
{
	str.open( name );
	if( !str ) {
		cerr << "can't open file: " << name << endl;
		return false;
	}
	return true;
}

int main( int argc, char *argv[] )
{
	if( argc != 3 ) {
		cerr << "usage: cmpspl <spline1> <spline2>" << endl;
		exit( 1 );
	}
	ifstream f1, f2;
	if( !open_file(f1, argv[1]) || !open_file(f2, argv[2]) ) {
		exit( 1 );
	}
	uint n, m;
	BSplineSurface	s1, s2;
	f1 >> n >> m >> s1;
	f2 >> n >> m >> s2;
	if( !f1 || !f2 ) {
		cerr << "can't read splines" << endl;
		exit( 1 );
	}

	Matrix<Point> &p1 = s1.get_points();
	Matrix<Point> &p2 = s2.get_points();
	if( p1.n() != p2.n() || p1.m() != p2.m() ) {
		cerr << "the two splines don't consist of the same amount of"
			" control points" << endl;
		exit( 1 );
	}
	cout << "average of squares of differences: ";
	cout << count_diff( p1, p2 ) << endl;

	return 0;
}
