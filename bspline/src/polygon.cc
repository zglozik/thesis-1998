
#include	"polygon.hh"
#include	<vtkMath.h>

//----------------------------------------------------------------------
//	bool Polygon::is_inside( double x, double y ) const
//----------------------------------------------------------------------

bool Polygon::is_inside( double x, double y ) const
{
	if( _points.size() < 3 ) return false;

	bool inside = false;
	uint intersection = 0;
	Points::iterator current = _points.begin();
	Points::iterator prev = _points.end();
	--prev;
	for( ; current != _points.end() && !inside; prev = current++ ) {
		if( is_on(x,y,(*prev).x,(*prev).y,
			  (*current).x, (*current).y) ) {
			inside = true;
		} else if( intersect(x,y,x,ty-1,(*prev).x,(*prev).y,
				     (*current).x,(*current).y) ) {
			intersection++;
		}
	}
	return (intersection % 2) || inside;
}

//----------------------------------------------------------------------
//	bool Polygon::is_outside( double x, double y ) const 
//----------------------------------------------------------------------

bool Polygon::is_outside( double x, double y ) const 
{
	return !is_inside(x,y);
}

//----------------------------------------------------------------------
//	bool Polygon::intersect( double a1, double a2, double b1, double b2,
//				double c1, double c2, double d1, double d2 )
//----------------------------------------------------------------------

bool Polygon::intersect( double a1, double a2, double b1, double b2,
			double c1, double c2, double d1, double d2 ) const
{
	double	A0[2], A1[2];
	double	*A[2] = { A0, A1 };
	double	x[2];

	A[0][0] = a1 - b1;
	A[0][1] = d1 - c1;
	A[1][0] = a2 - b2;
	A[1][1] = d2 - c2;
	x[0] = d1 - b1;
	x[1] = d2 - b2;

	bool inside;
	if( !vtkMath::SolveLinearSystem(A,x,2) ) {
		inside = false;
	} else {
		inside = less(0,x[0]) && less(x[0],1) &&
			less(0,x[1]) && less(x[1],1);
	}

	return inside;
}

//----------------------------------------------------------------------
//	bool Polygon::is_on( double x1, double x2, 
//			    double a1, double a2, double b1, double b2 ) const
//----------------------------------------------------------------------

bool Polygon::is_on( double x1, double x2, 
		    double a1, double a2, double b1, double b2 ) const
{
	double A[2] = { a1-b1, a2-b2 };
	double b[2] = { x1-b1, x2-b2 };
	bool on;

	if( (equal(A[0],0) && !equal(b[0],0)) ||
	    (equal(A[1],0) && !equal(b[1],0)) ) {
		on = false;
	} else if( equal(A[0],0) ) {
		on = equal(A[1],0) || (0<=b[1]/A[1] && b[1]/A[1]<=1);
	} else if( equal(A[1],0) ) {
		on = equal(A[0],0) || (0<=b[0]/A[0] && b[0]/A[0]<=1);
	} else {
		double x[2] = { b[0]/A[0], b[1]/A[1] };

		on = equal(x[0],x[1]) && 0<=x[0] && x[0]<=1;
	}

	return on;
}

//----------------------------------------------------------------------
/*
#include	<fstream.h>

int main()
{
	Polygon	polygon;

	ifstream input( "input" );
	if( !input ) {
		cerr << "can't open file\n";
		exit(1);
	}
	input >> polygon;
	input.close();

	cout << polygon.is_inside( 0.5, 0.5 ) << endl;
	cout << polygon.is_inside( 0.79, 0.6 ) << endl;
	cout << polygon.is_inside( 0, 0 ) << endl;
	cout << polygon.is_inside( 0, 1.7 ) << endl;
	cout << polygon.is_inside( 0, 0.7 ) << endl;
	cout << polygon.is_inside( 0.6, 0.8 ) << endl;

	return 0;
}
*/
