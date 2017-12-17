
#include	"matrix.hh"
#include	"bspline.hh"
#include	<math.h>
#include	<iostream.h>

int main () {
	Matrix<Point>	a(20, 20);
	
	for(int i = 0; i<20; i++) {
		for(int j = 0; j<20; j++) {
			a(i,j) = Point( double(i), double(j),  
					cos( ((double(i-10)/20)*(double(i-10)/20)+
						 (double(j-10)/20)*(double(j-10)/20)) * 16 * M_PI ) );
		}
	}

	cout << a;
	return 0;
}
