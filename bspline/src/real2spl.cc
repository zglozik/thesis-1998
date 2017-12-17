
#include	<approximation.hh>
#include	<fstream.h>
#include	<String.h>
#include	<unistd.h>

static int	u, v, U, V, n, m;
static double	max_err;
static char	excl_name[100];
static char	name[100];
static int	version;

static double approximate( Matrix<Point> *measured, BSplineSurface &surface,
			   const Exclusion *excl )
{
	Approximation	approx(*measured);
	Q1	q1;	// functional 1
	Q2	q2;	// functional 2
	Q3	q3;	// functional 3
	double	rv;

	switch( version ) {
	    case 1:
		rv = approx.approximate1( U, V, n, m, surface, true, excl );
		break;
	    case 2:
		rv = approx.approximate2( U, V, n, m, surface, true, excl );
		break;
	    case 3:
		rv = approx.approximate3( U, V, n, m, surface,
					 &Approximation::correct_parameter,
					 excl );
		break;
	    case 4:
		rv = approx.approximate3( U, V, n, m, surface,
					 &Approximation::correct_parameter2,
					 excl );
		break;
	    case 5:
		rv = approx.approximate3( U, V, n, m, surface,
					 &Approximation::correct_parameter3,
					 excl );
		break;
	    case 6:
		rv = approx.approximate4( U, V, n, m, surface, q1, max_err );
		break;
	    case 7:
		rv = approx.approximate4( U, V, n, m, surface, q2, max_err );
		break;
	    case 8:
		rv = approx.approximate4( U, V, n, m, surface, q3, max_err );
		break;
	    default:
		ERROR( "bad version parameter" );
	}
	return rv;
}

static bool parse_options( int argc, char *argv[] )
{
	extern int	optind;
	extern char	*optarg;

	excl_name[0] = '\0';
	u = 40;
	v = 40;
	U = 10;
	V = 10;
	n = 3;
	m = 3;
	version = 1;
	bool error = false;
	int c;
	while( (c=getopt(argc, argv, "12345678he:u:v:U:V:n:m:r:")) != EOF ) {
		switch( c ) {
		    case '1':
		    case '2':
		    case '3':
		    case '4':
		    case '5':
		    case '6':
		    case '7':
		    case '8':
			version = c - '0';
			break;
		    case 'h':
			error = true;
			break;
		    case 'e':
			strcpy(excl_name,optarg);
			break;
		    case 'u':
			u = atoi(optarg);
			break;
		    case 'v':
			v = atoi(optarg);
			break;
		    case 'U':
			U = atoi(optarg);
			break;
		    case 'V':
			V = atoi(optarg);
			break;
		    case 'n':
			n = atoi(optarg);
			break;
		    case 'm':
			m = atoi(optarg);
			break;
		    case 'r':
			max_err = atof(optarg);
			break;
		    case '?':
		    default:
			error = true;
			break;			
		}
	}
	if( optind+1 != argc ) {
		error = true;
	} else {
		strcpy(name,argv[optind]);
	}
	if( error ) {
		cerr << "usage: real2spl [-U num] [-V num] [-u num] [-v num] ";
		cerr << "[-n num] [-m num] [-r num] ";
		cerr << "[-e excl_file] [-1|-2|-3|-4|-5|-6|-7|-8] <file>\n";
		cerr << "\t-U num: number of control points in u ";
		cerr << "direction(10)\n";
		cerr << "\t-V num: number of control points in v ";
		cerr << "direction(10)\n";
		cerr << "\t-u num: triangles to draw in u direction(40)\n";
		cerr << "\t-v num: triangles to draw in v direction(40)\n";
		cerr << "\t-n num: order of b-splines in u direction(3)\n";
		cerr << "\t-m num: order of b-splines in v direction(3)\n";
		cerr << "\t-e excl_file: exclusion file\n";
		cerr << "\t-r number: maximum error allowed\n";
		cerr << "\t-1: using rectangular matrix(default)\n";
		cerr << "\t-2: using square matrix\n";
		cerr << "\t-3: using parameter corrections (first method)\n";
		cerr << "\t-4: using parameter corrections (second method)\n";
		cerr << "\t-5: using parameter corrections (third method)\n";
		cerr << "\t-6: also using smoothing functional 1\n";
		cerr << "\t-7: also using smoothing functional 2\n";
		cerr << "\t-8: also using smoothing functional 3\n";
		cerr << "\t<file>: real file\n";
	}
	return error;
}

int main( int argc, char **argv )
{
	if( parse_options(argc,argv) ) exit(1);

	String fname;
	fname = String(name) + String(".real");
	ifstream input( fname );
	if( !input ) {
		cerr << "can't open file: " << fname << endl;
		exit(1);
	}
	
	Matrix<Point> measured;
	input >> measured;
	if( !input ) {
		cerr << "can't read matrix from " << fname << endl;
		exit(1);
	}
	input.close();

	Exclusion *excl = 0;
	if( excl_name[0] ) {
		fname = String(excl_name) + String(".excl");
		input.open( fname );
		if( !input ) {
			cerr << "can't open file: " << fname << endl;
			exit(1);
		}
		Polygons polygons;
		input >> polygons;
		if( !input ) {
			cerr << "can't read polygons from " << fname << endl;
			exit(1);
		}
		input.close();
		excl = new Exclusion(polygons, measured.n(), measured.m());
	}

	fname = String(name) + String(".spl");
	ofstream out( fname );
	if( !out ) {
		cerr << "can't write file: " << fname << endl;
		exit(1);
	}
	cout << "Writing file: " << fname << "...\n";

	BSplineSurface surface;
	double diff = approximate( &measured, surface, excl );

	out << u << endl;
	out << v << endl;
	out << surface;
	out.close();
	cout << "done.\n";

	cout << "average of squares of differences: " << diff << endl;

	return 0;
}
