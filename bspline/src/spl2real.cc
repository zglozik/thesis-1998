
#include	<surfgen.hh>
#include	<String.h>
#include	<fstream.h>
#include	<unistd.h>

static double	maxerr;
static char	*name;

//----------------------------------------------------------------------
//	static bool parse_options( int argc, char *argv[] )
//----------------------------------------------------------------------

static bool parse_options( int argc, char *argv[] )
{
	extern char *optarg;
	extern int   optind;
	int c;
	bool error = false;
	
	maxerr = 1e-3;

	while( (c=getopt(argc,argv,"he:")) != EOF ) {
		switch( c ) {
		    case 'e':
			maxerr = atof(optarg);
			break;
		    case 'h':
		    case '?':
		    default:
			error = true;
			break;
		}
	}
	if( optind+1 != argc ) {
		error = true;
	} else {
		name = argv[optind];
	}
	if( error ) {
		cerr << "usage: spl2real [-e <number>] <file>\n";
		cerr << "\t-e num:\tmaximum error (1e-3)\n";
		cerr << "\t<file>:\tspline file\n";
	}
	return error;
}

//----------------------------------------------------------------------
//	main
//----------------------------------------------------------------------

int main(int argc, char **argv)
{
	SurfGen		surfgen;
	BSplineSurface	surface;
	uint	un, vn;
	
	if( parse_options(argc,argv) ) exit(1);

 	ifstream input( String(name)+String(".spl") );
	if( !input ) {
		cerr << "can't open file: " << String(name)+String(".spl");
		cerr << endl;
		exit(1);
	}
	input >> un;
	input >> vn;
	input >> surface;
	if( !input ) {
		cerr << "input is corrupted\n";
		exit( 1 );
	}
	input.close();

	surfgen.set_param( un, vn );
	surfgen.write_real( name, maxerr, surface );
	
	return 0;
}
