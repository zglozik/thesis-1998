
#include	<surfgen.hh>
#include	<fstream.h>
#include	<String.h>
#include	<unistd.h>
#include	"exclusion.hh"

static char	excl_name[100];
static char	name[100];

//----------------------------------------------------------------------
//	static void usage()
//----------------------------------------------------------------------

static void usage()
{
	cerr << "usage: spl2vtk [-e <excl>] <file>\n";
	cerr << "\t-e <excl>:\texclude points, <excl> is the exclusion file\n";
	cerr << "\t<file>:\tname of spline file w/o extension\n";
}

//----------------------------------------------------------------------
//	static bool parse_options( int argc, char *argv[] )
//----------------------------------------------------------------------

static bool parse_options( int argc, char *argv[] )
{
	extern int	optind;
	extern char	*optarg;

	excl_name[0] = '\0';
	bool error = false;
	int c;
	while( (c=getopt(argc, argv, "he:")) != EOF ) {
		switch( c ) {
		    case 'e':
			strcpy( excl_name, optarg );
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
		strcpy( name, argv[optind] );
	}
	if( error ) usage();
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

	Exclusion *excl = 0;
	if( excl_name[0] ) {
		input.open( String(excl_name)+String(".excl") );
		if( !input ) {
			cerr << "can't open file: ";
			cerr << String(excl_name)+String(".excl") << endl;
			exit(1);
		}
		Polygons polygons;
		input >> polygons;
		if( !input ) {
			cerr << "input is corrupted\n";
			exit(1);
		}
		input.close();
		excl = new Exclusion(polygons, un, vn);
	}

	surfgen.set_param( un, vn );
	surfgen.write_surface( name, surface, excl );
	delete excl;
	surfgen.write_controls( name, surface );
	
	return 0;
}
