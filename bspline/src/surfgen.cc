
#include	"surfgen.hh"

#include	<fstream.h>
#include	<vtkFloatPoints.h>
#include	<vtkFloatScalars.h>
#include	<vtkCellArray.h>
#include	<vtkPolyData.h>
#include	<vtkStructuredGrid.h>
#include	<../graphics/vtkPolyWriter.h>
#include	<../graphics/vtkStructuredGridWriter.h>
#include	<String.h>
#include	<matrix.hh>
#include	<unistd.h>
#include	<exclusion.hh>

//----------------------------------------------------------------------
//	void SurfGen::set_param( uint un, uint vn )
//----------------------------------------------------------------------

void SurfGen::set_param( uint un, uint vn )
{
	this->un = un;
	this->vn = vn;
}

//----------------------------------------------------------------------
//	bool SurfGen::write_surface( const char *name,
//				const BSplineSurface &s,
//				const Exclusion *excl ) 
//----------------------------------------------------------------------

bool SurfGen::write_surface( const char *name, const BSplineSurface &s,
			    const Exclusion *excl )
{
	String fname;

	bool error;
	vtkPolyData *spolygon = new vtkPolyData;	// for surface
	vtkPolyData *ppolygon = new vtkPolyData;	// for param. lines

	cout << "Scanning surface..."; cout.flush();

	error = make_points( spolygon, ppolygon, s );
	make_cells( spolygon, ppolygon, excl );

	fname = String(name) + String(".par");
	cout << "done.\nWriting file: " << fname << "..."; cout.flush();

	vtkPolyWriter *writer = new vtkPolyWriter;
	writer->SetFilename( (char *) (const char *) fname );
	writer->SetInput( ppolygon );
	writer->Write();

	fname = String(name) + String(".surf");
	cout << "done.\nWriting file: " << fname << "..."; cout.flush();

	error = error || make_scalars( spolygon, s, 0 );// minkowski curvature
	writer->SetFilename( (char *) (const char *) fname );
	writer->SetInput( spolygon );
	writer->Write();
/*
	fname = String(name) + String(".gauss");
	cout << "done.\nWriting file: " << fname << "..."; cout.flush();

	error = error || make_scalars( spolygon, s, 1 );// gauss curvature
	writer->SetFilename( (char *) (const char *) fname );
	writer->SetInput( spolygon );
	writer->Write();
*/
	writer->Delete();
	cout << "done.\n";
	if( error ) {
		cerr << endl << "error during computation of surface\n";
	}

	return error;
}

//----------------------------------------------------------------------
//	bool SurfGen::make_points(vtkPolyData *spolygon, vtkPolyData *ppolygon,
//				  const BSplineSurface &s )
//----------------------------------------------------------------------

bool SurfGen::make_points( vtkPolyData *spolygon, vtkPolyData *ppolygon,
			  const BSplineSurface &s )
{
	vtkFloatPoints *points = new vtkFloatPoints( un*vn );
	double u, v;
	float ap[3];
	Point p;

	for( uint ui = 0; ui < un; ui++ ) {
		u = s.umin_t() + ui * (s.umax_t() - s.umin_t()) / (un-1);
		for( uint vi = 0; vi < vn; vi++ ) {
			v = s.vmin_t() + vi * (s.vmax_t()-s.vmin_t()) / (vn-1);
			p = s( u, v );	// this op. is expensive!
			ap[0] = p[0];
			ap[1] = p[1];
			ap[2] = p[2];
			points->InsertNextPoint( ap );
		}
	}
	spolygon->SetPoints( points );
	ppolygon->SetPoints( points );
	points->Delete();

	return false;
}

//----------------------------------------------------------------------
//	bool SurfGen::make_scalars( vtkPolyData *polygon, 
//				const BSplineSurface &s, int curvature )
//----------------------------------------------------------------------

bool SurfGen::make_scalars( vtkPolyData *polygon, const BSplineSurface &s,
			   int curvature )
{
	vtkFloatScalars *scalars = new vtkFloatScalars( un*vn );
	double u, v;
	bool error = false;

	for( uint ui = 0; ui < un; ui++ ) {
		u = s.umin_t() + ui * (s.umax_t() - s.umin_t()) / (un-1);
		for( uint vi = 0; vi < vn; vi++ ) {
			v = s.vmin_t() + vi * (s.vmax_t()-s.vmin_t()) / (vn-1);
			if( !curvature ) {
				scalars->InsertNextScalar( 
					s.minkowski(u, v, error) );
			} else {
				scalars->InsertNextScalar( 
					s.gaussian(u, v, error) );
			}
		}
	}
	polygon->GetPointData()->SetScalars( scalars );
	scalars->Delete();

	return false;
}

//----------------------------------------------------------------------
//	void SurfGen::make_cells( vtkPolyData *spolygon, vtkPolyData *ppolygon,
//				const Exclusion *excl )
//----------------------------------------------------------------------

void SurfGen::make_cells( vtkPolyData *spolygon, vtkPolyData *ppolygon,
			 const Exclusion *excl )
{
	vtkCellArray *cells = new vtkCellArray;
	uint count = 0;
	bool pushing;
	for( uint ui = 0; ui < un-1; ui++ ) {
		pushing = false;
		for( uint vi = 0; vi < vn; vi++ ) {
			if( !excl || 
			   (excl->is_inside(ui,vi) && 
			    excl->is_inside(ui+1,vi)) ) {
				if( !pushing ) {
					pushing = true;
					count = 0;
					cells->InsertNextCell( 2*vn );
				}
				cells->InsertCellPoint( ui*vn + vi );
				cells->InsertCellPoint( (ui+1)*vn + vi );
				count += 2;
			} else if( pushing ) {
				pushing = false;
				cells->UpdateCellCount(count);
			}
		}
		if( pushing ) cells->UpdateCellCount(count);
	}	
	spolygon->SetStrips( cells );
	cells->Delete();

	cells = new vtkCellArray;
	for( uint i = 0; i <= 15; i++ ) {
		uint ui = uint(i * (un-1)/15.0 + 0.5);
		cells->InsertNextCell( vn );
		for( uint vi = 0; vi < vn; vi++ ) {
			cells->InsertCellPoint( ui*vn + vi );
		}
	}
	for( uint i = 0; i <= 15; i++ ) {
		uint vi = uint(i * (vn-1)/15.0 + 0.5);
		cells->InsertNextCell( un );
		for( uint ui = 0; ui < un; ui++ ) {
			cells->InsertCellPoint( ui*vn + vi );
		}
	}
	ppolygon->SetLines( cells );
	cells->Delete();
}

//----------------------------------------------------------------------
//	bool SurfGen::write_real( const char *name, double maxerr,
//				 const BSplineSurface &s )
//----------------------------------------------------------------------

bool SurfGen::write_real( const char *name, double maxerr,
			 const BSplineSurface &s )
{
	String fname = String(name) + String(".real");
	ofstream real( fname );
	if( !real ) {
		cerr << "can't write file: " << fname << endl;
		return false;
	}
	cout << "Writing file: " << fname << "..."; cout.flush();

	double u, v;
	Point p;
	Matrix<Point> measured(un, vn);

	for( uint ui = 0; ui < un; ui++ ) {
		u = s.umin_t() + ui * (s.umax_t() - s.umin_t()) / (un-1);
		for( uint vi = 0; vi < vn; vi++ ) {
			v = s.vmin_t() + vi * (s.vmax_t()-s.vmin_t()) / (vn-1);
			p = s( u, v );
			measured(ui,vi) = mess(p,maxerr);
		}
	}

	real << measured;
	real.close();
	cout << "done.\n";

	return true;
}

//----------------------------------------------------------------------
//	bool SurfGen::write_controls( const char *name,
//				const BSplineSurface &s )
//----------------------------------------------------------------------

bool SurfGen::write_controls( const char *name, const BSplineSurface &s )
{
	String fname = String(name) + String(".cp");
	cout << "Writing file: " << fname << "..."; cout.flush();

	vtkFloatPoints *points = new vtkFloatPoints( s.uK()*s.vK() );
	float ap[3];
	for( uint ui = 0; ui < s.uK(); ui++ ) {
		for( uint vi = 0; vi < s.vK(); vi++ ) {
			ap[0] = s.get_points()(ui,vi)[0];
			ap[1] = s.get_points()(ui,vi)[1];
			ap[2] = s.get_points()(ui,vi)[2];
			points->InsertNextPoint( ap );
		}
	}

	vtkStructuredGrid *grid = new vtkStructuredGrid;
	grid->SetDimensions( s.uK(), s.vK(), 1 );
	grid->SetPoints( points );
	points->Delete();

	vtkStructuredGridWriter *writer = new vtkStructuredGridWriter;
	writer->SetFilename( (char *) (const char *) fname );
	writer->SetInput( grid );
	writer->Write();
	writer->Delete();
	cout << "done.\n";

	return true;
}

//----------------------------------------------------------------------
//	Point SurfGen::mess( const Point &source, double maxerr )
//----------------------------------------------------------------------

Point SurfGen::mess( const Point &source, double maxerr )
{
	double alpha = M_PI * rand()/RAND_MAX;
	double beta = 2*M_PI * rand()/RAND_MAX;
	Point dir( sin(alpha)*sin(beta),
		  cos(alpha),
		  sin(alpha)*cos(beta) );
	return source + maxerr*rand()/RAND_MAX * dir;
}
