
#ifndef	_SURFGEN_HH_
#define	_SURFGEN_HH_

#include	"bspline.hh"

#include	<time.h>
#include	<stdlib.h>

class Exclusion;
class vtkPolyData;

class SurfGen {

    public:
	SurfGen() : un(50), vn(50) { srand(time(0)); }

	void	set_param( uint un, uint vn );
	bool	write_surface( const char *name, const BSplineSurface &s,
			      const Exclusion *excl = 0 );
	bool	write_real( const char *name, double maxerr,
			   const BSplineSurface &s );
	bool	write_controls( const char *name, const BSplineSurface &s );

    private:
	bool	make_points( vtkPolyData *spolygon, vtkPolyData *ppolygon, 
			    const BSplineSurface &s );
	bool	make_scalars( vtkPolyData *spolygon, const BSplineSurface &s,
			     int curvature );
	void	make_cells( vtkPolyData *spolygon, vtkPolyData *ppolygon,
			   const Exclusion *excl = 0 );
	Point	mess( const Point &source, double maxerr );

	uint	un;
	uint	vn;
};

#endif  // _SURGEN_HH_
