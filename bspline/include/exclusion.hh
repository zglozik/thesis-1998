
#ifndef	_EXCLUSION_HH_
#define	_EXCLUSION_HH_

#include	<debug.hh>
#include	<stdtypes.h>
#include	<list.hh>
#include	<set.hh>
#include	"polygon.hh"

typedef List<Polygon>	Polygons;

class Exclusion {

    public:
	Exclusion( const Polygons &polygons, uint n, uint m );
	Exclusion() { };

    public:
	void init( const Polygons &polygons, uint n, uint m );
	bool is_inside( uint i, uint j ) const;
	bool is_outside( uint i, uint j ) const;

	uint size() const;

    private:
	bool is_inside( const Polygons &polygons, double u, double v ) const;

	Set	points;
	uint	n, m;
};

istream &operator >> ( istream &in, Exclusion &e );

#endif // _EXCLUSION_HH_
