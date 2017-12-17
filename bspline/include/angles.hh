
#ifndef	_ANGLES_HH_
#define	_ANGLES_HH_

#include	<math.h>

inline double rad2deg( double a )
{
	return 180.0 * a / M_PI;
}

inline double deg2rad( double a )
{
	return M_PI * a / 180.0;
}

#endif				// _ANGLES_HH_
