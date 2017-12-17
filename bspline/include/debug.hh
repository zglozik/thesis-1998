
#ifndef	DEBUG_HH
#define	DEBUG_HH

#include	<stdlib.h>
#include	<stdio.h>

#ifdef	NEED_BOOL

typedef	int	bool;

#define	true	1
#define	false	0

#endif // NEED_BOOL

#define	ERROR( x )	{ fprintf( stderr, "error in file %s(%d): %s\n",\
				 __FILE__, __LINE__, x );		\
			  exit( 1 ); }

#ifdef	DEBUG

#include	<stdio.h>
#include	<assert.h>

#define	LOG( text )	fprintf( stderr, "function invoked: %s()"	\
				" in file: %s(%d)\n",			\
				text, __FILE__, __LINE__ )

#define	TEST_EXPR( expr )	assert( expr )

#else

#define	LOG( text )
#define	TEST_EXPR( expr )

#endif /* debug */

#endif /* DEBUG_HH */
