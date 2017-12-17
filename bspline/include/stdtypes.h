
#ifndef	STDTYPES_H
#define	STDTYPES_H

typedef	unsigned	uint;
typedef	unsigned long	ulong;
typedef	unsigned char	uchar;

#ifdef	NEED_BOOL
typedef enum { false, true }	bool;
#endif /* NEED_BOOL */

#define	ALLOC( size, type )	((type *) Alloc( size*sizeof(type)) )
#define	FREE( p )		(free(p), p = NULL)

#define	ARRAY_SIZE( a )		(sizeof(a)/sizeof(a[0]))

#endif /* STDTYPES_H */
