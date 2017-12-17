
#ifndef	_PAIR_HH_
#define	_PAIR_HH_

#include	"debug.hh"
#include	"funcs.hh"

template< class T1, class T2 >
struct Pair {
	typedef	T1	first_type;	
	typedef	T2	second_type;

	Pair() { }
	Pair( const T1 &a, const T2 &b ) : a(a), b(b) { }
virtual ~Pair() { }

	T1 &first()	{ return a; }
	T2 &second()	{ return b; }
	const T1 &first() const		{ return a; }
	const T2 &second() const	{ return b; }
	
    private:
	T1	a;
	T2	b;
};

template< class T1, class T2 >
inline Pair< T1, T2 > make_pair( const T1 &a, const T2 &b )
{
	return Pair<T1, T2>( a, b );
}

template< class T1, class T2 >
inline bool operator == ( const Pair<T1,T2> &p1, const Pair<T1,T2> &p2 )
{
	return p1.first() == p2.first() && p1.second() == p2.second();
}

template< class T1, class T2 >
inline bool operator < ( const Pair<T1,T2> &p1, const Pair<T1,T2> &p2 )
{
	return p1.first() < p2.first() || 
		p1.first() == p2.first() && p1.second() < p2.second();
}

#endif // _PAIR_HH_
