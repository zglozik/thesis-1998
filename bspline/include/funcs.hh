
#ifndef	_FUNCS_HH_
#define	_FUNCS_HH_

#define	DOUBLE_PRECISION	(1.e-14)

template< class ELEM >
inline ELEM square( const ELEM &elem )
{
	return elem * elem;
}

inline double abs( double x )
{
	return x < 0 ? -x : x;
}

inline bool equal( const double &a, const double &b )
{
	return abs( a - b ) < DOUBLE_PRECISION;
}

inline bool less( const double &a, const double &b )
{
	return a < b && !equal(a,b);
}

template< class T >
inline bool operator != ( const T &a, const T &b )
{
	return !(a == b);
}

template< class Container >
inline Container &empty( Container &c )
{
	while( !c.is_empty() ) c.destroy( c.begin() );
	return c;
}

template< class It1, class It2 >
inline It2 copy( It1 begin, It1 end, It2 dest )
{
	while( begin != end ) *dest++ = *begin++;

	return dest;
}

template< class It1, class Pred >
inline It1 find_if( It1 begin, It1 end, const Pred &pred )
{
	while( begin != end && !pred(*begin) ) ++begin;

	return begin;
}

template< class It1, class ELEM >
inline It1 find( It1 begin, It1 end, const ELEM &value )
{
	while( begin != end && *begin != value ) ++begin;

	return begin;
}

template< class It1, class ELEM >
inline It1 find_ref( It1 begin, It1 end, const ELEM &value )
{
	while( begin != end && &*begin != &value ) ++begin;

	return begin;
}

template< class Input, class Output, class Op >
inline Output transform( Input begin, Input end, Output result, const Op &op )
{
	while( begin != end ) *result++ = op(*begin++);

	return result;
}

template< class T >
inline const T &max( const T &a, const T &b )
{
	return a < b ? b : a;
}

template< class T >
inline const T &min( const T &a, const T &b )
{
	return a < b ? a : b;
}

template< class Container >
class insert_iterator {
    public:
	insert_iterator( Container &x )  : x( x ) { }

	insert_iterator &operator = ( const typename Container::value_type
				     &value )
	{
		x.insert( value );
		return *this;
	}
	insert_iterator	operator * ()		{ return *this; }
	insert_iterator	operator ++ ()		{ return *this; }
	insert_iterator	operator ++ (int)	{ return *this; }
    private:
	Container	&x;
};

template< class Container >
class back_insert_iterator {
    public:
	back_insert_iterator( Container &x ) : x( x ) { }

	back_insert_iterator &operator=( const typename Container::value_type
					&value )
	{
		x.push_back( value );
		return *this;
	}
	back_insert_iterator	operator * ()	{ return *this; }
	back_insert_iterator	operator ++ ()	{ return *this; }
	back_insert_iterator	operator ++ (int)	{ return *this; }
    private:
	Container	&x;
};

template< class ELEM >
class ostream_iterator {

    public:
	ostream_iterator( ostream &s ) : stream(s), string(NULL) { }
	ostream_iterator( ostream &s, char *str ) : stream(s), string(str) { }
	ostream_iterator &operator = ( const ELEM &value )
	{
		stream << value;
		if( string ) stream << string;
		return *this;
	}
	ostream_iterator operator * () { return *this; }
	ostream_iterator operator ++ ()	{ return *this; }
	ostream_iterator operator ++ (int) { return *this; }
    private:
	ostream	&stream;
	char	*string;
};

template< class Arg1, class Arg2, class Return >
struct binary_operation {
	typedef	Arg1	arg1st_type;
	typedef	Arg2	arg2nd_type;
	typedef	Return	return_type;

virtual	return_type operator()( const arg1st_type &, const arg2nd_type & )
	const = 0;
};

template< class Arg, class Return >
struct unary_operation {
	typedef	Arg	arg1st_type;
	typedef	Return	return_type;

virtual	return_type operator()( const arg1st_type & ) const = 0;
};

template< class ELEM >
struct identity : public unary_operation< ELEM, ELEM & > {
	ELEM &operator ()( const ELEM &a ) const { return (ELEM &)a; }
};

template< class ELEM >
struct is_less : public binary_operation< ELEM, ELEM, bool > {
	bool	operator()( const ELEM &a, const ELEM &b ) const
	{
		return a < b;
	}
};

template< class ELEM >
struct is_equal : public binary_operation< ELEM, ELEM, bool > {
	bool	operator()( const ELEM &a, const ELEM &b ) const
	{
		return a == b;
	}
};

#endif				// _FUNCS_HH_
