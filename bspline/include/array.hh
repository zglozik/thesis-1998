
#ifndef	_ARRAY_HH_
#define	_ARRAY_HH_

#include	"debug.hh"
#include	<iostream.h>

template< class ELEM >
class Array {

//	iterator

    public:
	typedef	ELEM	value_type;
	typedef	ELEM	*iterator;
 
//	ctors & dtors

    public:
	Array( int min_size = 0, int step_size = 10 );
	Array( const Array &array );
virtual	~Array();
	
//	opers

    public:
	iterator	begin()	{ return data; }
	iterator	end()	{ return data + size(); }
	const iterator	begin() const	{ return ((Array*)this)->begin(); }
	const iterator	end() const	{ return ((Array*)this)->end(); }
	ELEM &operator[](int n)	 
	{
		TEST_EXPR( 0 <= n && n < elements );
		return *(begin()+n);
	}
	const ELEM &operator[](int n) const { return (*(Array*)this)[n]; }
	ELEM	&front() { return *begin(); }
	ELEM	&back() { return *(end()-1); }
	const ELEM	&front() const { return *begin(); }
	const ELEM	&back() const { return *(end()-1); }

	Array	&operator = ( const Array &array );
	Array	&push_back( const ELEM &elem );

	void	resize( int new_size );
	bool	is_empty() const { return elements == 0; }
	int	size() const	{ return elements; }

//	private data

    private:
	ELEM	*data;
	int	elements;
	int	capacity;
	int	step_size;

};

template< class ELEM >
ostream &operator << ( ostream &o, const Array<ELEM> &a );

template< class ELEM >
istream &operator >> ( istream &i, Array<ELEM> &a );

//	include implementation
#include	"array.cc"

#endif				// _ARRAY_HH_
