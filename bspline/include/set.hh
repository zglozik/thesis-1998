
#ifndef	_SET_HH_
#define	_SET_HH_

#include	<debug.hh>
#include	<stdtypes.h>
#include	<iostream.h>

class Set {

    public:
	class iterator {

	friend	Set;
		iterator( const Set *set, int elem ) 
			: set(set), elem(elem) { }
	    public:
		iterator() : set(0), elem(0) { }

		bool operator == ( const iterator &b ) const
		{
			return set == b.set && elem == b.elem;
		}
		iterator &operator++()
		{
			elem++;
			while( elem <= set->max && !set->in(elem) ) elem++;
			if( elem > set->max ) elem = set->min-1;

			return *this;
		}
		iterator operator++( int )
		{
			iterator tmp = *this;
			++*this;
			return tmp;
		}
		iterator &operator--()
		{
			elem--;
			while( elem >= set->min && !set->in(elem) ) elem--;
			if( elem < set->min ) elem = set->min-1;

			return *this;
		}
		iterator operator--( int )
		{
			iterator tmp = *this;
			--*this;
			return tmp;
		}
		int operator *() const
		{
			return elem;
		}

	    private:
		const Set *set;
		int elem;
	}; // iterator

friend	iterator;
friend	ostream &operator << ( ostream &out, const Set &set );
friend	istream &operator >> ( istream &in, Set &set );

    public:
	Set( int min = 0, int max = 0 );
	Set( const Set &set );
virtual	~Set();

	iterator begin();
	iterator end();
	const iterator begin() const;
	const iterator end() const;

	Set &operator = ( const Set &set );

	void resize( int min, int max );

	bool in( int elem ) const;
	Set &put( int elem );
	Set &get( int elem );
	Set &empty();

	bool is_empty() const;
	uint size() const;
	
    private:
	void index( int elem, uint &i, uint &j ) const;

	uchar	*data;
	int	min;
	int	max;
};

ostream &operator << ( ostream &out, const Set &set );
istream &operator >> ( istream &in, Set &set );

#endif // _SET_HH_
