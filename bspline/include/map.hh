
#ifndef	_MAP_HH_
#define	_MAP_HH_

#include	"btree.hh"

template< class KEY, class ELEM >
class Map {

//	private classes

    private:
	class Pair {
	    public:
		Pair( const KEY &k, const ELEM &d ) : key(k), data(d) { }
		KEY	key;
		ELEM	data;
	};

	struct KeyOfValue {
		const KEY &operator()( const Pair &pair ) const
		{
			return pair.key;
		}
	};

	typedef	BTree< const KEY, Pair, KeyOfValue, is_less<KEY> > rep_type;

//	typedefs

	struct	Pair;
	struct	KeyOfValue;

    public:
	typedef	ELEM	value_type;
	typedef	KEY	key_type;
	
	class iterator {

	friend	Map;

	    protected:
		iterator( rep_type::iterator const &pos ) : pos(pos) { }

	    public:
		bool operator == ( const iterator &b ) const 
		{
			return pos == b.pos;
		}
		iterator &operator++()
		{
			++pos;
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
			--pos;
			return *this;
		}
		iterator operator--( int )
		{
			iterator tmp = *this;
			--*this;
			return tmp;
		}
		ELEM &operator * () { return (*pos).data; }

	    private:
		rep_type::iterator	pos;
	};

//	ctors & dtor

    public:
virtual	~Map() { };

//	opers

    public:
	iterator	begin()		{ return tree.begin(); }
	iterator	end()		{ return tree.end(); }
	const iterator	begin() const	{ return tree.begin(); }
	const iterator	end() const	{ return tree.end(); }

	Map	&operator = ( const Map &m ) { tree = m.tree; return *this; }

	ELEM		&operator[]( const KEY &key );
	const ELEM	&operator[]( const KEY &key ) const 
	{
		return ((Map*)this)->operator[]( key );
	}
	iterator find( const KEY &key ) const
	{ 
		return tree.find( key );
	}
	const KEY &key( iterator pos ) const
	{
		TEST_EXPR( pos != end() );
		return (*pos.pos).key;
	}
	
	bool destroy( const KEY &key );
	bool destroy( const iterator &pos )
	{
		return tree.destroy( pos.pos );
	}

	int	size() const { return tree.size(); }
	bool	is_empty() const { return tree.is_empty(); }

//	private data

	rep_type	tree;

};

//	include implementation
#include	"map.cc"

#endif // _MAP_HH_
