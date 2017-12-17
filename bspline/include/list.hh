
#ifndef	_LIST_HH_
#define	_LIST_HH_

#include	"debug.hh"
#include	<iostream.h>

template< class ELEM >
class List {

    public:
	typedef	ELEM	value_type;

//	implementation class

    private:
	struct Node {
		Node( const ELEM &data, Node *prev, Node *next ) 
			: data( data ), prev( prev ), next( next ) { }
		ELEM	data;
		Node	*prev, *next;
	};


//	iterator

    public:
	class iterator {
	    friend List;

		iterator( Node *x ) : node( x ) { TEST_EXPR( x != NULL ); }
	    public:
		iterator() : node(0) { }

		bool	 operator == ( const iterator &b ) const
		{
			return node == b.node;
		}
		iterator &operator++()
		{
			node = node->next;
			return *this;
		}
		iterator operator++( int )
		{
			iterator	tmp = *this;
			++*this;
			return tmp;
		}
		iterator &operator--()
		{
			node = node->prev;
			return *this;
		}
		iterator operator--( int )
		{
			iterator	tmp = *this;
			--*this;
			return tmp;
		}
		ELEM	&operator *() const 
		{
			return node->data;
		}

	    private:
		Node	*node;

	}; // iterator

//	ctors & dtor

    public:
	List();
	List( const List &list );
virtual	~List();

//	opers

    public:
	iterator	begin() { return head->next; }
	iterator	end()	{ return head; }
	const iterator	begin() const	{ return ((List*)this)->begin(); }
	const iterator	end() const	{ return ((List*)this)->end(); }

	List	&operator = ( const List &list );
	List	&push_front( const ELEM &elem );
	List	&push_back( const ELEM &elem );
	List	&insert( const ELEM &elem, const iterator &pos )
	{	
		Node *new_node = new Node( elem, pos.node->prev, pos.node );
		pos.node->prev->next = new_node;
		pos.node->prev = new_node;
		++elements;
		return *this;
	}

	bool	destroy( const iterator &pos )
	{
		if( pos == end() ) return false;

		Node *tmp_node = pos.node;
		tmp_node->prev->next = tmp_node->next;
		tmp_node->next->prev = tmp_node->prev;
		delete tmp_node;
		--elements;
		return true;
	}
	ELEM	&front() { return *begin(); }
	ELEM	&back() { return head->prev->data; }
	const ELEM	&front() const { return *begin(); }
	const ELEM	&back() const { return head->prev->data; }
	int	size() const { return elements; }
	bool	is_empty() const { return elements == 0; }

//	private data

    private:
	Node	*head;
	int	elements;

}; // List

template< class ELEM >
ostream &operator << ( ostream &o, const List<ELEM> &a );

template< class ELEM >
istream &operator >> ( istream &i, List<ELEM> &a );

//	include implementation
#include	"list.cc"

#endif				// _LIST_HH_
