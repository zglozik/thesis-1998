
#ifndef	_BTREE_HH_
#define	_BTREE_HH_

#include	"debug.hh"
#include	"funcs.hh"
#include	"pair.hh"

template< class KEY, class ELEM, class KeyOfElem, class Compare >
class BTree {

    public:
	struct Node {
		ELEM	data;
		int	state;	// it must be one of the values -1, 0, +1
		Node	*parent;
		Node	*left, *right;

		Node( const ELEM &data, Node *parent = NULL )
			: data(data), state(0), 
			parent(parent), left(0), right(0) { }
		Node *left_most()
		{
			Node *node = this;

			while( node->left ) node = node->left;
			return node;
		}
		Node *right_most()
		{
			Node *node = this;

			while( node->right ) node = node->right;
			return node;
		}
	};

//	typedefs

    public:
	typedef	ELEM	value_type;
	typedef	KEY	key_type;

	class	iterator;
	friend	iterator;

	class iterator {

	friend	Node;
	friend	BTree;

		iterator( Node *node ) : node(node) { }
	    public:
		iterator() { }

		bool operator == ( const iterator &b ) const
		{
			return node == b.node;
		}

		iterator &operator ++ ()
		{
			TEST_EXPR( node != NULL );
			if( node->right ) {
				node = node->right->left_most();
			} else {
				Node *prev = node;
				node = node->parent;

				while( node && node->right == prev ) {
					prev = node;
					node = node->parent;
				}
			}
			return *this;
		}

		iterator operator ++ (int)
		{
			iterator tmp = *this;
			++*this;
			return tmp;
		}

		iterator &operator -- ()
		{
			TEST_EXPR( node != NULL );
			if( node->left ) {
				node = node->left->right_most();
			} else {
				Node *prev = node;
				node = node->parent;

				while( node && node->left == prev ) {
					prev = node;
					node = node->parent;
				}
			}
			return *this;
		}

		iterator operator -- (int)
		{
			iterator tmp = *this;
			--*this;
			return tmp;
		}

		ELEM	&operator * () const
		{
			TEST_EXPR( node != NULL );
			return node->data;
		}

	    private:
		Node	*node;
	};

//	ctors & dtor

    public:
	BTree();
	BTree( const BTree &BTree );
virtual	~BTree();

//	opers

    public:
	BTree	&operator = ( const BTree &tree );

	iterator	begin()		{ return left_most; }
	iterator	end()		{ return NULL; }
	const iterator	begin() const	{ return ((BTree*)this)->begin(); }
	const iterator	end() const	{ return ((BTree*)this)->end(); }

	ELEM	&front() { return *begin(); }
	ELEM	&back() { return root->right_most()->data; }
	const ELEM	&front() const { return *begin(); }
	const ELEM	&back() const { return root->right_most()->data; }

	bool		destroy( const KEY &elem );

	int height( const iterator &pos ) const
	{
		return pos == end() ? 0 
			: 1 + max( height(left(pos)), height(right(pos)) );
	}

	iterator left( const iterator &pos ) const
	{
		TEST_EXPR( pos.node != NULL );
		return pos.node->left;
	}

	iterator right( const iterator &pos ) const
	{
		TEST_EXPR( pos.node != NULL );
		return pos.node->right;
	}

// ----------------------------------------------------------------------
//	insert( const ELEM &elem )
// ----------------------------------------------------------------------

	iterator insert( const ELEM &elem )
	{
		Node	*parent, **branch;

		((BTree*)this)->_find_place( parent, branch, 
					    key_of_value(elem) );
		
		Node *new_node = new Node( elem, parent );
		*branch = new_node;
		elements++;

		Node *prev = new_node;
		int h = +1;
		while( parent && h != 0 ) {
			int c = change( parent, prev );
			
			if( parent->state - c == 0 ) {
				h += state_changed( parent, c );
			} else {
				parent->state += c;
				if( parent->state == 0 ) h = 0;
				prev = parent;
				parent = parent->parent;
			}
		}
		left_most = root ? root->left_most() : 0;

		return new_node;
	}

// ----------------------------------------------------------------------
//	destroy( const iterator &pos )
// ----------------------------------------------------------------------

	bool destroy( const iterator &pos )
	{
		if( !pos.node )	return false;

		Node *node = pos.node;
		if( node->right ) {
			Node *next = node->right->left_most();
			node->data = next->data;
			node = next;
		} 

		Node *parent = node->parent;
		int	c = 0;	// from which direction we came
		if( parent ) c = -change( parent, node ); 

		Node **branch = parent_ref( node );
		*branch = node->left ? node->left : node->right;
		if( *branch ) (*branch)->parent = node->parent;

		int h = -1;
		while( parent && h != 0 ) {
			if( parent->state - c == 0 ) {
				h = state_changed( parent, c );
			} else {
				parent->state += c;
				if( parent->state != 0 ) h = 0;
			}
			if( parent->parent ) {
				c = -change( parent->parent, parent );
			}
			parent = parent->parent;
		}

		--elements;
		delete node;
		left_most = root ? root->left_most() : 0;

		return true;
	}

// ----------------------------------------------------------------------
//	find( const KEY &key ) const
// ----------------------------------------------------------------------

	iterator find( const KEY &key ) const
	{
		Node	*parent, **branch;

		_find_elem( parent, branch, key );
		return *branch;
	}

	int	size() const		{ return elements; }
	bool	is_empty() const	{ return elements == 0; }
	
//	private data

    private:

	Node **parent_ref( Node *node )
	{
		TEST_EXPR( node != NULL );
		Node	*parent = node->parent;
		Node	**branch;

		if( !parent ) {
			branch = &root;
		} else if( parent->left == node ) {
			branch = &parent->left;
		} else {
			branch = &parent->right;
		}
		return branch;
	}

	int change( Node *parent, Node *child ) const
	{
		return (parent->right == child) - (parent->left == child);
	}

// ----------------------------------------------------------------------
//	_find_place( Node *&parent, Node **&branch, const KEY &key ) const
// ----------------------------------------------------------------------

	void _find_place( Node *&parent, Node **&branch, const KEY &key ) const
	{
		parent = 0;
		Node *const *current = &root;

		while( *current ) {
			parent = *current;
			if( compare(key_of_value((*current)->data), key) ) {
				current = &(*current)->right;
			} else {
				current = &(*current)->left;
			}
		}
		branch = (Node **) current;
	};

// ----------------------------------------------------------------------
//	_find_elem( Node *parent, Node *const &branch, const KEY &key )
// ----------------------------------------------------------------------

	bool _find_elem( Node *&parent, Node **&branch,	const KEY &key ) const
	{
		Node	*const *current = &root;

		while( *current && key_of_value((*current)->data) != key ) {
			parent = *current;
			if( compare(key_of_value((*current)->data), key) ) {
				current = &(*current)->right;
			} else {
				current = &(*current)->left;
			}
		}
		branch = (Node **) current;
		return !*current;
	}
	
// ----------------------------------------------------------------------
//	state_changed( Node *&node, int diff )
// ----------------------------------------------------------------------

	int state_changed( Node *&node, int diff )
	{
		int	new_state = node->state + diff;
		Node	**branch;
		Node	*p1, *p2, *p3;
		Node	*parent = node->parent;
		int	result;

		branch = parent_ref( node );
		switch( new_state ) {
		    case -2:
			p3 = node;
			switch( node->left->state ) {
			    case 1:
				p1 = node->left;
				p2 = node->left->right;
				p1->parent = p2;
				p1->right = p2->left;
				if( p1->right ) p1->right->parent = p1;
				p3->parent = p2;
				p3->left = p2->right;
				if( p3->left ) p3->left->parent = p3;

				p1->state = p2->state > 0 ? -1 : 0;
				p3->state = p2->state < 0 ? +1 : 0;
				p2->state = 0;
				result = -1;
				break;

			    case 0:
			    case -1:
				p1 = node->left->left;
				p2 = node->left;
				p1->parent = p2;
				p3->parent = p2;
				p3->left = p2->right;
				if( p3->left ) p3->left->parent = p3;
				
				if( p2->state == 0 ) {
					p3->state = -1;
					p2->state = +1;
					result = 0;
				} else {
					p2->state = 0;
					p3->state = 0;
					result = -1;
				}
				break;

			    default:
				ERROR( "internal" );
			}
			break;

		    case 2:
			p1 = node;
			switch( node->right->state ) {
			    case -1:
				p2 = node->right->left;
				p3 = node->right;
				p1->parent = p2;
				p1->right = p2->left;
				if( p1->right ) p1->right->parent = p1;
				p3->parent = p2;
				p3->left = p2->right;
				if( p3->left ) p3->left->parent = p3;
				
				p1->state = p2->state > 0 ? -1 : 0;
				p3->state = p2->state < 0 ? +1 : 0;
				p2->state = 0;
				result = -1;
				break;

			    case 0:
			    case 1:
				p2 = node->right;
				p3 = node->right->right;
				p1->parent = p2;
				p1->right = p2->left;
				if( p1->right ) p1->right->parent = p1;
				p3->parent = p2;

				if( p2->state == 0 ) {
					p1->state = +1;
					p2->state = -1;
					result = 0;
				} else {
					p1->state = 0;
					p2->state = 0;
					result = -1;
				}
				break;

			    default:
				ERROR( "internal" );
			}
			break;

		    default:
			ERROR( "internal" );
		}
		node = p2;
		*branch = p2;
		p2->parent = parent;
		p2->left = p1;
		p2->right = p3;

		return result;		// change of height
	}

	Node	*root;
	Node	*left_most;
	int	elements;

	KeyOfElem	key_of_value;
	Compare		compare;
	
};

//	include implementation
#include	"btree.cc"

#endif // _BTREE_HH_
