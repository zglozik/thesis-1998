
// ----------------------------------------------------------------------
//	BTree()
// ----------------------------------------------------------------------

template< class KEY, class ELEM, class KeyOfElem, class Compare >
inline BTree < KEY,ELEM,KeyOfElem,Compare>::BTree() 
: root(NULL), left_most(NULL), elements(0)
{
}

// ----------------------------------------------------------------------
//	BTree( const BTree &tree ) 
// ----------------------------------------------------------------------

template< class KEY, class ELEM, class KeyOfElem, class Compare >
inline BTree<KEY,ELEM,KeyOfElem,Compare>::BTree( const BTree &tree ) 
: root(NULL), left_most(NULL), elements(0)
{
	copy( tree.begin(), tree.end(), insert_iterator<BTree>(*this) );
}

// ----------------------------------------------------------------------
//	~BTree()
// ----------------------------------------------------------------------

template< class KEY, class ELEM, class KeyOfElem, class Compare >
inline BTree<KEY,ELEM,KeyOfElem,Compare>::~BTree()
{
	while( !is_empty() ) destroy( begin() );
}

// ----------------------------------------------------------------------
//	operator = ( const BTree &tree )
// ----------------------------------------------------------------------

template< class KEY, class ELEM, class KeyOfElem, class Compare >
BTree<KEY,ELEM,KeyOfElem,Compare> &
BTree<KEY,ELEM,KeyOfElem,Compare>::operator = ( const BTree &tree )
{
	if( &tree != this ) {
		while( !is_empty() ) destroy( begin() );
		copy( tree.begin(), tree.end(), 
		     insert_iterator<BTree>(*this) );
	}

	return *this;
}


// ----------------------------------------------------------------------
//	destroy( const ELEM &elem )
// ----------------------------------------------------------------------

template< class KEY, class ELEM, class KeyOfElem, class Compare >
bool BTree<KEY,ELEM,KeyOfElem,Compare>::destroy( const KEY &elem )
{
	iterator pos = find(elem);
	return destroy( pos );
}

// ----------------------------------------------------------------------
/*
#include	<stdlib.h>
#include	<time.h>

int main()
{
	typedef	BTree< int, int, identity<int>, is_less<int> >	IntBTree;

	IntBTree	t1, t2;

	srand( time(NULL) );
	for( int i = 0; i < 100; i++ ) {
		t1.insert( rand() % 10 );
	}

	IntBTree t3( t1 );
	while( t1.destroy(5) );
	cout << "end destroy\n";

	t2 = t1;

	cout << "---" << endl;

	for( IntBTree::iterator n = t1.begin(); n != t1.end(); n++ ) {
		cout << *n;
		cout << " left: " << t1.height( t1.left(n) );
		cout << " right: " << t1.height( t1.right(n) );
		cout << " diff: ";
		cout << t1.height(t1.left(n)) - t1.height(t1.right(n)) << endl;
	}

//	copy( t1.begin(), t1.end(), ostream_iterator<int>(cout, "\n") );
	cout << "---" << endl;

	return 0;
}
*/
