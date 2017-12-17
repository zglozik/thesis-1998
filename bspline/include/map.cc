

template< class KEY, class ELEM >
ELEM &Map<KEY,ELEM>::operator[]( const KEY &key )
{
	iterator pos = tree.find( key );
	if( pos == end() ) {
		pos = tree.insert( Pair(key, ELEM()) );
	}

	return *pos;
}

template< class KEY, class ELEM >
bool Map<KEY,ELEM>::destroy( const KEY &key )
{
	return tree.destroy( key );
}

// ----------------------------------------------------------------------
/*
int main()
{
	Map< int, double >	floats;

	for( int i = 0; i < 30; i++ ) {
		floats[i] = i;
	}
	floats[4] = 4.0;
	floats[0] = 0.;
	floats.destroy( 5 );
	floats.destroy( 10 );

	copy( floats.begin(), floats.end(), 
	     ostream_iterator<double>(cout, "\n") );

	return 0;
}
*/
