
#include	"funcs.hh"

template< class ELEM >
Array<ELEM>::Array( int min_size, int step_size )
: elements(min_size), capacity(min_size), step_size(step_size)
{
	TEST_EXPR( min_size >= 0 && step_size > 0 );

	data = elements > 0 ? new ELEM [elements] : (ELEM *) NULL;
}

template< class ELEM >
Array<ELEM>::Array( const Array<ELEM> &array )
: elements(array.size()), capacity(array.size()), step_size(array.step_size)
{
	data = elements > 0 ? new ELEM [elements] : (ELEM *) NULL;

	copy( array.begin(), array.end(), begin() );
}

template< class ELEM >
Array<ELEM>::~Array()
{
	delete [] data;
}

template< class ELEM >
Array<ELEM> &Array<ELEM>::operator = ( const Array<ELEM> &array )
{
	if( &array != this ) {
		if( capacity != array.size() ) {
			delete [] data;
			capacity = array.size();
			data = capacity > 0 ? new ELEM [capacity] : 0;
		}
		elements = capacity;
		copy( array.begin(), array.end(), begin() );
	}
	
	return *this;
}

template< class ELEM >
Array<ELEM> &Array<ELEM>::push_back( const ELEM &elem )
{
	if( capacity == elements ) {
		int old_size = elements;
		resize( capacity + step_size );
		elements = old_size;
	}
	data[elements++] = elem;
	
	return *this;
}

template< class ELEM >
void Array<ELEM>::resize( int new_size )
{
	TEST_EXPR( new_size >= 0 );

	if( new_size != capacity ) {
		ELEM	*new_data = new ELEM [new_size];

		for( int i = 0; i != min(elements, new_size); i++ ) {
			new_data[i] = data[i];
		}
		delete [] data;
		data = new_data;
		elements = new_size;
		capacity = new_size;
	}
}

template< class ELEM >
ostream &operator << ( ostream &o, const Array<ELEM> &a )
{
	o << a.size() << endl;
	copy( a.begin(), a.end(), ostream_iterator<ELEM>(o, "\n") );
	return o;
}

template< class ELEM >
istream &operator >> ( istream &i, Array<ELEM> &a )
{
	int n = 0;

	i >> n;
	a.resize( n );
	for( int k = 0; k < n; k++ ) {
		i >> a[k];
	}
	return i;
}

//----------------------------------------------------------------------
/*
#include	<iostream.h>

int main()
{
	Array<int>		array( 5 ), array2;
	Array<int>::iterator	a;
	int			i;

	for( i = 0, a = array.begin(); a != array.end(); a++, i++ ) {
		*a = i;
	}
	for( i = 20; i != 31; i++ ) {
		array.push_back( i );
	}
	array[2] = 7;
	array2 = array;
	Array<int>	array3(0, 10);
	array3 = array2;
	array3.resize( 4 );

	copy( array.begin(), array.end(), ostream_iterator<int>(cout, "\n") );
	cout << "---\n";
	copy( array3.begin(), array3.end(), 
	     ostream_iterator<int>(cout, "\n") );

	return 0;
}
*/
