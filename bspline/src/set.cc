
#include	"set.hh"
#include	<string.h>
#include	<funcs.hh>

//----------------------------------------------------------------------
//	iterator operators
//----------------------------------------------------------------------

Set::iterator Set::begin()
{
	int elem = min;
	while( elem <= max && !in(elem) ) elem++;
	
	return elem <= max ? iterator(this,elem) : iterator(this,min-1);
}

Set::iterator Set::end()
{
	return iterator(this,min-1);
}

const Set::iterator Set::begin() const
{
	return ((Set*)this)->begin();
}

const Set::iterator Set::end() const
{
	return ((Set*)this)->end();
}

//----------------------------------------------------------------------
//	Set::Set( int min, int max )
//----------------------------------------------------------------------

Set::Set( int min, int max )
{
	data = 0;
	resize( min, max );
}

//----------------------------------------------------------------------
//	Set::Set( const Set &set )
//----------------------------------------------------------------------

Set::Set( const Set &set ) : min(set.min), max(set.max)
{
	uint size = (max-min) / sizeof(uchar) + 1;
	data = new uchar [size];
	memcpy( data, set.data, sizeof(uchar)*size );
}

//----------------------------------------------------------------------
//	Set::~Set()
//----------------------------------------------------------------------

Set::~Set()
{
	delete [] data;
}

//----------------------------------------------------------------------
//	Set &Set::operator = ( const Set &set )
//----------------------------------------------------------------------

Set &Set::operator = ( const Set &set )
{
	if( this != &set ) {
		delete [] data;
		min = set.min;
		max = set.max;
		uint size = (max-min) / sizeof(uchar) + 1;
		data = new uchar [size];
		memcpy( data, set.data, sizeof(uchar)*size );
	}
	return *this;
}

//----------------------------------------------------------------------
//	void Set::resize( int min, int max )
//----------------------------------------------------------------------

void Set::resize( int min, int max )
{
	TEST_EXPR( min <= max );

	if( this->min == min && this->max == max ) return;

	delete [] data;
	this->min = min;
	this->max = max;
	uint size = (max-min) / sizeof(uchar) + 1;
	data = new uchar [size];
	memset( data, 0, sizeof(uchar)*size );
}

//----------------------------------------------------------------------
//	bool Set::in( int elem ) const
//----------------------------------------------------------------------

bool Set::in( int elem ) const
{
	uint i, j;

	index( elem, i, j );
	return !!(data[i] & (1u << j));
}

//----------------------------------------------------------------------
//	Set &Set::put( int elem )
//----------------------------------------------------------------------

Set &Set::put( int elem )
{
	uint i, j;

	index( elem, i, j );
	data[i] |= 1u << j;

	return *this;
}

//----------------------------------------------------------------------
//	Set &Set::get( int elem )
//----------------------------------------------------------------------

Set &Set::get( int elem )
{
	uint i, j;

	index( elem, i, j );
	data[i] &= ~(1u << j);

	return *this;
}

//----------------------------------------------------------------------
//	Set &Set::empty()
//----------------------------------------------------------------------

Set &Set::empty()
{
	uint size = (max-min) / sizeof(uchar) + 1;

	memset( data, 0, size );
	return *this;
}

//----------------------------------------------------------------------
//	bool Set::is_empty() const
//----------------------------------------------------------------------

bool Set::is_empty() const
{
	return size() == 0;
}

//----------------------------------------------------------------------
//	uint Set::size() const
//----------------------------------------------------------------------

uint Set::size() const
{
	uint size = (max-min) / sizeof(uchar) + 1;
	uint num = 0;

	for( uint i = 0; i < size; i++ ) {
		for( uint j = 0; j < sizeof(uchar); j++ ) {
			num += !!(data[i] & (1u << j));
		}
	}

	return num;
}

//----------------------------------------------------------------------
//	void Set::index( int elem, uint &i, uint &j ) const
//----------------------------------------------------------------------

void Set::index( int elem, uint &i, uint &j ) const
{
	TEST_EXPR( min <= elem && elem <= max );

	i = (elem-min) / sizeof(uchar);
	j = (elem-min) % sizeof(uchar);
}

//----------------------------------------------------------------------
//	ostream &operator << ( ostream &out, const Set &set )
//----------------------------------------------------------------------

ostream &operator << ( ostream &out, const Set &set )
{
	out << set.min << ' ' << set.max << endl;
	out << set.size() << endl;
	copy( set.begin(), set.end(), ostream_iterator<int>(out," ") );
	return out;
}

//----------------------------------------------------------------------
//	istream &operator >> ( istream &in, Set &set )
//----------------------------------------------------------------------

istream &operator >> ( istream &in, Set &set )
{
	int	min, max;
	uint	size;
	in >> min >> max;
	in >> size;
	set.resize( min, max );

	for( uint i = 0; i < size; i++ ) {
		int elem;
		in >> elem;
		set.put( elem );
	}
	return in;
}

//----------------------------------------------------------------------
/*
int main()
{
	Set	set;
	Set	set2( 10, 300 );

	cin >> set;
	set2 = set;

	cout << set2.in( 15 ) << endl;
	cout << set2.in( 20 ) << endl;
	cout << set2.in( 12 ) << endl;
	set2.get( 15 );
	set2.get( 20 );
	set2.put( 30 );

	cout << set2.size() << endl;
	cout << set.size() << endl;
	cout << set2 << endl;
	set2.empty();
	cout << set2.is_empty() << endl;
	cout << set.is_empty() << endl;

	return 0;
}
*/
