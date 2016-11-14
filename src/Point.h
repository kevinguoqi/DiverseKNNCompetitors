/*
 * Point.h
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#ifndef POINT_H_
#define POINT_H_

#include <vector>

using namespace std;

class Point {

public:

	vector<double> attributes;


public:
	Point();
	virtual ~Point();

	Point( const Point& p2 )
	{
		this->attributes.reserve( p2.attributes.size() );
		copy( p2.attributes.begin(), p2.attributes.end(), back_inserter( this->attributes ));
	}

	/** Move constructor */
	Point (Point&& p2) noexcept  /* noexcept needed to enable optimizations in containers */
	{
		this->attributes.reserve( p2.attributes.size() );
		copy( p2.attributes.begin(), p2.attributes.end(), back_inserter( this->attributes ));
		p2.attributes.clear();
	}


public:
//	void setAttributes( vector<double>& attributes );


public:
	inline Point& operator=( const Point& p2 )
	{
		this->attributes.reserve( p2.attributes.size() );
		copy( p2.attributes.begin(), p2.attributes.end(), back_inserter( this->attributes ));
	    return *this;
	}

	/** Move assignment operator */
	Point& operator= (Point&& p2) noexcept
	{
		this->attributes.clear();
		this->attributes.reserve( p2.attributes.size() );
		copy( p2.attributes.begin(), p2.attributes.end(), back_inserter( this->attributes ));
		p2.attributes.clear();
		return *this;
	}

};

#endif /* POINT_H_ */
