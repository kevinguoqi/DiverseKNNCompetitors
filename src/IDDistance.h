/*
 * IDDistance.h
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#ifndef IDDISTANCE_H_
#define IDDISTANCE_H_

class ID_Distance {

public:
	int id;
	double spatialDistance;
	double diverseDistance;

public:
	ID_Distance(int id, double spatialDistance, double diverseDistance);
	virtual ~ID_Distance();



public:
	inline bool operator<(const ID_Distance& rhs)
	{
	    return this->spatialDistance < rhs.spatialDistance;
	}
};

#endif /* IDDISTANCE_H_ */
