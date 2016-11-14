/*
 * Dataset.cpp
 *
 *  Created on: Oct 4, 2016
 *      Author: zhengyuxin
 */

#include "Dataset.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

Dataset::Dataset(int dim) {
	this->dim = dim;
}

Dataset::~Dataset() {
	// TODO Auto-generated destructor stub
}


void Dataset::setDataSet( vector<Point>& dataset )
{
//	this->dataset.insert( this->dataset.begin(), dataset.begin(), dataset.end() );

	this->dataset.reserve( dataset.size() );
	copy( dataset.begin(), dataset.end(), back_inserter( this->dataset ));

}

void Dataset::setQuerySet( vector<Point>& queryset )
{
//	this->queryset.insert( this->queryset.begin(), queryset.begin(), queryset.end() );

	this->queryset.reserve( queryset.size() );
	copy( queryset.begin(), queryset.end(), back_inserter( this->queryset ));

}


void Dataset::setSpatialIndex( vector<int>& spatialIndex )
{
//	this->spatialIndex.insert( this->spatialIndex.begin(), spatialIndex.begin(), spatialIndex.end() );
	this->spatialIndex.reserve( spatialIndex.size() );
	copy( spatialIndex.begin(), spatialIndex.end(), back_inserter( this->spatialIndex ));
}

void Dataset::setDiverseIndex( vector<int>& diverseIndex )
{
//	this->diverseIndex.insert( this->diverseIndex.begin(), diverseIndex.begin(), diverseIndex.end() );
	this->diverseIndex.reserve( diverseIndex.size() );
	copy( diverseIndex.begin(), diverseIndex.end(), back_inserter( this->diverseIndex ));
}


// read data set
void Dataset::readDataSet( string& inputFile, int num = INT_MAX, char sep=',', bool skipFirstColumn = false )
{
	if( !skipFirstColumn ) {
		readFile( inputFile, dataset, num, sep );
	}
	else {
		readFileSkipFirstColumn( inputFile, dataset, num, sep );
	}

	cout << "data set size: " << dataset.size() << endl;
}


// read query set
void Dataset::readQuerySet( string& inputFile, int num = INT_MAX, char sep=',', bool skipFirstColumn = false )
{
	if( !skipFirstColumn ) {
		readFile( inputFile, queryset, num, sep );
	}
	else {
		readFileSkipFirstColumn( inputFile, queryset, num, sep );
	}
	cout << "query set size: " << queryset.size() << endl;
}



// read data from a file and add to an array or matrix
void Dataset::readFile( string& inputFile, vector<Point> &set, int num, char sep )
{
	ifstream input( inputFile );
	if ( input.is_open() )
	{
		string csvLine;
		int lineCount = 0;
		// read every line from the stream
		while( std::getline(input, csvLine) && lineCount < num )
		{
			Point p;
			istringstream csvStream(csvLine);
			string csvElement;
			for (int d=0; d < dim; d++)
			{
				// split the line into tokens
				std::getline(csvStream, csvElement, sep);
				p.attributes.push_back( atof( csvElement.c_str() ) );
			}
			set.push_back(p);
			lineCount++;
		}
	}
}



// read data from a file and add to an array or matrix
void Dataset::readFileSkipFirstColumn( string& inputFile, vector<Point> &set, int num, char sep )
{
	ifstream input( inputFile );
	if ( input.is_open() )
	{
		string csvLine;
		int lineCount = 0;
		// read every line from the stream
		while( std::getline(input, csvLine) && lineCount < num )
		{
			Point p;
			istringstream csvStream(csvLine);
			string csvElement;

			// the first column is id, skip
			std::getline(csvStream, csvElement, sep);

			for (int d=0; d < dim; d++)
			{
				// split the line into tokens
				std::getline(csvStream, csvElement, sep);
				p.attributes.push_back( atof( csvElement.c_str() ) );
			}
			set.push_back(p);
			lineCount++;
		}
	}
}






