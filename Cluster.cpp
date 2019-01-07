#include<iostream>
#include<vector>
#include<string>
#include<cstring>
#include<string.h>
#include <fstream>
#include <sstream>
#include<math.h>
#include <algorithm>
#include <chrono>
#include "Cluster.h"
#include "hash.h"

using namespace std;

Cluster::Cluster(int cluster_id, item_t item)
{
	this->_cluster_id = cluster_id;
	this->_centroid_id = item.id;
	this->_ai=0;
	for(int i=0; i<(int)item.vec.size(); i++){
		_centroid.push_back(item.vec[i]);	
	}
	//cout << "Centorid size " << _centroid.size() << endl;

}


bool Cluster::removeItem(unsigned int item_id)
{
	for(int i=0;i<(int)_items.size();i++){
		if(_items[i].id == item_id){
			_items.erase(_items.begin() + i);
			return true;
		}
	}
	return false;
}


vector_t Cluster::getCentroid()
{
	return _centroid;
}

void Cluster::setCentroidValue(int index, double value)
{
	_centroid[index] = value;
}

void Cluster::setCentroid(vector_t vec)
{
	_centroid = vec;
}

item_t Cluster::getItem(int index)
{
	return _items[index];
}



double Cluster::getTotalItems()
{
	return _items.size();
}


int Cluster::getID()
{
	return this->_cluster_id;
}

int Cluster::getCentroidSize()
{
	return (int)this->_centroid.size();
}

unsigned int Cluster::getCentroidID()
{
	return this->_centroid_id;
}

void Cluster::add(const item_t &item)
{
	int i = 0;
	for(i=0;i<_centroid.size();i++){
		_centroid[i] += item.vec[i];
	}
}

void Cluster::calculateFinal(int totalItems)
{
	unsigned int i;
	for(i=0;i<(int)_centroid.size();i++){
		_centroid[i] /= double(totalItems);
	}
}

/*insert item to cluster*/
void Cluster::insertItem(item_t item)
{
	this->_items.push_back(item);
}

double Cluster::computeAvgDistance(Metric metric) //computes a(i) for a cluster
{
	int i,j;
	double avgDist_a = 0;
	double totalAvg_a = 0;
	double ai=0;
	if(metric == euclidean) { 
		for(i=0;i<_items.size();i++) {
			avgDist_a = 0;
			for(j=0;j<_items.size();j++) {
				avgDist_a += euclideanNorm(_items[i].vec, _items[j].vec);
			}
			totalAvg_a = avgDist_a / double(_items.size());
			_ai_values.push_back(totalAvg_a);
		}
	}
	else{
		for(i=0;i<_items.size();i++) {
			avgDist_a = 0;
			for(j=0;j<_items.size();j++) {
				avgDist_a += (1-cosineSimilarity(_items[i].vec, _items[j].vec));
			}
			totalAvg_a = avgDist_a / double(_items.size());
			_ai_values.push_back(totalAvg_a);
		}
	}
	
	
	for(i=0;i<_ai_values.size();i++){
		_ai += _ai_values[i];
	}

	return (_ai / double(_ai_values.size()));
}
