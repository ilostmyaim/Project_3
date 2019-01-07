#pragma once

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <list>
#include "hash.h"

class Cluster {
private:
	int _cluster_id; //id of a cluster
	unsigned int _centroid_id; //id of centroid
	vector_t _centroid;
	std::vector<item_t> _items;
	std::vector<double> _ai_values; //used in silhouette
	std::vector<double> _cj; //k-dimensions, for k coins
	double _ai;
public:
	Cluster(int cluster_id, item_t item);
	void addItem(item_t item);
	bool removeItem(unsigned int item_id);
	vector_t getCentroid();
	void setCentroidValue(int index,double value);
	void setCentroid(vector_t vec);
	item_t getItem(int index);
	double getTotalItems();
	int getID();
	int getCentroidSize();
	unsigned int getCentroidID();
	void add(const item_t &item);
	void calculateFinal(int totalItems);
	void insertItem(item_t item);
	double computeAvgDistance(Metric metric); //computes a(i)
};

typedef struct {
	std::string text;
	int user_id;
	double score;
} tweet;


class User {
private:
	int _user_id;
	std::vector<tweet> _tweets;
public:

}