#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include <random>
#include <list>
#include <ctime>
#include <chrono>
#include <iomanip>
#include "LSH.h"
#include "hash.h"


using namespace std;
using namespace std::chrono;

using std::default_random_engine;
int seed = time(0);

default_random_engine generator(seed);

Hash::Hash(int b,int k)
{
	this->_tableSize = b;
	//cout << "Hashtable size is: " << _tableSize << endl;
	this->_w = 2; 
	for(int i = 0;i<k;i++){
		this->vec_t.push_back(random_offset());
		this->vec_v.push_back(random_vector());
	}
	_hashTable = new list<item_t>[_tableSize];

}

void Hash::insertItem(item_t item,unsigned int hashValue)
{
	_hashTable[hashValue].push_back(item);
}

/*range search*/
vector<item_t> Hash::traverseBucket(vector_t q, int cluster_id, long int hashValue,double R, double C=1,Metric metric=euclidean)
{
	double distance = 0;
	vector<item_t> items;
	
	for(auto &x: _hashTable[hashValue]){
		/*calculate Eucledian deistance*/
		if(metric == euclidean) { 
			distance = euclideanNorm(q,x.vec);
			if( distance < C*R){
				if(x.cluster_id == -1){
					x.cluster_id = cluster_id;
					x.distance = distance;
					items.push_back(x);
				}
				else if((x.cluster_id != -1) && (cluster_id != x.cluster_id)){ //if more than two centroid exist in bucket
					if(distance < x.distance){
						x.cluster_id = cluster_id;
						x.distance = distance;
						items.push_back(x);
					}
				}
			}
		}
		else if(metric == cosine) {
			distance = 1 -cosineSimilarity(q,x.vec);
			if( distance < C*R){ 
				if(x.cluster_id = -1){
					x.cluster_id = cluster_id;
					items.push_back(x);
				}
				else if((x.cluster_id != -1) && (cluster_id != x.cluster_id)){ //if more than two centroid exist in bucket
					if(distance < x.distance){
						x.cluster_id = cluster_id;
						x.distance = distance;
						items.push_back(x);
					}
				}
			}
		}	
	}
	return items;
}
/*Nearest Neighbor serach*/
double Hash::nearestNeighborTraverse(vector_t q, long int hashValue, int L, Metric metric)
{
	double distance = 0, app_distance = 30000000;
	double min_distance = 30000000;
	int count_retrieved = 0; // check how many items have been retrieved
	int nearest_id,app_id = 0; //NN and approximate NN id
	vector_t nearest_vec;
	double approximation_factor = 0.8;//approximation factor for approximate NN
	auto start_appNN = high_resolution_clock::now();//start stopwatch for approximate NN
	auto start_NN = high_resolution_clock::now();//start stopwatch for NN
	cout.precision(9);
	
	//check if bucket is empty first
	if(_hashTable[hashValue].empty())
		return -1;
	//for every node in _hashTable[hashValue]
	for(auto x: _hashTable[hashValue]){
		//calculate Eucledian deistance
		if(metric == euclidean){ 
			distance = euclideanNorm(q,x.vec);
			count_retrieved++;
			if( distance <= min_distance){ //get minimum distance 
				min_distance = distance;
				nearest_id = x.id;
				nearest_vec = x.vec;
			}
			if( distance <= (1+approximation_factor)*app_distance){// get approximate minimum distance
				app_distance = distance;
				app_id = x.id;

			}

			if(count_retrieved > 3*L){//stop if number of items retrieved is larger than 3*L
				break;
			}
		}
		else if (metric == cosine){
			min_distance = 4;
			app_distance = 4;
			distance = 1 - cosineSimilarity(q, x.vec);//according to theory
			count_retrieved++;
			if( distance <= min_distance){
				min_distance = distance;
				nearest_id = x.id;
				nearest_vec = x.vec;
			}
			if( distance <= (1-approximation_factor)*app_distance){
				app_distance = distance;
				app_id = x.id;

			}
			if(count_retrieved > 3*L){
				break;
			}
			
		}	
	}
	auto stop_appNN = high_resolution_clock::now();//stop app NN stopwatch
	auto stop_NN = high_resolution_clock::now(); //stop NN stopwatch
	//calculate timefor approximate NN
	duration<double> duration_appNN = (stop_appNN - start_appNN);
	//calculate time for Nearest neighbor
	duration<double> duration_NN = (stop_NN - start_NN);
	//output
	cout << "Nearest neighbor: " << nearest_id << endl;
	cout << "distanceLSH: " << app_distance << endl;
	cout << "distanceTrue: " << min_distance << endl;
	cout << "tLSH: " << duration_appNN.count() << std::fixed << endl;
	cout << "tTrue: " << duration_NN.count() << std::fixed << endl;
	double ret_val= double(app_distance/min_distance); //max ratio 
	cout << "ret_val: " << ret_val << endl;
	cout << "********************************" << endl;

	
	return (double)(ret_val);
}

//displays the buckets of a hashtable
void Hash::displayHash()
{
	for(int i = 0;i<this->_tableSize;i++){
		cout << "Bucket["<<i<<"]"<<endl;
		for(auto x: _hashTable[i]){
			cout << " --> " << x.id <<",cid: "<<x.cluster_id;
		}
		cout << endl;
	}
}

int Hash::getTableSize()
{
	return this->_tableSize;
}

//Euclidean Norm
double euclideanNorm(vector_t u, vector_t v)
{
	int i;
	double sum = 0;
	for(i = 0; i < (int)v.size(); i++){
		sum = sum + pow(u[i]-v[i],2);
	}
	return sqrt(sum);
}

//Cosine Similarity
double cosineSimilarity(vector_t x, vector_t y)
{
	double dot = 0, denom_x = 0, denom_y = 0;
	dot = innerProduct(x, y);
	denom_x = innerProduct(x, x);
	denom_y = innerProduct(y, y);

	return ( dot / (sqrt(denom_x) * sqrt(denom_y)) );

}

double innerProduct(vector_t u, vector_t v)
{
	int i;
	double sum = 0;
	for(i=0; i < (int)v.size(); i++ ){
		sum = sum + u[i]*v[i];
	}

	return sum;
}



//generates g function 
double Hash::hash(vector_t p)
{
	double value = 0;
	double sum = 0;
	uniform_int_distribution<int> distribution(0,500);//initialize for R value
	for(int i=0;i<(int)this->vec_v.size();i++){ 
		value = ((innerProduct(p, this->vec_v[i]) + vec_t[i]) / double(_w));
		sum+= (distribution(generator) * value);
	}
	return sum;
}

//hash function for cosine similarity
double Hash::cosineHash(vector_t p)
{
	double value = 0;
	string hi;
	double hash_value = 0;
	for(int i = 0; i<(int)this->vec_v.size();i++){
		value = ((innerProduct(p, this->vec_v[i])));//according to theory
		if(value >= 0){ //according to theory
			hi.append("1");
		}
		else{
			hi.append("0");
		}
	}
	//get hash value from bitstring
	hash_value = stol(hi, nullptr, 2);
	
	return hash_value ;

}

double Hash::hashCUBE(vector_t p)
{
	double value = 0;
	double sum = 0;
	string hi;
	int one_zero = 0;
	int found = 0;
	uniform_int_distribution<int> distribution(0,200);//initialize for ri value
	uniform_int_distribution<int> CUBEdistribution(0,1);
	for(int i=0;i<(int)this->vec_v.size();i++){
		found = 0;
		value = ((innerProduct(p, this->vec_v[i]) + vec_t[i]) / _w);
		for(int j=0;j<this->value_vec.size();j++){
			if(value == this->value_vec[j]){ //if value already exists in value_vec
				//proceed and append the appropriate one or zero to hi string
				found = 1;
				one_zero = this->one_zero_vec[j];
				if(one_zero == 0){
					hi.append("0");
				}
				else{
					hi.append("1");
				}
				break;// exit search in value_vec
			}
		}	
		if(found == 0){// if value doesn't exist in value_vec,then add value to value_vec
			this->value_vec.push_back(value);
			one_zero = CUBEdistribution(generator); //coin tosh
			this->one_zero_vec.push_back(one_zero);//add one_zero to one_zero_vec
			if(one_zero == 0){
				hi.append("0");
			}
			else{
				hi.append("1");
			}
		}
	}
	value = stol(hi, nullptr, 2);

	return value;

}

vector_t Hash::random_vector()
{
	normal_distribution<double> distribution(0.0, 1.0); //according to theory

	vector_t vec_1;
	int i;
	double num = 0;
	for(i=0; i<DIMENSION; i++){
		num = distribution(generator);
		vec_1.push_back(num);
	}
	
	return vec_1;
}

double Hash::random_offset()
{
	uniform_int_distribution<int> distribution(0,_w); //according to theory
	double t;
	t = distribution(generator);
	return t;
}

//displays a vector
void print_vector(vector_t v)
{
	for(int i = 0;i<(int)v.size();i++){
		cout << v[i] << " ";
	}
	cout << endl;
}
