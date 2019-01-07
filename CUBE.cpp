#include<iostream>
#include<vector>
#include<string>
#include<string.h>
#include <fstream>
#include <sstream>
#include<math.h>
#include <algorithm>
#include <chrono>
#include "CUBE.h"
#include "hash.h"
#include "LSH.h"



using namespace std;
using namespace std::chrono;


CUBE::CUBE(int k, int m, int probes, string inputFile, string outputFile, string queryFile,Metric metric)
{
	_k = k;
	_MC = m;
	_probes = probes;
	_inputFile = inputFile;
	_queryFile = queryFile;
	_outputFile = outputFile;
	
	/*read number of vectors from input file*/
	ifstream inFile(_inputFile);
	int count = 0;
	string line;
	if(metric == euclidean) {
		if(inFile.is_open()) {
			while(getline(inFile, line))
				++count;
		}
		_hashTableSize = log2(count);

	}
	else if(metric == cosine){ 
		_hashTableSize = pow(2,_k);
	}



	_hashTable = new Hash(_hashTableSize,_k);

}

void CUBE::setQueryFileName(string newName)
{
	this-> _queryFile = newName;
}

int CUBE::get_k()
{
	return this->_k;
}

int CUBE::get_M()
{
	return this->_MC;
}

int CUBE::get_probes()
{
	return this->_probes;
}

string CUBE::get_inputFile()
{
	return this->_inputFile;
}

string CUBE::get_outputFile()
{
	return this->_outputFile;
}


//displays CUBE structure
void CUBE::displayCUBE()
{
	cout << "display CUBE" << endl;
	_hashTable->displayHash();

}

//Range search in given Radius
vector< vector<item_t> > CUBE::rangeSearch(vector_t q, int cluster_id,double R, double C=1, Metric metric = euclidean)
{
	int i_l;
	string hash_string;
	long int hash_value;
	long double actualHashValue=0;
	long int tmp_hash;
	int tmp_probes;
	tmp_probes = this->_probes; 
	vector< vector<item_t> > items; 
	if(metric == euclidean){ 
		hash_value = _hashTable->hashCUBE(q);
		actualHashValue = ((hash_value % M) + M) % _hashTableSize;
		items.push_back(_hashTable->traverseBucket(q,cluster_id ,actualHashValue, R, C , metric));
		
		tmp_probes--;
		tmp_hash = actualHashValue - 1;
		while(tmp_probes > 0){ //start checking left neighbor buckets first
			if(tmp_hash >= 0){
				items.push_back(_hashTable->traverseBucket(q, cluster_id ,actualHashValue, R, C , metric));
				tmp_probes--;
				tmp_hash--;
			}
			else{
				break;
			}
		}
		tmp_hash = actualHashValue + 1;
		while(tmp_probes > 0){ //start checking right neighbor buckets,if there are no more left neighbor buckets
			if(tmp_hash < this->_hashTableSize){
				items.push_back(_hashTable->traverseBucket(q,cluster_id , actualHashValue, R, C , metric));
				tmp_probes--;
				tmp_hash++;
			}
			else{
				break;
			}
		}
	}
	else{//cosine
		hash_value = _hashTable->cosineHash(q);
		actualHashValue = hash_value;
		items.push_back(_hashTable->traverseBucket(q,cluster_id , actualHashValue, R, C , metric));
		tmp_probes--;
		tmp_hash = actualHashValue - 1;
		while(tmp_probes > 0){ //start checking left neighbor buckets
			if(tmp_hash >= 0){
				items.push_back(_hashTable->traverseBucket(q, cluster_id ,actualHashValue, R, C , metric));
				tmp_probes--;
				tmp_hash--;
			}
			else{
				break;
			}
		}
		tmp_hash = actualHashValue + 1;
		while(tmp_probes > 0){ //start checking righ neighbor buckets
			if(tmp_hash < this->_hashTableSize){
				items.push_back(_hashTable->traverseBucket(q, cluster_id ,actualHashValue, R, C , metric));
				tmp_probes--;
				tmp_hash++;
			}
			else{
				break;
			}
		}
	}
	return items;
}

void CUBE::insertAllItems(vector<item_t> &items,Metric metric)
{
	int j=0;
	long int hash_value;
	long double actualHashValue=0;
	if(metric == euclidean){
		for(j=0;j<items.size();j++){
		 hash_value = _hashTable->hashCUBE(items[j].vec);
		 actualHashValue = ((hash_value % M) + M) % _hashTableSize;
		 _hashTable->insertItem(items[j], actualHashValue);
		}
	}
	else{
		for(j=0;j<items.size();j++){
		 hash_value = _hashTable->cosineHash(items[j].vec);
		 actualHashValue = ((hash_value % M) + M) % _hashTableSize;
		 _hashTable->insertItem(items[j], actualHashValue);
		}
	}
}

//NN and approximate NN search
double CUBE::nearestNeighbor(vector_t q,Metric metric)
{
	int i_l;
	string hash_string;
	long int hash_value;
	long double actualHashValue=0; 
	long int tmp_hash;
	int tmp_probes = this->_probes;
	double max_ratio = -5;
	double ratio = 0;
	if(metric == euclidean){ 
		hash_value = _hashTable->hashCUBE(q);// get hash value
		actualHashValue = ((hash_value % M) + M) % _hashTableSize;
		//search NN in bucket
		ratio = (double)_hashTable->nearestNeighborTraverse(q, actualHashValue, this->_MC/3,metric);
		if(ratio > max_ratio){
			max_ratio = ratio;
		}
		tmp_probes--;
		tmp_hash = actualHashValue - 1;
		while(tmp_probes > 0){ //start checking left neighbor buckets
			if(tmp_hash >= 0){
				ratio = _hashTable->nearestNeighborTraverse(q,tmp_hash, this->_MC/3,metric);
				if(ratio != -1){ 
					if(ratio > max_ratio){
						max_ratio = ratio;
					}
				}
				tmp_probes--;
				tmp_hash--;
			}	
			else{
				break;
			}
		}
		tmp_hash = actualHashValue + 1;
		while(tmp_probes > 0){ //start checking right neighbor buckets
			if(tmp_hash < this->_hashTableSize){
				ratio = _hashTable->nearestNeighborTraverse(q,tmp_hash, this->_MC/3,metric);
				if(ratio != -1){ 
					if(ratio > max_ratio){
						max_ratio = ratio;
					}
				}
				tmp_probes--;
				tmp_hash++;
			}
			else{
				break;
			}
		}

	
	}
	else{//cosine
		hash_value = _hashTable->cosineHash(q);
		actualHashValue = hash_value;
		ratio = _hashTable->nearestNeighborTraverse(q, actualHashValue, this->_MC/3,metric);
		if(ratio > max_ratio){
			max_ratio = ratio;
		}
		tmp_probes--;
		tmp_hash = actualHashValue - 1;
		while(tmp_probes > 0){ //start checking left neighbor buckets
			if(tmp_hash >= 0){
				ratio = _hashTable->nearestNeighborTraverse(q,tmp_hash, this->_MC/3,metric);
				if(ratio != -1){ 
					if(ratio > max_ratio){
						max_ratio = ratio;
					}
				}
				tmp_probes--;
				tmp_hash--;
			}
			else{
				break;
			}
		}
		tmp_hash = actualHashValue + 1;
		while(tmp_probes > 0){ //start checking righ neighbor buckets
			if(tmp_hash < this->_hashTableSize){
				ratio = _hashTable->nearestNeighborTraverse(q,tmp_hash, this->_MC/3,metric);
				if(ratio != -1){ 
					if(ratio > max_ratio){
						max_ratio = ratio;
					}
				}
				tmp_probes--;
				tmp_hash++;
			}
			else{
				break;
			}
		}
	}
	//output max ratio
	cout << "nearestNeighbor max_ratio: " << max_ratio << endl;
	return max_ratio;
}



//initialize command line parameters
void initParametersCube(int* k, int* m, int* probes, std::string &input_file, std::string & output_file, std::string & query_file,string &met,int argc, char** argv)
{
	int i;
	//int k_flag = 0, m_flag = 1, input_flag = 1, output_flag = 1, query_flag = 1;
	for(i = 1; i<argc; i++){
		if(strcmp(argv[i], "-o") == 0){
			i++;
			output_file = argv[i];
			//output_flag = 1;
		}
		else if(strcmp(argv[i],"-q") == 0){
			i++;
			query_file = argv[i];
			//query_flag = 1;
		}
		else if(strcmp(argv[i],"-d") == 0){
			i++;
			input_file = argv[i];
			//input_flag = 1;
		}
		else if(strcmp(argv[i],"-k") == 0){
			i++;
			*k=stoi(argv[i]);
			//k_flag = 1;
		}
		else if(strcmp(argv[i],"-L") == 0){
			i++;
			*m=stoi(argv[i]);
			//m_flag = 1;
		}
		else if(strcmp(argv[i],"-metric") == 0){
			i++;
			met = argv[i];
		}
		else if(strcmp(argv[i],"-probes") == 0){
			i++;
			*probes=stoi(argv[i]);
		}
	}
		
}