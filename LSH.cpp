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
#include "LSH.h"
#include "hash.h"



using namespace std;
using namespace std::chrono;

LSH::LSH(int k, int L, string inputFile, string outputFile, string queryFile,Metric metric)
{
	/*initialize parameters*/
	_k = k;
	_L = L;

	_inputFile = inputFile;
	cout << "Input file " << _inputFile << endl;
	_outputFile = outputFile;
	_queryFile = queryFile;
	_arrayOfHashTables = new Hash*[_L]; //create array of pointers to hashtables
	ifstream inFile(_inputFile);
	int count = 0;
	string line;
	if(metric == euclidean) { //calculate number of input vectors to define table size
		if(inFile.is_open()) {
			while(getline(inFile, line))
				++count;
		}
		cout << "Count " << count << endl;
		_hashTableSize = count / 250;

	}
	else if(metric == cosine){ //for cosine ,table size is 2^k
		_hashTableSize = pow(2,_k);
	}
	for(int i = 0;i<_L;i++){ 
		_arrayOfHashTables[i] = new Hash(_hashTableSize,_k); //create hash tables	
	}
}

void LSH::setQueryFileName(string newName)
{
	this-> _queryFile = newName;
}

LSH::~LSH()
{
	for(int i = 0;i<_L;i++){
		delete _arrayOfHashTables[i];
	}
	delete[] _arrayOfHashTables;
}

void LSH::insertAllItems(vector<item_t> &items,Metric metric)
{
	int j=0;
	int i_l;
	long int hash_value;
	long double actualHashValue=0;
	if(metric == euclidean){
		for(i_l=0;i_l<_L;i_l++){ //fill all hash tables with items
			for(j=0;j<items.size();j++){ 
				hash_value = _arrayOfHashTables[i_l]->hash(items[j].vec);
				actualHashValue = ((hash_value % M) + M) % _hashTableSize;
				_arrayOfHashTables[i_l]->insertItem(items[j],actualHashValue);
			}
		}
	}
	else{
		for(i_l=0;i_l<_L;i_l++){
			for(j=0;j<items.size();j++){ 
				hash_value = _arrayOfHashTables[i_l]->cosineHash(items[j].vec);
				actualHashValue = hash_value;
				_arrayOfHashTables[i_l]->insertItem(items[j],actualHashValue);
			}
		}
	}
}

//displays LSH structure
void LSH::displayLSH()
{
	cout << "display lsh" << endl;
	for(int i = 0;i<this->_L;i++){	
		cout << "Printing hashtable["<<i<<"]"<<endl;
		_arrayOfHashTables[i]->displayHash();

	}
}

//Range search in given Radius
vector< vector<item_t> > LSH::rangeSearch(vector_t q, int cluster_id, double R, double C=1, Metric metric = euclidean)
{
	int i_l;
	string hash_string;
	long int hash_value;
	long double actualHashValue=0;
	vector< vector<item_t> > items; 
	if(metric == euclidean){ 
		for(i_l = 0; i_l < _L; i_l++){ //for each array hash q vector and execute range search
			hash_value = _arrayOfHashTables[i_l]->hash(q);
			actualHashValue = ((hash_value % M) + M) % _hashTableSize;
			items.push_back(_arrayOfHashTables[i_l]->traverseBucket(q, cluster_id,actualHashValue, R, C , metric));
		}
	}
	else{
		for(i_l = 0; i_l < _L; i_l++){
			hash_value = _arrayOfHashTables[i_l]->cosineHash(q);
			actualHashValue = hash_value;
			items.push_back(_arrayOfHashTables[i_l]->traverseBucket(q, cluster_id, actualHashValue, R, C , metric));
		}
	}
	return items;

}

//NN and approximate NN search
double LSH::nearestNeighbor(vector_t q,Metric metric)
{
	int i_l;
	string hash_string;
	long int hash_value;
	long double actualHashValue=0;
	double max_ratio = -5;
	double ratio = 0;
	if(metric == euclidean){ 
		for(i_l = 0; i_l < _L; i_l++){
			hash_value = _arrayOfHashTables[i_l]->hash(q);
			actualHashValue = ((hash_value % M) + M) % _hashTableSize;
			ratio = (double)_arrayOfHashTables[i_l]->nearestNeighborTraverse(q, actualHashValue, this->_L,metric);
			if(ratio != -1){ 
				if(ratio > max_ratio){
					max_ratio = ratio;
				}
			}
		}
	}
	else{
		for(i_l = 0; i_l < _L; i_l++){
			hash_value = _arrayOfHashTables[i_l]->cosineHash(q);
			actualHashValue = hash_value;
			ratio = (double)_arrayOfHashTables[i_l]->nearestNeighborTraverse(q, actualHashValue, this->_L,metric);
			if(ratio != -1){ 
				if(ratio > max_ratio){
					max_ratio = ratio;
				}
			}
		}
	}
	cout << "nearestNeighbor max_ratio: " << max_ratio << endl;
	return max_ratio;
}

int LSH::sizeofLSH()
{
	int memorySum = 0;
	for(int i = 0;i<_L;i++){
		memorySum += sizeof(_arrayOfHashTables[i]);
	}
	return memorySum;
}


//initializes command line parameters
void initParameters(int* k, int* L, std::string &input_file, std::string & output_file, std::string & query_file,string &met,int argc, char** argv)
{
	int i;
	//int k_flag = 0, L_flag = 1, input_flag = 1, output_flag = 1, query_flag = 1;
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
			*L=stoi(argv[i]);
			//L_flag = 1;
		}
		else if(strcmp(argv[i],"-metric") == 0){
			i++;
			met = argv[i];
		}
	}
}


