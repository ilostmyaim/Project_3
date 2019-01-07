OUT = cluster
OBJS = cluster_.o LSH.o CUBE.o hash.o Cluster.o KMeans.o 
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG) -std=c++11
LFLAGS = -Wall $(DEBUG)
HEADERS=$(wildcard *.h)
 
%.o: %.c $(HEADERS)
	$(CC) $(LFLAGS) -c $< -o $@

all:cluster

cluster: Cluster.o KMeans.o hash.o cluster_.o LSH.o CUBE.o
	$(CC) $(LFLAGS) Cluster.o KMeans.o hash.o cluster_.o LSH.o CUBE.o -o $@

clean: 
	rm -rf *.o $(OUT)

