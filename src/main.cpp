#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <list>
#include <utility>

using namespace std;

int V, E;

int main(int argc, char const *argv[])
{
	if(argc != 3)
	{
		printf("Usage: ./main <infile> <outfile>\n");
		exit(1);
	}

	ifstream fin; ofstream fout;

	fin.open(argv[1], ios::in);
	fout.open(argv[2], ios::out);

	// Read number of vertices and edges
	fin>>V>>E;
	// printf("Number of vertices = %d, Number of edges = %d", V, E);
	// Adj list
	vector<int> adj[V];

	// Read edges
	int l;
	string line;
	// Read pending newline
	getline(fin, line);
	for(l = 0; l < V; l++)
	{
		getline(fin, line);
		for(int i = 0; line[i] != '\0'; i++)
		{
			if(line[i] != ' ')
				adj[l].push_back(atoi(&line[i]));
		}
	}

	// Verify
	for (int i = 0; i < V; ++i)
	{		
		for(vector<int>::iterator it = adj[i].begin(); it != adj[i].end(); it++)
		{
			printf("%d ", *it);
		}
		printf("\n");
	}

	return 0;
}