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

#define mp make_pair
#define pb push_back

int n, m;
int num_procs; // Number of processors, assume n = 2^(2r)


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
	fin>>n>>m;
	// printf("Number of vertices = %d, Number of edges = %d", V, E);
	// Adj list
	vector<int> adj[n];

	// Read edges
	int l;
	string line;
	// Read pending newline
	getline(fin, line);
	for(l = 1; l <=n; l++)
	{
		getline(fin, line);
		for(int i = 0; line[i] != '\0'; i++)
		{
			if(line[i] != ' ')
				adj[l].push_back(atoi(&line[i]));
		}
	}

	// Verify
	for (int i = 1; i <=n; ++i)
	{		
		for(vector<int>::iterator it = adj[i].begin(); it != adj[i].end(); it++)
		{
			printf("%d ", (*it));
		}
		printf("\n");
	}


	omp_set_num_threads(num_procs);

	int sqrt_p = (int)(sqrt(num_procs));
	int nthrds;
	#pragma omp parallel default(none) shared(nthrds, num_procs, sqrt_p)
	{
		int pi, pj;

		int tid = omp_get_thread_num();

		#pragma omp single
		{
			nthrds = omp_get_num_threads();
			if(nthrds != num_procs)
			{
				printf("Number of threads not OK; quitting\n");
				exit(1);
			}
		}

		// 2D Index of current processor
		pi = tid / sqrt_p; pj = tid % sqrt_p;

	}

int p=num_procs;
int sqrtp=sqrt(p);
int sp=sqrtp;
int phase=1;

// numbered 1 to sqrt(p)
#define pii pair<int,pair<int,int>>

vector <vector< pii > > Pe[sqrtp+1];
vector <int> V[sqrtp+1];
int g[n+1];

for(int i=1;i<=n;++i)
{
	int group=rand()%sqrtp+1;
	V[group].push_back(i);
	g[i]=group;
}

for(int i=1;i<=sqrtp;++i)
{
	for(int j=1;j<=sqrtp;++j)
	{	
		for(int k=0;k<V[i].size();++k)
		{
			int v1=V[i][k];
			for(int u=0;u<adj[v1].size();++u)
			{
				int v2=adj[v1][u];
				if(g[v2]==j)
				{
					Pe[i][j].pb(mp(1,mp(v1,v2)));
				}
			}
		}

	}
}

int maxn=n+1;
while(sp>1)
{
	while()//no significant decrease
	{
		vector<int> M[sp+1];
		vector<pair<int> > Me[sp + 1];
		int next[maxn+1];

		for(int i=1;i<=sp;++i)
		{
			sort(Pe[i][i].begin(), Pe[i][i].end(), greater<pii>());
			int l=Pe[i][i].size();
			int taken	
			for(int j=0;j<l;++j)
			{
				pair<int, int> temp=Pe[i][i][j].second;
				int v1=temp.second, v2=temp.first; 
				if(taken[v1]==0 && taken[v2]==0)
				{
					M[i].pb(v1);
					M[i].pb(v2);
					Me[i].pb(mp(v1,v2));

					//critical
					next[v1]=maxn;
					next[v2]=maxn;
					++maxn;
				}
			}
		}


		for(int i=1;i<=sp;++i)
		{
			for(int j=1;j<=sp;++j)
			{
				if(i!=j)
				{

				}
				else
				{

				}
			}
		}


	}

	//Do folding
	sp>>=1;	
}

	return 0;
}