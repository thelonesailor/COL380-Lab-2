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
#include <algorithm>

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

pair<int,int> composed[2*n+2];

int maxn=n+1;
while(sp>1)
{
	while()//no significant decrease
	{
		// vertices in matching
		vector<int> M[sp+1];
		// edges in matching
		vector<pair<int,int> > Me[sp + 1];
		int next[maxn+1]={};

		for(int i=1;i<=sp;++i)
		{
			sort(Pe[i][i].begin(), Pe[i][i].end(), greater<pii>());
			int l=Pe[i][i].size();
			int taken[maxn+1];	
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
					composed[maxn]=mp(v1,v2);

					++maxn;
					taken[v1]=taken[v2]=1;
				}
			}
		}


		for(int i=1;i<=sp;++i)
		{
			for(int j=1;j<=sp;++j)
			{
				if(i!=j)
				{
					vector< pii > temp=Pe[i][j],t2;
					int s=temp.size();
					for(int k=0;k<s;++k)
					{
						int v1 = temp[k].second.first, v2 = temp[k].second.second;
						if(next[v1]==0)
						{next[v1]=v1;}
						if(next[v2]==0)
						{next[v2]=v1;}
						int a=next[v1],b=next[v2];
						// assert (a!=b);

						int f=0;
						for (int u = 0; u < t2.size(); ++u)
						{
							int a1 = t2[u].second.first, b1 = t2[u].second.second;
							if(a==a1 && b==b1)
							{
								t2[u].first+=temp[k].first;
								f=1;
								break;
							}
							else if (a == b1 && b == a1)
							{
								t2[u].first += temp[k].first;
								f=1;
								break;
							}
						}
						if(f==0)
						t2.pb(mp(temp[k].first,mp(a,b)));	
					}

					Pe[i][j]=vector <pii> (t2);
				}
				else
				{
					vector<pii> temp = Pe[i][j], t2;
					int s = temp.size();
					for (int k = 0; k < s; ++k)
					{
						int v1 = temp[k].second.first, v2 = temp[k].second.second;
						if (next[v1] == 0)
						{
							next[v1] = v1;
						}
						if (next[v2] == 0)
						{
							next[v2] = v1;
						}
						int a = next[v1], b = next[v2];
						if(a==b)
						continue;

						int f = 0;
						for (int u = 0; u < t2.size(); ++u)
						{
							int a1 = t2[u].second.first, b1 = t2[u].second.second;
							if ((a == a1 && b == b1) || (a == b1 && b == a1))
							{
								t2[u].first += temp[k].first;
								f = 1;
								break;
							}
							else if (a == b1 && b == a1)
							{
								t2[u].first += temp[k].first;
								f = 1;
								break;
							}
						}
						if (f == 0)
						t2.pb(mp(temp[k].first, mp(a,b)));

					}
					Pe[i][j] = vector<pii>(t2);
				}
			}
		}


	}

	//Do folding
	sp>>=1;	
}

	return 0;
}