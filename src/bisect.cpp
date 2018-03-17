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
#include <unordered_map>
#include <algorithm>

using namespace std;

#define mp make_pair
#define pb push_back

#define A 0
#define B 1

int n, m;
int num_procs; // Number of processors, assume n = 2^(2r)

bool sortbysecdesc(const pair<pair<int, int>,int> &a,
                   const pair<pair<int, int>,int> &b)
{
       return a.second>b.second;
}

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

//--------------------------------------------------------------------------------//

	/*	Compute first bisection 
	 *	This will be parallelized later
	 */
	
	int Vm, Em; 	// New V and E	

	vector<int> adjm[Vm+1];							// Adjacency list
	vector< pair<int, pair<int, int> > > edge_wt;	// Edge weights

	// NOTE - assuming, edge pair:first element <= second element

	unordered_map< pair<int, int>, int > edge_map;  // Store edge weights
	// unordered_map<int, int> partition;				// Store partition - 0 or 1
	vector<int> partition[Vm+1];			// No need of a map

	// Create map of edges to weights
	for(int i = 0; i < edge_wt.size(); i++)
	{
		unordered_map[edge_wt[i].second] = edge_wt[i].first;
	}	

	// Store D values of vertices
	int D[Vm+1];
	// Total external cost
	int T = 0;
	// Compute D[v]
	for(int i = 1; i <= Vm; i++)
	{
		int I = 0;
		int E = 0;
		int part_v = partition[i];
		// Search in j's adjacency list
		for(int j = 0; j < adjm[i].size(); j++)
		{
			int j2 = adm[i][j];
			int part_w = partition[j2];
			pair<int, int> e = make_pair(i, j2);
			if(i > j2)
				e = make_pair(j2, i);
			int cost = edge_map[e];
			if(part_v == part_w)
				I += cost;
			else
				E += cost;					
		}
		D[i] = E - I;
		T += E;
	}

	// T has been counted twice, so divide by two - CHECK
	T /= 2;

	// Store markings of edges - marked is 1, unmarked is 0
	vector<int> mark_v(Vm+1, 0);

	// Store two partitions of vertices
	unordered_set<int> set0;
	unordered_set<int> set1;

	for(int i = 1; i <= Vm; i++)
	{
		int cset = partition[i];
		if(cset == 1)
			set0.insert(i);
		else
			set1.insert(i);
	}

	// int G = 0;	// Total gain	
	// gmap is a map
	unordered_map <pair<int, int>, int > gmap;
	unordered_map <int, pair<int, int> > r_gmap; // Inverse gmap
	
	int gmax = 0;

	// Kernighan-Lin algorithm follows. Can be converted to GGP, GGGP etc
	do
	{
		// int max_gain = 0; int imax = 0; int jmax = 0;
		// int count = 0;
		// Make gains for all pairs of vertices BETWEEN set0 and set1
		unordered_set<int>::iterator it1, it2;
		// Make a linked list
		list<int> ord_g_list;
		vector<int> ord_g;
		vector<pair<pair<int, int>, int> > g_edge;	// Vector of pair-to-g values
		for(it1 = set0.begin(); it1 != set0.end(); it1++)
		{
			for(it2 = set1.begin(); it2 != set1.end(); it2++)
			{
				int i = *it1; int j = *it2;
				if(i > j)
					i = *it2;
					j = *it1;
				pair<int, int> ij = make_pair(i, j);
				if(edge_map.find(ij) != edge_map.end())
				{
					gmap[ij] = D[i] + D[j] - 2*edge_map[ij];
					r_gmap[gmap[ij]] = ij;
					ord_g.push_back(gmap[ij]);
					g_edge.push_back(make_pair(ij, gmap[ij]));
					/*if(count == 0)
					{
						max_gain = gmap[ij];
						imax = i; jmax = j;
					}
					else
					{
						if(gmap[ij] > max_gain)
						{
							max_gain = gmap[ij];
							imax = i; jmax = j;
						}
					}
					count++;*/					
				}
				else
				{

				}				
			}
		}
		
		// Sort
		sort(ord_g.begin(), ord_g.end());
		// Sort by g values
		sort(g_edge.begin(), g_edge.end(), sortbysecdesc);
		// Push all sorted elements into list (doubly linked)
		// for(vector<int>::iterator i = 0; i != ord_g.end(); ++i)
			// ord_g_list.push_back(*i);

		int num_marked = 0;
		int stop = set0.size();
		if(set1.size() < stop)
			stop = set1.size();
		int G = 0;
		// See - where in the vector we are
		int vec_cur = 0;
		int g_val;
		// Set of swap pairs, with their gains
		vector<pair<int, int>, int > set_pairs;		
		while(num_marked <= 2*stop)
		{
			int chk = 0;
			int a, b;
			
			// Check for a proper g value that can be used, i.e. two unmarked vertices that can be removed
			while(chk == 0)
			{
				g_val = g_edge[vec_cur].second; // First unmarked value WILL be highest
				// pair<int, int> p_v = r_gmap[g_val];
				a = g_edge[vec_cur].first.first; b = g_edge[vec_cur].first.second;
				// if(mark_v[a]|mark_v[b] == 0)	// Only if both are unmarked
				if(!(mark_v[a]|mark_v[b]))
				{
					chk = 1;
					G += g_val;
				}
				vec_cur++;
			}

			// Assume - i is in partition 0, j in 1
			if(partition[a] == B)	// Means partition[b] == A
			{
				int tmp = a;
				a = b;
				b = a;
			}

			// Mark a and b
			mark_v[a] = 1; mark_v[b] = 1;
			// Add to set_pairs
			set_pairs.push_back(make_pair(make_pair(a, b), g_val));
			// Update D values
			for(int i = 1; i <= Vm; i++)
			{
				if(mark_v[i] == 1)
					continue;	// If marked, continue
				if(partition[i] == A)
				{
					pair<int, int> e1, e2;
					if(i <= a)
						e1 = make_pair(i, a);
					else
						e1 = make_pair(a, i);
					if(i <= b)
						e2 = make_pair(i, b);
					else
						e2 = make_pair(b, i);
					D[i] += 2*(edge_map[e1] - edge_map[e2]);
				}
				else
				{
					pair<int, int> e1, e2;
					if(i <= a)
						e1 = make_pair(i, a);
					else
						e1 = make_pair(a, i);
					if(i <= b)
						e2 = make_pair(i, b);
					else
						e2 = make_pair(b, i);
					D[i] += 2*(edge_map[e2] - edge_map[e1]);	
				}
			}

			// Update the values in the list, after updating D's
			for(int lt = cur_vec; lt < g_edge.size(); lt++)
			{
				int a1 = g_edge[lt].first.first; int b1 = g_edge[lt].first.second;
				if(a1 > b1)
				{
					int tmp = a1;
					a1 = b1;
					b1 = tmp;
				}
				g_edge[lt].second = D[a1] +D[b1] - 2*edge_map[g_edge[lt].first];
				gmap[g_edge[lt].first] = g_edge[lt].second; // Update gmap also; maybe not needed
			}

			// Sort again, beginning from cur_vec
			sort(g_edge.begin()+cur_vec, g_edge.end(), sortbysecdesc);

			num_marked += 2;
		}
		// Now, find gmax
		int g_sum =set_pairs[0].second;
		gmax = set_pairs[0].second;
		int kmax = 0;
		int cnt = 1;
		for(vector<pair<pair<int, int>, int> >::iterator it = set_pairs.begin()+1; it != set_pairs.end(); ++it)
		{
			g_sum += *it.second;
			if(g_sum > gmax)
			{
				gmax = g_sum;
				kmax = cnt;
			}
			cnt++;
		}

		// Swap all values till done
		if(gmax > 0)
		{
			for(int i = 0; i < kmax; ++i)
			{
				int ak = set_pairs[i].first.first; int bk = set_pairs[i].first.second;
				set0.erase(ak);	set0.insert(bk);
				set1.erase(bk); set1.insert(ak);
			}
		}

	}while(gmax > 0);

//--------------------------------------------------------------------------------//

	return 0;
}