//************************************************************//
//Algorithm for RECURSIVE bisection of coarsened graph. Serial//
//************************************************************//

#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
// #include <list>
// #include <utility>
// #include <unordered_map>
// #include <unordered_set>
// #include <algorithm>

using namespace std;

#define mp make_pair
#define pb push_back

#define A 0
#define B 1

int n, m;
int num_procs; // Number of processors, assume n = 2^(2r)
// Number of partitions - is a power of two
int num_partitions;


// Hash function
struct HASH{
  size_t operator()(const pair<int,int>&x)const{
    return hash<long long>()(((long long)x.first)^(((long long)x.second)<<32));
  }
};

bool sortbysecdesc(const pair<pair<int, int>,int> &a,
                   const pair<pair<int, int>,int> &b)
{
       return a.second>b.second;
}

// k - current partition number. Use it to index into set_partition
void bisect(int k, unordered_set<int> set_partition[], unordered_map< pair<int, int>, int, HASH >& edge_map, vector<int> adjm[], vector<int>& v_wt);

int check_pow_2(int n)
{
	if(n <= 0)
		return 1;
	if(n == 1)
		return 0;
	int x = n;
	while(x > 1)
	{
		if(x%2 == 1)
			return 1;
		x /= 2;
	}
	return 0;
}

int main(int argc, char const *argv[])
{
	if(argc != 4)
	{
		printf("Usage: ./bisect <infile> <outfile> <num_partitions>\n");
		exit(1);
	}

	num_partitions = atoi(argv[3]);
	// Check - power of two
	if(check_pow_2(num_partitions) == 1)
	{
		printf("ERROR - number of partitions must be a power of 2\n");
		exit(1);
	}

	ifstream fin; ofstream fout;

	fin.open(argv[1], ios::in);
	fout.open(argv[2], ios::out);	

	// Read number of vertices and edges
	fin>>n>>m;
	// printf("Number of vertices = %d, Number of edges = %d", V, E);
	// Adj list
	vector<int> adj[n+1];

	// NOTE - assuming, edge pair:first element <= second element
	unordered_map< pair<int, int>, int, HASH > edge_map;  // Store edge weights
	// Vertex weights also
	vector<int> v_wt(n+1);

	// Read edges
	int l;
	string line;
	// Read pending newline
	getline(fin, line);
	for(l = 1; l <=n; l++)
	{
		// Set vertex weight
		v_wt[l] = 1;
		getline(fin, line);
		for(int i = 0; line[i] != '\0'; i++)
		{
			int ival = 0;
			if(line[i] != ' ')
			{
				ival = atoi(&line[i]);
				adj[l].push_back(ival);
				int t1 = l; int t2 = ival;
				if(t1 > t2)
				{
					int tmp = t1;
					t1 = t2;
					t2 = tmp;
				}
				// Insert new value OR make existing value 1
				edge_map[make_pair(t1, t2)] = 1;
			}
			

		}
	}

	// Verify
	// for (int i = 1; i <=n; ++i)
	// {		
	// 	for(vector<int>::iterator it = adj[i].begin(); it != adj[i].end(); it++)
	// 	{
	// 		printf("%d ", (*it));
	// 	}
	// 	printf("\n");
	// }

	printf("Number of vertices = %d, number of edges = %d\n", n, m);
	// printf("Following are edges with their weights:\n");
	// for(auto it = edge_map.begin(); it != edge_map.end(); it++)
		// printf("(%d, %d) -> %d\n", it->first.first, it->first.second, it->second);
	
	// exit(0);

	/*
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
	*/

//--------------------------------------------------------------------------------//

	/*	Compute first bisection 
	 *	This will be parallelized later
	 */		
	
	int Vm, Em; 	// New V and E	

	// vector<int> adjm[Vm+1];							// Adjacency list
	// vector< pair<int, pair<int, int> > > edge_wt;	// Edge weights	
	// Create map of edges to weights
	// for(int i = 0; i < edge_wt.size(); i++)
	// {
	// 	// edge_map[edge_wt[i].second] = edge_wt[i].first;
	// 	pair< pair<int, int>, int> ins (edge_wt[i].second, edge_wt[i].first);
	// 	edge_map.insert( ins );
	// }

	// Sets to store partitions
	unordered_set<int> set_partition[2*num_partitions-1];
	
	// Store the initial set of vertices
	for(int i = 1; i <= n; ++i)
		set_partition[0].insert(i);

	// recursively call bisect
	bisect(0, set_partition, edge_map, adj, v_wt);

	// Print all partitions
	printf("Printing %d partitions\n", num_partitions);
	for(int i = num_partitions - 1; i < 2*num_partitions - 1; i++)
	{
		printf("Partition %d is size %d\n", i - num_partitions + 1, set_partition[i].size());
		// for(auto it = set_partition[i].begin(); it != set_partition[i].end(); it++)
			// printf("%d ", *it);
		// printf("\n");
	}

	/*

	// Local to function
	unordered_map<int, int> partition;				// Store partition - 0 or 1	

	// Store D values of vertices - this is now unordered_map
	unordered_map<int, int> D;
	// Total external cost
	int T = 0;
	// Store markings of vertices - marked is 1, unmarked is 0
	unordered_map<int, int > mark_v;

	// Compute D[v]. Iterate over vertices in set	
	for(unordered_set<int>::iterator it = S0.begin(); it != S0.end(); it++)		
	{
		int i = *it;
		mark_v[i] = 0;	// Initialize mark_v to zero
		int I = 0;
		int E = 0;
		int part_v = partition[i];
		// Search in j's adjacency list
		for(int j = 0; j < adjm[i].size(); j++)
		{
			int j2 = adjm[i][j];
			if(S0.find(j2) == S0.end())	// Should be in same S0 partition/induced subgraph
				continue;
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

	// Store two partitions of vertices
	unordered_set<int> set0;
	unordered_set<int> set1;

	for(unordered_set<int>::iterator it = S0.begin(); it != S0.end(); it++)
	{
		int cset = partition[*it];
		if(cset == 1)
			set0.insert(*it);
		else
			set1.insert(*it);
	}

	// int G = 0;	// Total gain	
	// gmap is a map
	unordered_map <pair<int, int>, int, HASH > gmap;
	// unordered_map <int, pair<int, int> > r_gmap; // Inverse gmap
	
	int gmax = 0;

	// Kernighan-Lin algorithm follows. Can be converted to GGP, GGGP etc
	do
	{		
		unordered_set<int>::iterator it1, it2;
		// Make a linked list
		// list<int> ord_g_list;
		// vector<int> ord_g;
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
					g_edge.push_back(make_pair(ij, gmap[ij]));
					// r_gmap[gmap[ij]] = ij;
					// ord_g.push_back(gmap[ij]);
					// if(count == 0)
					// {
					// 	max_gain = gmap[ij];
					// 	imax = i; jmax = j;
					// }
					// else
					// {
					// 	if(gmap[ij] > max_gain)
					// 	{
					// 		max_gain = gmap[ij];
					// 		imax = i; jmax = j;
					// 	}
					// }
					// count++;					
				}
				else
				{
					gmap[ij] = D[i] + D[j];
					g_edge.push_back(make_pair(ij, gmap[ij]));
				}				
			}
		}
		
		// Sort
		// sort(ord_g.begin(), ord_g.end());
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
		vector< pair< pair<int, int>, int > > set_pairs;		
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
			for(unordered_set<int>::iterator it = S0.begin(); it != S0.end(); it++)
			{
				int i = *it;
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
					int em1 = 0; int em2 = 0;
					if(edge_map.find(e1) != edge_map.end())
						em1 = edge_map[e1];
					if(edge_map.find(e2) != edge_map.end())
						em2 = edge_map[e2];
					D[i] += 2*(em1 - em2);
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
					int em1 = 0; int em2 = 0;
					if(edge_map.find(e1) != edge_map.end())
						em1 = edge_map[e1];
					if(edge_map.find(e2) != edge_map.end())
						em2 = edge_map[e2];
					D[i] += 2*(em2 - em1);
				}
			}

			// Update the values in the list, after updating D's
			for(int lt = vec_cur; lt < g_edge.size(); lt++)
			{
				int a1 = g_edge[lt].first.first; int b1 = g_edge[lt].first.second;
				if(a1 > b1)
				{
					int tmp = a1;
					a1 = b1;
					b1 = tmp;
				}
				int em = 0;
				if(edge_map.find(g_edge[lt].first) != edge_map.end())
					em = edge_map[g_edge[lt].first];
				g_edge[lt].second = D[a1] + D[b1] - 2*em;
				gmap[g_edge[lt].first] = g_edge[lt].second; // Update gmap also; maybe not needed
			}

			// Sort again, beginning from vec_cur
			sort(g_edge.begin()+vec_cur, g_edge.end(), sortbysecdesc);

			num_marked += 2;
		}
		// Now, find gmax
		int g_sum =set_pairs[0].second;
		gmax = set_pairs[0].second;
		int kmax = 0;
		int cnt = 1;
		for(vector<pair<pair<int, int>, int> >::iterator it = set_pairs.begin()+1; it != set_pairs.end(); ++it)
		{
			g_sum += it->second;
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

	*/


//--------------------------------------------------------------------------------//

	return 0;
}

// k - current partition number, S - set of ones in the current set
void bisect(int k, unordered_set<int> set_partition[], unordered_map< pair<int, int>, int, HASH >& edge_map, vector<int> adjm[], vector<int> &v_wt)
{
	if(k >= num_partitions-1 && k < 2*num_partitions-1)
		return;	// Base case
	if(num_partitions == 1)
		return;
	// Local to function
	unordered_map<int, int> partition;				// Store partition - 0 or 1

	int tot_wt = 0;
	for(auto it = set_partition[k].begin(); it != set_partition[k].end(); it++)
	{
		tot_wt += v_wt[*it];
		// Also set partition initially to zero
		partition[*it] = 0;
	}
	
	// COMPUTE PARTITION - do a BFS
	queue<int> q;
	// Get start node
	int u = *(set_partition[k].begin());
	partition[u] = 1;
	q.push(u);
	// Total weight
	int bfs_wt = 1;
	// Start BFS
	int chk = 0;
	while(bfs_wt <= tot_wt/2)
	{
		int f = q.front();
        q.pop();
 
        // cout << f << " ";
 
        // Enqueue all adjacent of f and mark them visited 
        for (auto i = adjm[f].begin(); i != adjm[f].end(); i++) {
        	// Ensure - vertex is in the set
            if (partition[*i] == 0 && set_partition[k].find(*i)!=set_partition[k].end()) {
                q.push(*i);
                partition[*i] = 1;
                bfs_wt += v_wt[*i];
                if(bfs_wt > tot_wt/2)
                {
                	chk = 1;
                	break;
                }

            }
        }
        if(chk == 1)
        	break;
	}
	// printf("\n");
	// Store D values of vertices - this is now unordered_map
	unordered_map<int, int> D;
	// Total external cost
	int T = 0;
	// Store markings of vertices - marked is 1, unmarked is 0
	unordered_map<int, int > mark_v;

	

	// T has been counted twice, so divide by two - CHECK
	// T /= 2;	

	// Store two partitions of vertices
	// unordered_set<int> set0;
	// unordered_set<int> set1;

	int k1 = 2*k+1; int k2 = 2*k+2;

	for(unordered_set<int>::iterator it = set_partition[k].begin(); it != set_partition[k].end(); it++)
	{
		int cset = partition[*it];
		if(cset == A)	// Check
			set_partition[k1].insert(*it);
		else
			set_partition[k2].insert(*it);
	}

	// int G = 0;	// Total gain	
	// gmap is a map
	unordered_map <pair<int, int>, int, HASH > gmap;
	// unordered_map <int, pair<int, int> > r_gmap; // Inverse gmap
	
	int gmax = 0;

	// Kernighan-Lin algorithm follows. Can be converted to GGP, GGGP etc
	int num_steps = 0;
	do
	{		
		// Compute D[v]. Iterate over vertices in set	
		for(auto it = set_partition[k].begin(); it != set_partition[k].end(); it++)		
		{
			int i = *it;
			mark_v[i] = 0;	// Initialize mark_v to zero
			int I = 0;
			int E = 0;
			int part_v = partition[i];
			// Search in j's adjacency list
			for(int j = 0; j < adjm[i].size(); j++)
			{
				int j2 = adjm[i][j];
				if(set_partition[k].find(j2) == set_partition[k].end())	// Should be in same set_partition[k] partition/induced subgraph
					continue;
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
			// T += E;
		}

		// Computed D[v]

		// unordered_set<int>::iterator it1, it2;
		// Make a linked list
		// list<int> ord_g_list;
		// vector<int> ord_g;
		vector<pair<pair<int, int>, int> > g_edge;	// Vector of pair-to-g values
		for(auto it1 = set_partition[k1].begin(); it1 != set_partition[k1].end(); it1++)
		{
			for(auto it2 = set_partition[k2].begin(); it2 != set_partition[k2].end(); it2++)
			{
				int i = *it1; int j = *it2;
				
				pair<int, int> ij = make_pair(i, j);
				pair<int, int> search_edge = make_pair(i, j);
				if(i > j)
					search_edge = make_pair(j, i);

				if(edge_map.find(search_edge) != edge_map.end())
				{
					gmap[ij] = D[i] + D[j] - 2*edge_map[search_edge];
					g_edge.push_back(make_pair(ij, gmap[ij]));
					// r_gmap[gmap[ij]] = ij;
					// ord_g.push_back(gmap[ij]);
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
					gmap[ij] = D[i] + D[j];
					g_edge.push_back(make_pair(ij, gmap[ij]));
				}				
			}
		}
		
		// Sort
		// sort(ord_g.begin(), ord_g.end());
		// Sort by g values
		sort(g_edge.begin(), g_edge.end(), sortbysecdesc);
		// Push all sorted elements into list (doubly linked)
		// for(vector<int>::iterator i = 0; i != ord_g.end(); ++i)
			// ord_g_list.push_back(*i);

		int num_marked = 0;
		int stop = set_partition[k1].size();
		if(set_partition[k2].size() < stop)
			stop = set_partition[k2].size();
		int G = 0;
		// See - where in the vector we are
		int vec_cur = 0;
		int g_val;
		// Set of swap pairs, with their gains
		vector< pair< pair<int, int>, int > > set_pairs;

		while(num_marked <= 2*stop)
		{
			// printf("Iteration %d\n", num_marked);
			int chk = 0;
			int a, b;		
			// Check for a proper g value that can be used, i.e. two unmarked vertices that can be removed
			while(chk == 0)
			{
				g_val = g_edge[vec_cur].second; // First unmarked value WILL be highest
				// pair<int, int> p_v = r_gmap[g_val];
				a = g_edge[vec_cur].first.first; b = g_edge[vec_cur].first.second;
				// if(mark_v[a]|mark_v[b] == 0)	// Only if both are unmarked
				if((mark_v[a] == 0 && mark_v[b]) == 0)
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

			num_marked += 2;
			if(num_marked > 2*stop)
				break;

			// Update D values	
			pair<int, int> e1, e2;
			int em1, em2;
			for(unordered_set<int>::iterator it = set_partition[k].begin(); it != set_partition[k].end(); it++)
			{
				int i = *it;
				if(mark_v[i] == 1)
					continue;	// If marked, continue				
				if(i <= a)
					e1 = make_pair(i, a);
				else
					e1 = make_pair(a, i);
				if(i <= b)
					e2 = make_pair(i, b);
				else
					e2 = make_pair(b, i);				
				em1 = 0; em2 = 0;
				if(edge_map.find(e1) != edge_map.end())
					em1 = edge_map[e1];
				if(edge_map.find(e2) != edge_map.end())
					em2 = edge_map[e2];

				if(partition[i] == A)
					D[i] += 2*(em1 - em2);
				else
					D[i] += 2*(em2 - em1);
			}

			// Update the values in the list, after updating D's
			int em;
			for(int lt = vec_cur; lt < g_edge.size(); lt++)
			{
				int a1 = g_edge[lt].first.first; int b1 = g_edge[lt].first.second;
				// If either is marked, continue
				if(mark_v[a1] == 1 || mark_v[b1] == 1)
					continue;
				pair<int, int> x = make_pair(a1, b1);
				if(a1 > b1)
				{
					int tmp = a1;
					a1 = b1;
					b1 = tmp;
				}
				em = 0;
				if(edge_map.find(make_pair(a1, b1)) != edge_map.end())
					em = edge_map[make_pair(a1, b1)];
				g_edge[lt].second = D[a1] + D[b1] - 2*em;
				gmap[x] = g_edge[lt].second; // Update gmap also; maybe not needed
			}

			// Sort again, beginning from vec_cur
			sort(g_edge.begin()+vec_cur, g_edge.end(), sortbysecdesc);			
		}

		// Now, find gmax
		int g_sum = set_pairs[0].second;
		gmax = set_pairs[0].second;
		int kmax = 0;
		int cnt = 1;
		for(vector<pair<pair<int, int>, int> >::iterator it = set_pairs.begin()+1; it != set_pairs.end(); ++it)
		{
			g_sum += it->second;
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
				// Swap ak and bk
				set_partition[k1].erase(ak); set_partition[k1].insert(bk);				
				set_partition[k2].erase(bk); set_partition[k2].insert(ak);
				partition[ak] = 1; partition[bk] = 0;
			}
		}
		num_steps++;	
	}while(gmax > 0 && num_steps <= 10);
	
	bisect(2*k+1, set_partition, edge_map, adjm, v_wt);
	bisect(2*k+2, set_partition, edge_map, adjm, v_wt);
}