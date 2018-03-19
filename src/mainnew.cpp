#include <bits/stdc++.h>
using namespace std;

#define pii pair<int,pair<int,int> >
#define mp make_pair
#define pb push_back

// Hash function
struct HASH
{
	size_t operator()(const pair<int, int> &x) const
	{
		return hash<long long>()(((long long)x.first) ^ (((long long)x.second) << 32));
	}
};

int num_partitions;
// k - current partition number, S - set of ones in the current set
void bisect(int k, vector<set<int>> &set_partition, unordered_map<pair<int, int>, int, HASH> &edge_map, vector<vector<int>> &adjm, vector<int> &v_wt)
{
	if (k >= num_partitions - 1 && k < 2 * num_partitions - 1)
		return; // Base case
	if (num_partitions == 1)
		return;
	if (set_partition[k].size() <= 1) // Partitions are too small
		return;

	// Local to function
	unordered_map<int, int> partition; // Store partition - 0 or 1

	int tot_wt = 0;
	for (auto it = set_partition[k].begin(); it != set_partition[k].end(); it++)
	{
		tot_wt += v_wt[*it];
		// Also set partition initially to zero
		partition[*it] = 0;
	}

	// COMPUTE PARTITION - do a BFS

	int bfs_wt = 0;
	// Get start node. Code adapted from GeeksForGeeks - https://www.geeksforgeeks.org/bfs-using-stl-competitive-coding/
	for (auto ut = set_partition[k].begin(); ut != set_partition[k].end(); ut++)
	{
		int u = *ut;
		if (partition[u] == 1)
			continue;
		partition[u] = 1;
		// Total weight
		bfs_wt += v_wt[u];
		if (bfs_wt >= tot_wt / 2) // Check - did we cross the bound?
			break;
		queue<int> q;
		q.push(u);
		// Start BFS
		int chk = 0;

		while (!(q.empty()))
		{
			int f = q.front();
			q.pop();

			// printf("f is %d\n", f);
			// Enqueue all adjacent of f and mark them visited
			for (auto i = adjm[f].begin(); i != adjm[f].end(); i++)
			{
				// Ensure - vertex is in the set
				if (partition[*i] == 0 && set_partition[k].find(*i) != set_partition[k].end())
				{
					q.push(*i);
					partition[*i] = 1;
					// cout<<*i<<" ";
					bfs_wt += v_wt[*i];
					if (bfs_wt >= tot_wt / 2)
					{
						chk = 1;
						break;
					}
				}
			}
			if (chk == 1)
				break;
		}
		// Check for bfs_wt again
		if (chk == 1)
			break;
	}

	// T has been counted twice, so divide by two - CHECK
	// T /= 2;

	int k1 = 2 * k + 1;
	int k2 = 2 * k + 2;

	for (auto it = set_partition[k].begin(); it != set_partition[k].end(); it++)
	{
		int cset = partition[*it];
		if (cset == 0) // Check
			set_partition[k1].insert(*it);
		else
			set_partition[k2].insert(*it);
	}

	bisect(2 * k + 1, set_partition, edge_map, adjm, v_wt);
	bisect(2 * k + 2, set_partition, edge_map, adjm, v_wt);
	return;

}
	int n, m;
	int num_procs = 16; // Number of processors, assume n = 2^(2r)


	int main(int argc, char const *argv[])
	{
		if (argc != 4)
		{
			printf("Usage: ./main <infile> <outfile> <k>\n");
			exit(1);
		}

		ifstream fin;
		ofstream fout;

		// fin.open(argv[1], ios::in);
		// fout.open(argv[2], ios::out);
		fin.open("sample_test_file", ios::in);
		// fin.open("in2.txt", ios::in);
		fout.open("out1.txt", ios::out);

		// Read number of vertices and edges
		fin >> n >> m;
		// printf("Number of vertices = %d, Number of edges = %d", V, E);

		// Adj list
		cout << n << ' ' << m << endl;
		int a[n];
		vector<vector<int>> adj(n + 1);

		// cout << n << ' ' << m << endl;
		int sz = n + m;

		// Read edges
		int l;
		string line;
		string temp;
		// Read pending newline
		getline(fin, line);
		for (l = 1; l <= n; l++)
		{
			getline(fin, line);
			int i = 0;
			for (i = 0; line[i] == ' '; i++)
			{
			}

			for (; line[i] != '\0'; i++)
			{
				if (line[i] != ' ')
					temp+=line[i];
				else
					{adj[l].push_back(atoi(temp.c_str()));temp="";}

			}

			if(temp.length()>0)
			{
			adj[l].push_back(atoi(temp.c_str()));
			temp = "";}
	}
cout<<"Graph Read\n";

	// Verify
	// for (int i = 1; i <=n; ++i)
	// {
	// 	printf("%d   ", (i));
	// 	for(vector<int>::iterator it = adj[i].begin(); it != adj[i].end(); it++)
	// 	{
	// 		printf("%d ", (*it));
	// 	}
	// 	printf("\n");
	// }


	// omp_set_num_threads(num_procs);

	// int sqrt_p = (int)(sqrt(num_procs));
	// int nthrds;
	// #pragma omp parallel default(none) shared(nthrds, num_procs, sqrt_p)
	// {
	// 	int pi, pj;

	// 	int tid = omp_get_thread_num();

	// 	#pragma omp single
	// 	{
	// 		nthrds = omp_get_num_threads();
	// 		if(nthrds != num_procs)
	// 		{
	// 			printf("%d\n",nthrds);
	// 			printf("Number of threads not OK; quitting\n");
	// 			exit(1);
	// 		}
	// 	}

	// 	// 2D Index of current processor
	// 	pi = tid / sqrt_p; pj = tid % sqrt_p;

	// }

int p=num_procs;
int sqrtp=sqrt(p);
int sp=sqrtp;

// numbered 1 to sqrt(p)
vector<vector<vector<pii> >  >Pe;
Pe.resize(sqrtp + 1);

for(int i=0;i<=sqrtp;++i)
{
	Pe[i].resize(sqrtp+1);
}

vector<int>	W(2*sz+1, 1);

// for(int i=1;i<=n;++i)
// {W[i]=1;}

vector <int> V[sqrtp+1];
set<int> vp[sqrtp+1];
vector<int> g(n+1);


// for(int i=1;i<=n;++i)
// {
// 	int group=rand()%sqrtp+1;
// 	V[group].push_back(i);
// 	vp[group].insert(i);
// 	g[i]=group;

// }
//--------------------------------------------------------------------------------//
vector<int> visited(n + 1);
for (int i = 1; i <= n; i++)
{
	visited[i] = 0;
	g[i] = sqrtp;
}

int tot_wt = n / sqrtp;

int vt = 1;
// Get start node. Code adapted from GeeksForGeeks - https://www.geeksforgeeks.org/bfs-using-stl-competitive-coding/
// for(auto ut = set_partition[k].begin(); ut != set_partition[k].end(); ut++)
for (int cur_partition = 1; cur_partition < sqrtp; ++cur_partition)
{
	int bfs_wt = 0;
	while (bfs_wt < tot_wt)
	{
		while (visited[vt] == 1)
			vt++;
		int u = vt;
		visited[u] = 1;
		// partition[u] = cur_partition;
		g[u] = cur_partition;
		vp[cur_partition].insert(u);
		V[cur_partition].push_back(u);
		// Total weight
		bfs_wt++;
		if (bfs_wt >= tot_wt)
			break;
		queue<int> q;
		q.push(u);
		int chk = 0;

		// Start BFS
		while (!(q.empty()))
		{
			int f = q.front();
			q.pop();

			// printf("f is %d\n", f);
			// Enqueue all adjacent of f and mark them visited
			for (auto i = adj[f].begin(); i != adj[f].end(); i++)
			{
				// Ensure - vertex is not visited yet
				if (visited[*i] == 0)
				{
					q.push(*i);
					g[*i] = cur_partition;
					vp[cur_partition].insert(*i);
					V[cur_partition].push_back(*i);
					visited[*i] = 1;
					// cout<<*i<<" ";
					bfs_wt++;
					if (bfs_wt >= tot_wt)
					{
						chk = 1;
						break;
					}
				}
			}
			if (chk == 1)
				break;
		}
		if (chk == 1)
			break;
	}
}

// Insert all remaining vt elements
// while(vt <= n)
for (int i = 1; i <= n; ++i)
{
	if (visited[i] == 0)
	{
		vp[sqrtp].insert(i);
		V[sqrtp].push_back(i);
	}
	// vt++;
}

// int sx=0;
// for(int i=1;i<=sqrtp;++i)
// {
// 	sx+=V[i].size();
// }
// cout<<sx<<endl;exit(1);
//--------------------------------------------------------------------------------//

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
				// if(v2==v1)
				// {continue;}
				if(g[v2]==j)
				{
					if(v1==0 || v2==0 || v1==v2)
					{cout<<"error\n";}

					if(v1<v2)
					Pe[i][j].pb(mp(1,mp(v1,v2)));
					else
					Pe[i][j].pb(mp(1, mp(v2, v1)));
				}
			}
		}

	}
}

vector<pair<int,int> > composed(2*sz+2);

int T=2000;

vector<vector< vector<int> > > made;
vector<vector<vector<pii> > >edges;

made.resize(11);
edges.resize(11);

for(int i=0;i<10;++i)
{
	made[i].resize(T+1);
	edges[i].resize(T+1);
}

// for (int i = 1; i <= sp; ++i)
// {
// 	int temp = 0;
// 	for (auto it = vp[i].begin(); it != vp[i].end(); it++)
// 	{
// 		temp += W[*it];
// 	}
// 	cout << temp << ' ';
// }
// cout << endl;

int phase=1;
int maxn=n+1;

// for (int i = 1; i < maxn+5; ++i)
// {
// 	cout << W[i] << ' ';
// }
// cout << endl;

//complexity: log(srtp)*T*p*m
while(sp>1)
{
	int t=0;
	cout<<"phase= "<<phase<<endl;
	while(t<T)//no significant decrease
	{++t;
		// cout<<t<<endl;
		// vertices in matching
		vector<int> M[sp+1];
		// edges in matching
		vector<pair<int,int> > Me[sp + 1];

		vector<int> next(maxn+1);
		#pragma omp parallel for shared(Pe,M,made)
		for(int i=1;i<=sp;++i)
		{
			// cout<<"sorting started"<<endl;
			// sort(Pe[i][i].begin(), Pe[i][i].end(), greater<pii>());


			vector<vector<pair<int,int>>> adj_list(2*sz + 1);
			// unordered_map<pair<int, int>, int, HASH> edge_map2; // Store edge weights

			for (auto it = Pe[i][i].begin(); it != Pe[i][i].end(); it++)
			{
				int a = it->second.first,b = it->second.second;
				if(a==b)
				{cout<<"error\n";continue;}
				adj_list[a].push_back(mp(b,it->first));
				adj_list[b].push_back(mp(a,it->first));

				// pair<int, int> e;
				// if (a > b)
				// swap(a,b);
				// e = make_pair(a,b);
			}


			// cout<<"sorted"<<endl;

			// int l=Pe[i][i].size();
			vector<bool> taken(maxn+1,0);	
			for(int j=1;j<maxn;++j)
			{
				if(adj_list[i].size()==0)
				continue;
				// pair<int, int> temp=Pe[i][i][j].second;
				int v1 = j,v2,mw=-1;
				for(int u=0;u<adj_list[v1].size();++u)
				{
					int w = adj_list[v1][u].second;
					if (w>mw)
					{mw=w;
						v2 = adj_list[v1][u].first;}
				}

				if(v1>v2)
				{swap(v1,v2);}
				if(taken[v1]==0 && taken[v2]==0)
				{
					M[i].pb(v1);
					M[i].pb(v2);
					Me[i].pb(mp(v1,v2));

					//critical
				#pragma omp critical
				{	next[v1]=maxn;
					next[v2]=maxn;
					W[maxn]=W[v1]+W[v2];
					// cout<<v1<<" "<<v2<<endl;
					composed[maxn]=mp(v1,v2);

					// if (vp[i].find(v1) == vp[i].end() || vp[i].find(v2) == vp[i].end())
					// {
					// 	cout << "error2\n";
					// }
					vp[i].erase(v1);
					vp[i].erase(v2);
					vp[i].insert(maxn);
					// cout<<"vp["<<i<<"] changed\n";
					made[phase][t].pb(maxn);
					// cout<<maxn<<endl;
					++maxn;
					//keep track of which made when
				}
					
					taken[v1]=taken[v2]=1;
				}
			}
		}

		// for(int i=1;i<=sp;++i)
		// {
		// 	for(int j=1;j<=sp;++j)
		// 	{
		// 		edges[phase][t].insert(edges[phase][t].end(),Pe[i][j].begin(), Pe[i][j].end());
		// 	}
		// }

// cout<<"big loop\n";
#pragma omp parallel for collapse(2)
		for(int i=1;i<=sp;++i)
		{
			for(int j=1;j<=sp;++j)
			{
				if(i!=j)
				{
					vector< pii > temp=Pe[i][j],t2;
					int s=temp.size();
					map<pair<int,int>,int> mpp;
					for(int k=0;k<s;++k)
					{
						int v1 = temp[k].second.first, v2 = temp[k].second.second;
						
						if(next[v1]==0)
						{next[v1]=v1;}
						if(next[v2]==0)
						{next[v2]=v2;}
						int a=next[v1],b=next[v2];
						if(a>b)
						{swap(a,b);}

						if (a == b)
							continue;


						pair<int,int> te;
						te.first=a;te.second=b;
						if(mpp[te]>0)
						{
							int i = mpp[te]-1;
							t2[i].first += temp[k].first;
						}
						else
						{t2.pb(mp(temp[k].first,mp(a,b)));
						mpp[mp(a,b)]=t2.size();}

						// int f=0;
						// for (int u = 0; u < t2.size(); ++u)
						// {
						// 	int a1 = t2[u].second.first, b1 = t2[u].second.second;
						// 	if(a==a1 && b==b1)
						// 	{
						// 		t2[u].first+=temp[k].first;
						// 		f=1;
						// 		break;
						// 	}
						// 	else if (a == b1 && b == a1)
						// 	{
						// 		t2[u].first += temp[k].first;
						// 		f=1;
						// 		break;
						// 	}
						// }

						// if(f==0)
						// {t2.pb(mp(temp[k].first,mp(a,b)));
						// mpp[mp(a,b)]=t2.size()-1;}

					}

					Pe[i][j].clear();
					Pe[i][j].insert(Pe[i][j].begin(), t2.begin(), t2.end());
				}
				else
				{
					vector<pii> temp = Pe[i][j], t2;
					int s = temp.size();
					map<pair<int, int>, int> mpp;
					for (int k = 0; k < s; ++k)
					{
						int v1 = temp[k].second.first, v2 = temp[k].second.second;
						// if(v1==v2)
						// {cout<<""}
						if (next[v1] == 0)
						{
							next[v1] = v1;
						}
						if (next[v2] == 0)
						{
							next[v2] = v2;
						}
						int a = next[v1], b = next[v2];
						if(a>b)
						{swap(a,b);}

						if(a==b)
						continue;

						pair<int, int> te = mp(a, b);
						if (mpp[te] > 0)
						{
							int i = mpp[te]-1;
							t2[i].first += temp[k].first;
						}
						else
						{
							// cout<<"mpp "<<mpp[te]<<endl;
							if(a==b)
							{cout<<"error"<<endl;}
							t2.pb(mp(temp[k].first, mp(a, b)));
							mpp[mp(a, b)] = t2.size();
						}

						// int f = 0;
						// for (int u = 0; u < t2.size(); ++u)
						// {
						// 	int a1 = t2[u].second.first, b1 = t2[u].second.second;
						// 	if ((a == a1 && b == b1) || (a == b1 && b == a1))
						// 	{
						// 		t2[u].first += temp[k].first;
						// 		f = 1;
						// 		break;
						// 	}
						// 	else if (a == b1 && b == a1)
						// 	{
						// 		t2[u].first += temp[k].first;
						// 		f = 1;
						// 		break;
						// 	}
						// }
						// if (f == 0)
						// t2.pb(mp(temp[k].first, mp(a,b)));

					}
					Pe[i][j].clear();
					Pe[i][j].insert(Pe[i][j].begin(), t2.begin(), t2.end());
					// Pe[i][j] = t2;
				}
			}
	
	// for(int i=1;i<=sp;++i)
	// {

	// 	int temp=0;
	// 	for (auto it = vp[i].begin(); it != vp[i].end(); it++)
	// 	{
	// 		temp+=W[*it];
	// 	}
	// 	cout<<temp<<' ';
	// }
	// cout << endl;

	// for(int i=1;i<maxn;++i)
	// {cout<<W[i]<<' ';}
	// cout<<endl;

		}

}

	//Do folding

	sp>>=1;
	++phase;

#pragma omp parallel for collapse(2)
	for(int i=1;i<=sp;++i)
	{
		for(int j=1;j<=sp;++j)
		{
			Pe[i][j].insert(Pe[i][j].end(), Pe[i][j+sp].begin(),Pe[i][j+sp].end());
			Pe[i][j].insert(Pe[i][j].end(), Pe[i + sp][j].begin(), Pe[i + sp][j].end());
			Pe[i][j].insert(Pe[i][j].end(), Pe[i + sp][j + sp].begin(), Pe[i + sp][j + sp].end());
			
			if(i==j)
			{
				vp[i].insert(vp[i + sp].begin(), vp[i + sp].end());
			}
		}
	}

}

#pragma omp parallel for collapse(2)
for (int i = 1; i <= sp; ++i)
{
	for (int j = 1; j <= sp; ++j)
	{
		if (i != j)
		{
			vector<pii> temp = Pe[i][j], t2;
			int s = temp.size();
			map<pair<int, int>, int> mpp;
			for (int k = 0; k < s; ++k)
			{
				int v1 = temp[k].second.first, v2 = temp[k].second.second;

				int a = v1, b = v2;
				if (a > b)
				{
					swap(a, b);
				}

				if (a == b)
					continue;

				pair<int, int> te;
				te.first = a;
				te.second = b;
				if (mpp[te] > 0)
				{
					int i = mpp[te] - 1;
					t2[i].first += temp[k].first;
				}
				else
				{
					t2.pb(mp(temp[k].first, mp(a, b)));
					mpp[mp(a, b)] = t2.size();
				}

				// int f=0;
				// for (int u = 0; u < t2.size(); ++u)
				// {
				// 	int a1 = t2[u].second.first, b1 = t2[u].second.second;
				// 	if(a==a1 && b==b1)
				// 	{
				// 		t2[u].first+=temp[k].first;
				// 		f=1;
				// 		break;
				// 	}
				// 	else if (a == b1 && b == a1)
				// 	{
				// 		t2[u].first += temp[k].first;
				// 		f=1;
				// 		break;
				// 	}
				// }

				// if(f==0)
				// {t2.pb(mp(temp[k].first,mp(a,b)));
				// mpp[mp(a,b)]=t2.size()-1;}
			}

			Pe[i][j].clear();
			Pe[i][j].insert(Pe[i][j].begin(), t2.begin(), t2.end());
		}
		else
		{
			vector<pii> temp = Pe[i][j], t2;
			int s = temp.size();
			map<pair<int, int>, int> mpp;
			for (int k = 0; k < s; ++k)
			{
				int v1 = temp[k].second.first, v2 = temp[k].second.second;
				
				int a = v1, b = v2;
				if (a > b)
				{
					swap(a, b);
				}

				if (a == b)
					continue;

				pair<int, int> te = mp(a, b);
				if (mpp[te] > 0)
				{
					int i = mpp[te] - 1;
					t2[i].first += temp[k].first;
				}
				else
				{
					// cout<<"mpp "<<mpp[te]<<endl;
					if (a == b)
					{
						cout << "error" << endl;
					}
					t2.pb(mp(temp[k].first, mp(a, b)));
					mpp[mp(a, b)] = t2.size();
				}

				// int f = 0;
				// for (int u = 0; u < t2.size(); ++u)
				// {
				// 	int a1 = t2[u].second.first, b1 = t2[u].second.second;
				// 	if ((a == a1 && b == b1) || (a == b1 && b == a1))
				// 	{
				// 		t2[u].first += temp[k].first;
				// 		f = 1;
				// 		break;
				// 	}
				// 	else if (a == b1 && b == a1)
				// 	{
				// 		t2[u].first += temp[k].first;
				// 		f = 1;
				// 		break;
				// 	}
				// }
				// if (f == 0)
				// t2.pb(mp(temp[k].first, mp(a,b)));
			}
			Pe[i][j].clear();
			Pe[i][j].insert(Pe[i][j].begin(), t2.begin(), t2.end());
			// Pe[i][j] = t2;
		}
	}
}


	// vector <pair<int,int> > e;
	// e=Pe[1][1];
	set<int> s;
	l = Pe[1][1].size();
	for (int i = 0; i < l; ++i)
	{
		int a = Pe[1][1][i].second.first, b = Pe[1][1][i].second.second;
		// s.insert(a);
		// s.insert(b);
		cout << a << " " << b << "  " << Pe[1][1][i].first << endl;
		if (a == b)
		{
			// cout<<"equal error\n";
			// cout << a << ' ' << b << '\n';
		}
	}
	s = vp[1];
	cout << "T= " << T << endl;
	cout << maxn << endl;
	cout << s.size() << " nodes" << endl;
	cout << Pe[1][1].size() << " edges" << endl;

	// take k as input
	int k = atoi(argv[3]), i = 0;
	num_partitions = k;
	vector<set<int>> part(k + 1);

	//calculate any partitioning

	// for (auto it = s.begin(); it != s.end(); it++)
	// {
	// 	part[i].insert(*it);
	// 	++i;
	// 	if(i==k)i=0;
	// 	// cout<<" W= "<<W[*it]<<endl;
	// }

	// for (int i = 0; i < k; ++i)
	// {
	// 	cout << part[i].size() << endl;
	// }
	// cout<<endl;
	//till here

	//input given: W,list of edges P[1][1],vp[1]

	vector<pair<int, pair<int, int>>> edges2;
	edges2.insert(edges2.begin(), Pe[1][1].begin(), Pe[1][1].end());
	vector<vector<int>> adj_list(maxn + 1);
	unordered_map<pair<int, int>, int, HASH> edge_map2; // Store edge weights

	for (auto it = edges2.begin(); it != edges2.end(); it++)
	{
		int a = it->second.first, b = it->second.second;
		adj_list[a].push_back(b);
		adj_list[b].push_back(a);

		if (a > b)
		{
			swap(a, b);
		}
		edge_map2[mp(a, b)] = it->first;
	}

	vector<set<int>> set_partition(2 * num_partitions - 1);

	// // Store the initial set of vertices
	for (auto i = vp[1].begin(); i != vp[1].end(); ++i)
		set_partition[0].insert(*i);

	// recursively call bisect
	bisect(0, set_partition, edge_map2, adj_list, W);

	//output taken: part.
	for (int i = 0; i < k; ++i)
	{
		part[i] = set_partition[i + num_partitions - 1];
	}

	for (int i = 0; i < k; ++i)
	{
		cout << part[i].size() << endl;
	}
	cout << endl;

	vector<int> partition(maxn + 1);
	for (int i = 0; i < k; ++i)
	{
		for (auto it = part[i].begin(); it != part[i].end(); it++)
		{
			// std::cout << *it;
			partition[*it] = i;
		}
	}

	cout << phase << endl;
	cout << "Uncoarsening started" << endl;
	int sum = 0;
	for (int i = phase - 1; i > 0; --i)
	{
		for (int t = T; t > 0; --t)
		{
			int l = made[i][t].size();
			sum += l;
			for (int z = 0; z < l; ++z)
			{
				int rem = made[i][t][z], f = 0;
				for (int j = 0; j < k; ++j)
				{

					std::set<int>::iterator itr;
					itr = part[j].find(rem);
					if (itr != part[j].end())
					{
						f = 1;
						part[j].erase(rem);
						// pair<int,int> temp=composed[rem];
						partition[composed[rem].first] = j;
						partition[composed[rem].second] = j;
						part[j].insert(composed[rem].first);
						part[j].insert(composed[rem].second);
					}
				}
				if (f == 0)
				{
					cout << rem << "could not be removed " << i << " " << t << "\n";
				}
			}
			//input given

		l2:;
			int ma = -1, mi = n + 2;
			int mai = -1, mii = -1;

			for (int j = 0; j < k; ++j)
			{
				int temp = part[j].size();
				if (temp > ma)
				{
					ma = temp;
					mai = j;
				}
				if (temp < mi)
				{
					mi = temp;
					mii = j;
				}
			}

			float a = ma - mi;
			a /= n;

			// printf("%f\n",a);

			if (a >= 0.05 && (ma - mi) > 1)
			{
				int temp = *(part[mai].begin());
				partition[temp] = mii;
				part[mii].insert(temp);
				part[mai].erase(temp);

				goto l2;
			}

			//output taken
		}
	}

	cout << "all made=" << sum << endl;
	int sumxx = 0;
	for (int i = 0; i < k; ++i)
	{
		sumxx += part[i].size();
		cout
			<< part[i].size() << endl;
	}
	cout << sumxx << endl;

	for (int i = 0; i < k; ++i)
	{
		for (auto it = part[i].begin(); it != part[i].end(); it++)
		{
			// std::cout << *it;
			// partition[*it] = i;
			if (*it > n)
			{
				cout << "error " << *it << endl;
			}
		}
	}

	for (int i = 1; i <= n; ++i)
	{
		fout << partition[i] << ' ';
	}

	return 0;
}