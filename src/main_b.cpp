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
#include <set>
#include <map>

using namespace std;

#define pii pair<int,pair<int,int> >
#define mp make_pair
#define pb push_back

int n, m;
int num_procs=64; // Number of processors, assume n = 2^(2r)


int main(int argc, char const *argv[])
{
	if(argc != 4)
	{
		printf("Usage: ./main <infile> <outfile> <k>\n");
		exit(1);
	}

	ifstream fin; ofstream fout;

	fin.open(argv[1], ios::in);
	fout.open(argv[2], ios::out);
	// fin.open("sample_test_file", ios::in);
	// fin.open("in2.txt", ios::in);
	// fout.open("out1.txt", ios::out);


	// Read number of vertices and edges
	fin>>n>>m;
	// printf("Number of vertices = %d, Number of edges = %d", V, E);

	// Adj list
	cout<<n<<' '<<m<<endl;
	int a[n];
	vector<vector<int> > adj(n+1);
	
	// cout << n << ' ' << m << endl;
	int sz=n+m;

	// Read edges
	int l;
	string line;
	string temp;
	// Read pending newline
	getline(fin, line);
	for(l = 1; l <=n; l++)
	{
		getline(fin, line);
		int i=0;
		for (i = 0; line[i] == ' '; i++)
		{}

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


for(int i=1;i<=n;++i)
{
	int group=rand()%sqrtp+1;
	V[group].push_back(i);
	vp[group].insert(i);
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
				if(v2==v1)
				{continue;}
				if(g[v2]==j)
				{
					// if(v1==0 || v2==0)
					// {cout<<"error\n";}

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

int T=10;

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
		// vertices in matching
		vector<int> M[sp+1];
		// edges in matching
		vector<pair<int,int> > Me[sp + 1];

		vector<int> next(maxn+1);
		// #pragma omp parallel for shared(Pe,M,made)
		for(int i=1;i<=sp;++i)
		{
			sort(Pe[i][i].begin(), Pe[i][i].end(), greater<pii>());
			int l=Pe[i][i].size();
			vector<bool> taken(maxn+1,0);	
			for(int j=0;j<l;++j)
			{
				// pair<int, int> temp=Pe[i][i][j].second;
				int v1 = Pe[i][i][j].second.first, v2 = Pe[i][i][j].second.second;
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
					++maxn;
					//keep track of which made when
				}
					

					taken[v1]=taken[v2]=1;
				}
			}
		}

		for(int i=1;i<=sp;++i)
		{
			for(int j=1;j<=sp;++j)
			{
				edges[phase][t].insert(edges[phase][t].end(),Pe[i][j].begin(), Pe[i][j].end());
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
						{t2.pb(mp(temp[k].first,mp(a,b)));
						mpp[mp(a,b)]=t2.size()-1;}

					}

					Pe[i][j].clear();
					Pe[i][j].insert(Pe[i][j].begin(), t2.begin(), t2.end());
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
							next[v2] = v2;
						}
						int a = next[v1], b = next[v2];
						if(a>b)
						{swap(a,b);}

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

// vector <pair<int,int> > e;
// e=Pe[1][1];
set<int> s;
// l = Pe[1][1].size();
//  for (int i = 0; i < l; ++i)
// {
// 	int a = Pe[1][1][i].second.first, b = Pe[1][1][i].second.second;
// 	s.insert(a);
// 	s.insert(b);

// 	if(a==b)
// 	{
// 		// cout<<"equal error\n";
// 		// cout << a << ' ' << b << '\n';
// 	}

// }
s=vp[1];
cout<<"T= "<<T<<endl;
cout<<l<<endl;
cout << maxn << endl;
cout << s.size() << " nodes"<<endl;


// take k as input
int k = atoi(argv[3]), i = 0;
vector<set<int> > part(k+1);

//calculate any partitioning

for (auto it = s.begin(); it != s.end(); it++)
{
	part[i++].insert(*it);
	if(i==k)i=0;
	cout<<" W= "<<W[*it]<<endl;
	
}
// exit(1);

for (int i = 0; i < k; ++i)
{
	cout << part[i].size() << endl;
}
cout<<endl;
//till here



//input given: W,


//output taken: part.


vector<int> partition(maxn+1);
for(int i=0;i<k;++i)
{
for (auto it = part[i].begin(); it != part[i].end(); it++)
	{
		// std::cout << *it;
		partition[*it]=i;
	}
}

cout<<phase<<endl;
int sum=0;
for(int i=phase-1;i>=1;--i)
{
	for(int t=T;t>0;--t)
	{
		int l=made[i][t].size();
		sum+=l;
		for(int z=0;z<l;++z)
		{
			int rem = made[i][t][z],f=0;
			for(int j = 0; j < k; ++j)
			{

				std::set<int>::iterator itr;
				itr=part[j].find(rem);
				if (itr != part[j].end())
				{
					f=1;
					part[j].erase(rem);
					// pair<int,int> temp=composed[rem];
					partition[composed[rem].first] = j;
					partition[composed[rem].second] = j;
					part[j].insert(composed[rem].first);
					part[j].insert(composed[rem].second);
				}

			}
			if(f==0)
			{cout<<rem<<"could not be removed "<<i<<" "<<t<<"\n";}
		}
//input given

l2:;
int ma=-1,mi=n+2;
int mai=-1,mii=-1;

	for (int j = 0; j < k; ++j)
	{
		int temp=part[j].size();
		if(temp>ma)
		{ma=temp;mai=j;}
		if(temp<mi)
		{mi=temp;mii=j;}
	}

float a=ma-mi;
a/=n;

// printf("%f\n",a);

if(a>=0.05 && (ma-mi)>1)
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

cout<<"all made="<<sum<<endl;
int sumxx=0;
for(int i=0;i<k;++i)
{
	sumxx += part[i].size();
				 cout
			 << part[i].size() << endl;
}
cout<<sumxx<<endl;

for (int i = 0; i < k; ++i)
{
	for (auto it = part[i].begin(); it != part[i].end(); it++)
	{
		// std::cout << *it;
		// partition[*it] = i;
		if(*it>n)
		{cout<<"error "<<*it<<endl;}
	}
}

for(int i=1;i<=n;++i)
{
fout<<partition[i]<<' ';	
}

return 0;
}