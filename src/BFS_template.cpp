int sqrtp; // Number of diagonal processors
unordered_map<int, int> partition;

vector<int> g(n+1);
set<int> vp[sqrtp+1];
vector<int> V[sqrtp+1];

vector<vector<int> > adj(n+1);
//--------------------------------------------------------------------------------//
vector<int> visited(n+1);
for(int i = 1; i <= n; i++)
{
	visited[i] = 0;
	g[i] = sqrtp-1;
}

int tot_wt = n/sqrtp;

int vt = 1;
	// Get start node. Code adapted from GeeksForGeeks - https://www.geeksforgeeks.org/bfs-using-stl-competitive-coding/
	// for(auto ut = set_partition[k].begin(); ut != set_partition[k].end(); ut++)
	for(int cur_partition = 0; cur_partition < sqrtp-1; cur_partition++)
	{		
		int bfs_wt = 0;
		while(bfs_wt < tot_wt)
		{
			while(visited[vt] == 1)
				vt++;
			int u = vt;
			visited[u] = 1;
			// partition[u] = cur_partition;
			g[u] = cur_partition;
			vp[cur_partition].insert(u);
			V[cur_partition].push_back(u);
			// Total weight
			bfs_wt++;
			if(bfs_wt >= tot_wt)
				break;
			queue<int> q;
			q.push(u);
			int chk = 0;
			
			// Start BFS
			while(!(q.empty()))
			{
				int f = q.front();
		        q.pop();
		 
		 		// printf("f is %d\n", f);
		        // Enqueue all adjacent of f and mark them visited 
		        for (auto i = adj[f].begin(); i != adj[f].end(); i++) {
		        	// Ensure - vertex is not visited yet
		            if (visited[*i] == 0) {
		                q.push(*i);
		                g[*i] = cur_partition;
		                vp[cur_partition].insert(*i);
		                V[cur_partition].push_back(*i);
		                visited[*i] = 1;
		                // cout<<*i<<" ";
		                bfs_wt++;
		                if(bfs_wt >= tot_wt)
		                {
		                	chk = 1;
		                	break;
		                }

		            }
		        }
		        if(chk == 1)
		        	break;
			}
			if(chk == 1)
				break;
		}
	}

// Insert all remaining vt elements
// while(vt <= n)
for(int i = 1; i <= n; ++i)
{
	if(visited[i] == 0)
	{
		vp[cur_partition].insert(i);
		V[cur_partition].push_back(i);
	}
	// vt++;
}

//--------------------------------------------------------------------------------//