#include<iostream>
#include<stdlib.h>
 
using namespace std;
 
// A function to generate random graph.
void GenerateRandGraphs(int NOE, int NOV)
{
	int i, j, edge[NOE][2], count;
	i = 0;
	// Build a connection between two random vertex.
	while(i < NOE)
	{
		edge[i][0] = rand()%NOV+1;
		edge[i][1] = rand()%NOV+1;
 
		if(edge[i][0] == edge[i][1])
			continue;
		else
		{
			for(j = 0; j < i; j++)
			{
				if((edge[i][0] == edge[j][0] && edge[i][1] == edge[j][1]) || (edge[i][0] == edge[j][1] && edge[i][1] == edge[j][0]))
					i--;
			}
		}
		i++;
	}
 
	// Print the random graph.
	// cout<<"\nThe generated random random graph is: ";
	for(i = 0; i < NOV; i++)
	{
		count = 0;
		// cout<<"\n\t"<<i+1<<"-> { ";
		for(j = 0; j < NOE; j++)
		{
			if(edge[j][0] == i+1)
			{
				cout<<edge[j][1]<<" ";
				count++;
			}
			else if(edge[j][1] == i+1)
			{
				cout<<edge[j][0]<<" ";
				count++;
			}
			else if(j == NOE-1 && count == 0)
				cout<<"Isolated Vertex!";
		}
		// cout<<" }\n";
		cout<<endl;
	}	
}
 
int main(int argc, char const *argv[])
{

	int n, i, e, v;
	
	if(argc != 3)
	{
		printf("Usage: ./randg <n> <e>\n");
		exit(1);
	} 

	// cout<<"Random graph generation: ";
 
	// Randomly assign vertex and edges of the graph.
	// v = 11+rand()%10;
	v = atoi(argv[1]);
	// cout<<"\nThe graph has "<<v<<" vertexes.";
	// printf("\nThe graph has %d vertices.\n", v);
	// exit(1);
	// e = rand()%((v*(v-1))/2);
	// cout<<argv[2];
	e = (atoi(argv[2]));
	// %((v*(v-1))/2);
	// cout<<"The graph has "<<e<<" edges.";
	// printf("The graph has %d edges.\n", e);
 	
 	cout<<v<<" "<<e<<endl;

	// A function to generate a random undirected graph with e edges and v vertexes.
	GenerateRandGraphs(e, v);
}