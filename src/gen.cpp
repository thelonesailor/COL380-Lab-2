#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <utility>
#include <algorithm>

using namespace std;

int main(int argc,char const *argv[])
{
	int n=atoi(argv[1]),m=n*2;

	int a[n+1][n+1];
	int k=0;
	while(k<m)
	{
		int i=rand()%n+1;
		int j = rand() % n + 1;
		if(a[i][j]==0)
		{a[i][j]=1;++k;}
	}

	cout<<n<<' '<<m<<endl;
	for(int i=1;i<=n;++i)
	{
		for(int j=1;j<=n;++j)
		{
			if(a[i][j] && i!=j)
		cout<<j<<' ';

		}
		cout<<endl;
	}
	return 0;
}