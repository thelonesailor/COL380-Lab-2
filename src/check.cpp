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
using namespace std;

#define pii pair<int, pair<int, int>>
#define mp make_pair
#define pb push_back


int main()
{
	ifstream fin;
	ofstream result;

	fin.open("sample_test_file", ios::in);
	result.open("out1.txt", ios::out);

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
				temp += line[i];
			else
			{
				adj[l].push_back(atoi(temp.c_str()));
				temp = "";
			}
		}

		adj[l].push_back(atoi(temp.c_str()));
		temp = "";
	}

	
}