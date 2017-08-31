// BAModel.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <iostream>
#include <ctime>

using namespace std;

void Initiate1(int n0, int &edge_number, int *edge1, int *edge2, vector<vector<int>> &adjacencylist) {
	vector<int>M;
	adjacencylist.push_back(M);
	for (int s = 0; s < (n0 - 1); s++) {
		adjacencylist.push_back(M);
		edge1[edge_number] = s;
		edge2[edge_number] = s + 1;
		edge_number += 1;
		adjacencylist[s].push_back(s + 1);
		adjacencylist[s + 1].push_back(s);
	}
	edge1[edge_number] = n0 - 1;
	edge2[edge_number] = 0;
	adjacencylist[n0 - 1].push_back(0);
	adjacencylist[0].push_back(n0 - 1);
	edge_number += 1;
}

void Initiate2(int n0, int &edge_number, int *edge1, int *edge2, vector<vector<int>> &adjacencylist) {
	vector<int>M;
	adjacencylist.push_back(M);
	for (int s = 0; s < n0; s++) {
		for (int t = s + 1; t < n0; t++) {
			adjacencylist.push_back(M);
			edge1[edge_number] = s;
			edge2[edge_number] = t;
			edge_number += 1;
			adjacencylist[s].push_back(t);
			adjacencylist[t].push_back(s);
		}
	}
}

void ChoseRandom(const int m, int &node_number, int &edge_number, int *nodes, int *nodes_new, int *edge1, int *edge2, vector<vector<int>> &adjacencylist, int L, int graph) {
	int i = 0;
	bool notin;
	while (i < m) {
		notin = true;
		int node = 0;
		if (graph == 1) {
			int node_arg = int(rand() / (double)(RAND_MAX + 1) * edge_number);
			int randnum = rand() % 2;
			if (randnum == 0) {
				node = edge1[node_arg];
			}
			else {
				node = edge2[node_arg];
			}
			//cout << node_arg << ' ' << randnum << ' ' << edges.size() << ' ' << node << endl;
		}
		if (graph == 2) {
			int node_arg = int(rand() / (double)(RAND_MAX + 1) * node_number);
			node = nodes[node_arg];
			}
		if (graph == 3) {
			int node_arg = int(rand() / (double)(RAND_MAX + 1) * node_number);
			node = nodes[node_arg];
			for (int l = 0; l < L; l++) {
				int node_arg = int(rand() / (double)(RAND_MAX + 1) *  adjacencylist[node].size());
				node = adjacencylist[node][node_arg];
				}
			}
		
		for (int n = 0; n < m; n++) {
			if (nodes_new[n] == node) {
				notin = false;
				}
			}
		if (notin == true) {
			nodes_new[i] = node;
			i += 1;
		}
	}
	//cout << node_number << ' ' << nodes_new[0] << ' ' << nodes_new[1] << ' ' << nodes_new[2] << ' ' << nodes_new[3] << ' ' << nodes_new[4] << endl;
}

void Iterate(int n0, int N, const int m, int &node_number, int *nodes, int *nodes_new, int &edge_number, int *edge1, int *edge2, vector<vector<int>> &adjacencylist, int L, int graph) {
	vector<int> M;
	for (int n = n0; n < N; n++) {
		node_number += 1;
		nodes[node_number] = node_number;
		adjacencylist.push_back(M);

		for (int n = 0; n < m; n++) {
			nodes_new[n] = -1;
		}

		ChoseRandom(m, node_number, edge_number, nodes, nodes_new, edge1, edge2, adjacencylist, L, graph);
		for (int t = 0; t < m; t++) {
			//cout << node_number << ' ' << nodes_new[t] << endl;
			edge1[edge_number] = node_number;
			edge2[edge_number] = int(nodes_new[t]);
			adjacencylist[node_number].push_back(int(nodes_new[t]));
			adjacencylist[int(nodes_new[t])].push_back(node_number);
			//cout << edge1[edge_number] << ' ' << edge2[edge_number] << endl;
			edge_number += 1;
			}		
	}
}

void BAModel(int initiate, int &deg_num, int n0, int N, const int m, int *nodes, int *nodes_new, int max_edge, int *edge1, int *edge2, int L, int graph, int *degrees, int *maxdeg, int &r){

	int edge_number = 0;

	vector<vector<int>> adjacencylist;

	for (int e = 0; e < max_edge; e++) {
		edge1[e] = -1;
		edge2[e] = -1;
	}

	if (initiate == 1) {
		Initiate1(n0, edge_number, edge1, edge2, adjacencylist);
	}
	else {
		Initiate2(n0, edge_number, edge1, edge2, adjacencylist);
	}

	int node_number = n0 - 1;

	for (int n = 0; n < n0; n++) {
		nodes[n] = n;
	}

	Iterate(n0, N, m, node_number, nodes, nodes_new, edge_number, edge1, edge2, adjacencylist, L, graph);

	int max = 0;

	for (int n = 0; n < N; n++) {
		degrees[deg_num] = adjacencylist[n].size();
		//cout << n << ' ' << adjacencylist[n].size() << endl;
		if (adjacencylist[n].size() > max) {
			max = adjacencylist[n].size();
		}
		deg_num += 1;
	}

	maxdeg[r] = max;

}

extern "C" __declspec(dllexport) void iterations(int N, int n0, int m, int L, int graph, int initiate, int R, int *nodes, int *nodes_new, int max_edge, int *edge1, int *edge2, int *degrees, int *maxdeg) {
//void main() {
	//const int N = 5;
	//const int n0 = 4;
	//const int m = 3;
	//const int R = 2;
	//int L = 1;
	//int graph = 1;
	//int initiate = 2;
	//int nodes[N];
	//int degrees[N];
	//int nodes_new[m];
	//int edge1[10000];
	//int edge2[10000];
	//int max_edge = 10000;
	//int maxdeg[R];

	srand(time(0));

	int deg_num = 0;

	for (int r = 0; r < R; r++) {
		BAModel(initiate, deg_num, n0, N, m, nodes, nodes_new, max_edge, edge1, edge2, L, graph, degrees, maxdeg, r);
	}

}