#include "gephi.h"
#include <fstream>
#include <sstream>
using namespace std;

// output_location should end with
void outputKKCore(string output_location, BiGraph& g, int left_k, int right_k) {
	string edge_output_location = output_location + "edge.txt";
	string node_output_location = output_location + "node.txt";
	ofstream edge_output(edge_output_location.c_str());
	ofstream node_output(node_output_location.c_str());
	vector<bool> left_KKCore;
	vector<bool> right_KKCore;
	left_KKCore.resize(g.getV1Num());
	right_KKCore.resize(g.getV2Num());
	fill_n(left_KKCore.begin(), left_KKCore.size(), false);
	fill_n(right_KKCore.begin(), right_KKCore.size(), false);
	for (int i = 0; i < g.getV1Num(); i++) {
		if (g.left_index[i].size() > left_k) {
			if (g.left_index[i][left_k] >= right_k) {
				left_KKCore[i] = true;
			}
		}
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		if (g.right_index[i].size() > right_k) {
			if (g.right_index[i][right_k] >= left_k) {
				right_KKCore[i] = true;
			}
		}
	}
	node_output << "Id;Side" << endl;
	stringstream sstm;
	for (int i = 0; i < left_KKCore.size(); i++) {
		if (left_KKCore[i]) {
			sstm.str("");
			sstm << "left_" << i;
			node_output << sstm.str() << ";left" << endl;
		}
	}
	for (int i = 0; i < right_KKCore.size(); i++) {
		if (right_KKCore[i]) {
			sstm.str("");
			sstm << "right_" << i;
			node_output << sstm.str() << ";right" << endl;
		}
	}
	edge_output << "Source;Target;Type" << endl;
	for (int u = 0; u < g.getV1Num(); u++) {
		if (!left_KKCore[u]) continue;
		vector<vid_t> neigh = g.getV1Neighbors(u);
		int degree = g.getV1Degree(u);
		for (int j = 0; j < degree; j++) {
			vid_t v = neigh[j];
			if (right_KKCore[v]) {
				edge_output << "left_" << u << ";right_" << v << ";Undirected"<< endl;
			}
		}
	}
	edge_output.close();
	node_output.close();
}

void printGraph(BiGraph& g) {
	string output_location = "../gephi/";
	string edge_output_location = output_location + "edge.csv";
	string node_output_location = output_location + "node.csv";
	ofstream edge_output(edge_output_location.c_str());
	ofstream node_output(node_output_location.c_str());
	node_output << "Id;Side" << endl;
	stringstream sstm;
	for (int i = 0; i < g.num_v1; i++) {
		sstm.str("");
		sstm << "left_" << i;
		node_output << sstm.str() << ";left" << endl;
	}
	for (int i = 0; i < g.num_v2; i++) {
		sstm.str("");
		sstm << "right_" << i;
		node_output << sstm.str() << ";right" << endl;
	}
	edge_output << "Source;Target;Type" << endl;
	for (int u = 0; u < g.getV1Num(); u++) {
		vector<vid_t> neigh = g.getV1Neighbors(u);
		int degree = g.getV1Degree(u);
		for (int j = 0; j < degree; j++) {
			vid_t v = neigh[j];
			edge_output << "left_" << u << ";right_" << v << ";Undirected" << endl;
		}
	}
	edge_output.close();
	node_output.close();
}