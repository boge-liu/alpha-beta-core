#include "preprocessing.h"
using namespace std;

void preprocessing_raw_graph(string input_file_location, string output_file_location) {
	ifstream rawgraph;
	rawgraph.open(input_file_location.c_str());
	map<string, int> left_oldid2newidmap;
	map<string, int> right_oldid2newidmap;
	string s, t;
	int left = 0;
	int right = 0;
	BiGraph g;
	while (rawgraph >> s >> t) {
		if (left_oldid2newidmap.find(s) == left_oldid2newidmap.end()) {
			left_oldid2newidmap[s] = left;
			g.neighborHash_v1.resize(g.neighborHash_v1.size() + 1);
			left++;
		}
		if (right_oldid2newidmap.find(t) == right_oldid2newidmap.end()) {
			right_oldid2newidmap[t] = right;
			g.neighborHash_v2.resize(g.neighborHash_v2.size() + 1);
			right++;
		}
		g.neighborHash_v1[left_oldid2newidmap[s]].insert(right_oldid2newidmap[t]);
		g.neighborHash_v2[right_oldid2newidmap[t]].insert(left_oldid2newidmap[s]);
	}
	rawgraph.close();

	g.dir = output_file_location;
	g.print(true);
}

void preprocessing_raw_graph_live_journal(string input_file_location, string output_file_location) {
	ifstream rawgraph;
	rawgraph.open(input_file_location.c_str());
	int s, t;
	int left = 0;
	int right = 0;
	BiGraph g;
	g.init(11514053, 11514053);
	int UMAXID = 0;
	int VMAXID = 0;
	while (rawgraph >> s >> t) {
		if (UMAXID < s) UMAXID = s;
		if (VMAXID < t) VMAXID = t;
		g.addEdge(s-1, t-1);
	}
	rawgraph.close();

	g.num_v1 = UMAXID;
	g.num_v2 = VMAXID;
	g.dir = output_file_location;
	g.print();
}

void preprocessing_raw_graph_amazon_ratings(string input_file_location, string output_file_location) {
	ifstream rawgraph;
	rawgraph.open(input_file_location.c_str());
	map<string, int> left_oldid2newidmap;
	map<string, int> right_oldid2newidmap;
	string s, t, m, n;
	int left = 0;
	int right = 0;
	while (rawgraph >> s >> t >> m >> n) {
		if (left_oldid2newidmap.find(s) == left_oldid2newidmap.end()) {
			left_oldid2newidmap[s] = left;
			left++;
		}
		if (right_oldid2newidmap.find(t) == right_oldid2newidmap.end()) {
			right_oldid2newidmap[t] = right;
			right++;
		}
	}
	rawgraph.close();

	rawgraph.open(input_file_location.c_str());
	BiGraph g;
	g.init(left_oldid2newidmap.size(), right_oldid2newidmap.size());
	while (rawgraph >> s >> t >> m >> n) {
		g.addEdge(left_oldid2newidmap[s], right_oldid2newidmap[t]);
		g.num_edges++;
	}
	g.dir = output_file_location;
	g.print();
}

void preprocessing_raw_graph_with_duplicate_edges(string input_file_location, string output_file_location) {
	ifstream rawgraph;
	rawgraph.open(input_file_location.c_str());
	int s, t, p, q; double b;
	int left = 0;
	int right = 0;
	BiGraph g;
	g.init(35444383, 35444383);
	int UMAXID = 0;
	int VMAXID = 0;
	int pres = -1; int pret = -1;
	while (rawgraph >> s >> t >> p >> q) {
		if (pres == s && pret == t) continue;
		pres = s;
		pret = t;
		if (UMAXID < s) UMAXID = s;
		if (VMAXID < t) VMAXID = t;
		g.addEdge(s - 1, t - 1);
	}
	rawgraph.close();

	g.num_v1 = UMAXID;
	g.num_v2 = VMAXID;
	g.dir = output_file_location;
	g.print();
}

void preprocessing_raw_graph_with_duplicate_edges_real(string input_file_location, string output_file_location) {
	ifstream rawgraph;
	rawgraph.open(input_file_location.c_str());
	int s, t, p, q; double b;
	int left = 0;
	int right = 0;
	BiGraph g;
	g.init(35444383, 35444383);
	int UMAXID = 0;
	int VMAXID = 0;
	int pres = -1; int pret = -1;
	while (rawgraph >> s >> t) {
		if (pres == s && pret == t) continue;
		pres = s;
		pret = t;
		if (g.isEdge(s, t)) continue;
		if (UMAXID < s) UMAXID = s;
		if (VMAXID < t) VMAXID = t;
		g.addEdge(s, t);
	}
	rawgraph.close();

	g.num_v1 = UMAXID;
	g.num_v2 = VMAXID;
	g.dir = output_file_location;
	g.print();
}