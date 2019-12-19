#include "incremental_test.h"
using namespace std;

bool isIndexEqual(BiGraph& g1, BiGraph& g2) {
	if (g1.left_index.size() != g2.left_index.size()) return false;
	for (int i = 0; i < g1.left_index.size(); i++) {
		if (g1.left_index[i].size() != g2.left_index[i].size()) return false;
		for (int j = 0; j < g1.left_index[i].size(); j++) {
			if (g1.left_index[i][j] != g2.left_index[i][j]) {
				return false;
			}
		}
	}
	if (g1.right_index.size() != g2.right_index.size()) return false;
	for (int i = 0; i < g1.right_index.size(); i++) {
		if (g1.right_index[i].size() != g2.right_index[i].size()) return false;
		for (int j = 0; j < g1.right_index[i].size(); j++) {
			if (g1.right_index[i][j] != g2.right_index[i][j]) {
				return false;
			}
		}
	}
	return true;
}

void random_incremental_correctness_test(BiGraph& g, int num_insert_edges) {
	cout << "start random incremental correctness test" << endl;
	g.left_index.clear();
	g.right_index.clear();
	int left_id_range = g.getV1Num() - 1;
	int right_id_range = g.getV2Num() -1;
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis_left(0, left_id_range);
	uniform_int_distribution<int> dis_right(0, right_id_range);
	static_construct_index_with_peeling_algorithm_acc(g);
	for (int i = 0; i < num_insert_edges;) {
		vid_t u = dis_left(gen);
		vid_t v = dis_right(gen);
		if (g.isEdge(u, v)) continue;
		i++;
		//dynamic_KKCore_incremental(g, u, v);
		//g.addEdge(u, v);
		//static_construct_index_with_peeling_algorithm_acc(g);
		//u = 4;
		//v = 4;
		i++;
		BiGraph copy = g;
		copy.left_index.clear();
		copy.right_index.clear();
		copy.addEdge(u, v);
		static_construct_index_with_peeling_algorithm_acc(copy);
		dynamic_KKCore_incremental_basic(g, u, v);
		if (!isIndexEqual(g, copy)) {
			cout << "error" << endl;
			g.dir = "..//error//random1_";
			g.print();
		}
	}
	cout << "finished" << endl;
}