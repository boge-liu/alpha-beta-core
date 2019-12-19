#include "mischievous.h"
#include "random.h"
#include "static.h"
#include "dynamic.h"
#include "gephi.h"
#include "preprocessing.h"
#include "incremental_test.h"
#include <ctime>
#include "experiment.h"
#include "paper.h"
using namespace std;
// m and q
void mischievous() {
	random_device rd;
	mt19937 gen(rd());
	int left_size = 3000; int right_size = 2000; int num_edge = 2000000;
	int round = 0;
	while (true) {
		round++;
		BiGraph g;
		g.init(left_size, right_size);
		uniform_int_distribution<int> dis_right(0, right_size - 1);
		uniform_int_distribution<int> dis_left(0, left_size - 1);
		for (int i = 0; i < num_edge;) {
			int u = dis_left(gen);
			int v = dis_right(gen);
			if (!g.isEdge(u, v)) {
				g.addEdge(u, v);
				i++;
			}
		}
		coreIndexKCore(g);
		int grid = 0;
		for (int alpha = 1; alpha <= g.v1_max_degree; alpha++) {
			int max = 0;
			for (vid_t u = 0; u < g.num_v1; u++) {				
				if (g.left_index[u].size() <= alpha) continue;
				if (g.left_index[u][alpha] > max) max = g.left_index[u][alpha];
			}
			grid += max;
		}		
		if (grid > g.num_edges)
			cout << "strange!" << endl;
		cout << grid << endl;
		//if (round % 1000000 == 0)
		//	cout << "round: " << round << endl;
	}
}

void mischievousReal() {
	for (int i = 0; i < DATASET.size(); i++) {
		string graph_location = "../realGraph/" + DATASET[i];
		BiGraph g(graph_location);
		cout << "start " << ALIAS[i] << " index_kcore construction!" << endl;
		int km = coreIndexKCore(g);
		int grid = 0;
		for (int alpha = 1; alpha <= g.v1_max_degree; alpha++) {
			int max = 0;
			for (vid_t u = 0; u < g.num_v1; u++) {
				if (g.left_index[u].size() <= alpha) continue;
				if (g.left_index[u][alpha] > max) max = g.left_index[u][alpha];
			}
			grid += max;
		}
		if (grid > g.num_edges)
			cout << "strange!" << endl;
		cout << grid << endl;
	}
}

void test_bicore_index_correctness(string graph_location, int alpha, int beta) {
	BiGraph g(graph_location);
	time_t start, end, time;
	cout << "start index construction!" << endl;

	start = clock();
	coreIndexKCore(g);
	end = clock();
	time = end - start;
	cout << "run time: " << time << endl;

	vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
	cout << "start bicore index construction!" << endl;
	start = clock();
	build_bicore_index(g, bicore_index_u, bicore_index_v);
	end = clock();
	time = end - start;
	cout << "buidling time: " << time << endl;

	vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2;
	left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2);

	start = clock();
	retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left1, right1, alpha, beta);
	end = clock();
	time = end - start;
	cout << "bicore index: " << time << endl;
	start = clock();
	findabcore_peeling(alpha, beta, g, left2, right2);
	end = clock();
	time = end - start;
	cout << "peeling: " << time << endl;
	if (!test_abcore_equal(left1, right1, left2, right2)) {
		cout << "error" << endl;
		return;
	}
}