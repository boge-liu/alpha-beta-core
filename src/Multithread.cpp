#include "Multithread.h"

using namespace std;

vector<vector<atomic_int>> atomic_index;

//ofstream result("result_multithread.txt");

void multithread_compute_beta_max(int alpha_start, int alpha_end, BiGraph& g) {
	for (int alpha = alpha_start; alpha <= alpha_end; alpha++) {
		alphaCopyPeel_for_kcore(alpha, g);
	}
}

void multithread_compute_alpha_max(int beta_start, int beta_end, BiGraph& g) {
	inv(g);
	for (int alpha = beta_start; alpha <= beta_end; alpha++) {
		alphaCopyPeel_for_kcore(alpha, g);
	}
	inv(g);
}

void multithread() {

}

void initialize_bigraph(BiGraph& g) {
	int left_degree_max = 0;
	for (int i = 0; i < g.getV1Num(); i++) {
		if (left_degree_max < g.getV1Degree(i)) left_degree_max = g.getV1Degree(i);
	}
	int right_degree_max = 0;
	for (int i = 0; i < g.getV2Num(); i++) {
		if (right_degree_max < g.getV2Degree(i)) right_degree_max = g.getV2Degree(i);
	}
	// init g's max degree and index
	g.v1_max_degree = left_degree_max;
	g.v2_max_degree = right_degree_max;
	g.left_index.resize(g.getV1Num());
	g.right_index.resize(g.getV2Num());
	g.left_delete.resize(g.getV1Num());
	g.right_delete.resize(g.getV2Num());
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
}

void merge_result(BiGraph& target_g, vector<BiGraph>& bigraphsA, vector<BiGraph>& bigraphsB) {
	for (int i = 0; i < bigraphsA.size(); i++) {
		BiGraph& g_tmp = bigraphsA[i];
		for (vid_t u = 0; u < target_g.num_v1; u++) {
			for (int alpha = 1; alpha < target_g.left_index[u].size(); alpha++) {
				if (target_g.left_index[u][alpha] < g_tmp.left_index[u][alpha]) {
					target_g.left_index[u][alpha] = g_tmp.left_index[u][alpha];
				}
			}
		}
		for (vid_t v = 0; v < target_g.num_v2; v++) {
			for (int beta = 1; beta < target_g.right_index[v].size(); beta++) {
				if (target_g.right_index[v][beta] < g_tmp.right_index[v][beta]) {
					target_g.right_index[v][beta] = g_tmp.right_index[v][beta];
				}
			}
		}
	}
	for (int i = 0; i < bigraphsB.size(); i++) {
		BiGraph& g_tmp = bigraphsB[i];
		for (vid_t u = 0; u < target_g.num_v1; u++) {
			for (int alpha = 1; alpha < target_g.left_index[u].size(); alpha++) {
				if (target_g.left_index[u][alpha] < g_tmp.left_index[u][alpha]) {
					target_g.left_index[u][alpha] = g_tmp.left_index[u][alpha];
				}
			}
		}
		for (vid_t v = 0; v < target_g.num_v2; v++) {
			for (int beta = 1; beta < target_g.right_index[v].size(); beta++) {
				if (target_g.right_index[v][beta] < g_tmp.right_index[v][beta]) {
					target_g.right_index[v][beta] = g_tmp.right_index[v][beta];
				}
			}
		}
	}
}

void multithread_index_construction(BiGraph& g, int threads) {
	cout << "multithread thread_num: " << threads << endl;
	initialize_bigraph(g);
	vector<BiGraph> bigraphsA;
	vector<BiGraph> bigraphsB;
	vector<thread> threadgroupsA;
	vector<thread> threadgroupsB;
	int km = maxkcorenumber_test_optimal(g);
	int gap = ceil(km * 1.0 / (threads / 2));
	for (int i = 0; i < threads / 2; i++) {
		bigraphsA.push_back(g);
		bigraphsB.push_back(g);
	}
	auto start = chrono::system_clock::now();
	for (int i = 0; i < threads / 2; i++) {
		threadgroupsA.push_back(thread(multithread_compute_beta_max, i * gap + 1, MIN(km, i * gap + gap), ref(bigraphsA[i])));
	}
	for (int i = 0; i < threads / 2; i++) {
		threadgroupsB.push_back(thread(multithread_compute_alpha_max, i * gap + 1, MIN(km, i * gap + gap), ref(bigraphsB[i])));
	}
	for (int i = 0; i < threadgroupsA.size(); i++) {
		threadgroupsA[i].join();
	}
	for (int i = 0; i < threadgroupsB.size(); i++) {
		threadgroupsB[i].join();		
	}
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	merge_result(g, bigraphsA, bigraphsB);
	cout << "running time: " << elapsed_seconds.count() << endl;
}

void multithread_test(BiGraph& g, int threads) {
	BiGraph tmp = g;
	auto start = chrono::system_clock::now();
	int delta = coreIndexKCore(g);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	//result << "delta: " << delta << endl;
	cout << "nonmultithread running time: "  << elapsed_seconds.count() << endl;
	//result << "nonmultithread running time: " << elapsed_seconds.count() << endl;
	start = chrono::system_clock::now();
	multithread_index_construction(tmp, threads);
	end = chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "multithread running time: " << elapsed_seconds.count() << endl;
	//result << "multithread running time: " << elapsed_seconds.count() << endl;
	if (isIndexEqual(g, tmp)) {
		cout << "correct" << endl;
	}
	else {
		cout << "error" << endl;
	}
	//result.close();
}

void pure_multithread_test(BiGraph& g, int threads, ofstream& output_file) {
	auto start = chrono::system_clock::now();
	multithread_index_construction(g, threads);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "multithread total running time: " << elapsed_seconds.count() << endl;
	//result << "multithread total running time: " << elapsed_seconds.count() << endl;
}