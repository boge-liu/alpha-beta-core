#include "experiment.h"
#include "paper.h"
#include "incremental_test.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include "random.h"

using namespace std;

string OUTPUT_LOC = "../experiments/";
string INPUT_LOC = "../realGraph/";
vector<Datasets> scala_dataset = { Datasets::TR, Datasets::WT, Datasets::OG };
int DYNAMIC_TEST_NUM = 10;
bool DYNAMIC_WITH_LIMIT = false;

void timetest_basic() {
	string result_location = OUTPUT_LOC + "indexConsTime_Basic.txt";
	ofstream result(result_location);
	result << "dataset" << "," << "alias" << "," << "time" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		string graph_location = INPUT_LOC + DATASET[i];
		BiGraph g(graph_location);
		time_t start, end, time;
		cout << "start " << ALIAS[i] << " index_basic construction!" << endl;
		start = clock();
		coreIndexBasic(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		result << DATASET[i] << "," << ALIAS[i] << "," << time << endl;
	}
	result.close();
}

void timetest_basic_single(string graph_name, string graph_alias) {
	string result_location = OUTPUT_LOC + "indexConsTime_Basic_single_" + graph_alias + ".txt";
	ofstream result(result_location);
	result << "dataset" << "," << "alias" << "," << "time" << endl;
	string graph_location = INPUT_LOC + graph_name;
	BiGraph g(graph_location);
	time_t start, end, time;
	cout << "start " << graph_alias << " index_basic construction!" << endl;
	start = clock();
	coreIndexBasic(g);
	end = clock();
	time = end - start;
	cout << "run time: " << time << endl;
	result << graph_name << "," << graph_alias << "," << time << endl;
	result.close();
}

void timetest_swap() {
	string result_location = OUTPUT_LOC + "indexConsTime_Swap.txt";
	ofstream result(result_location);
	result << "dataset" << "," << "alias" << "," << "time" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		string graph_location = INPUT_LOC + DATASET[i];
		BiGraph g(graph_location);
		time_t start, end, time;
		cout << "start " << ALIAS[i] << " index_swap construction!" << endl;
		start = clock();
		coreIndexSwap(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		result << DATASET[i] << "," << ALIAS[i] << "," << time << endl;
	}
	result.close();
}

void timetest_swap_single(string graph_name, string graph_alias) {
	string result_location = OUTPUT_LOC + "indexConsTime_Swap_single_" + graph_alias + ".txt";
	ofstream result(result_location);
	result << "dataset" << "," << "alias" << "," << "time" << endl;
	string graph_location = INPUT_LOC + graph_name;
	BiGraph g(graph_location);
	time_t start, end, time;
	cout << "start " << graph_alias << " index_swap construction!" << endl;
	start = clock();
	coreIndexSwap(g);
	end = clock();
	time = end - start;
	cout << "run time: " << time << endl;
	result << graph_name << "," << graph_alias << "," << time << endl;
	result.close();
}

void timetest_kcore() {
	string result_location = OUTPUT_LOC + "indexConsTime_KCore.txt";
	ofstream result(result_location);
	result << "dataset" << "," << "alias" << "," << "time" << "," << "km" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		string graph_location = INPUT_LOC + DATASET[i];
		BiGraph g(graph_location);
		time_t start, end, time;
		cout << "start " << ALIAS[i] << " index_kcore construction!" << endl;
		start = clock();
		int km = coreIndexKCore_with_pruning_rules(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		result << DATASET[i] << "," << ALIAS[i] << "," << time << "," << km << endl;
	}
	result.close();
}

void timetest_kcore_single(string graph_name, string graph_alias) {
	string result_location = OUTPUT_LOC + "indexConsTime_KCore_single_" + graph_alias + ".txt";
	ofstream result(result_location);
	result << "dataset" << "," << "alias" << "," << "time" << "," << "km" << endl;
	string graph_location = INPUT_LOC + graph_name;
	BiGraph g(graph_location);
	time_t start, end, time;
	cout << "start " << graph_alias << " index_kcore construction!" << endl;
	start = clock();
	int km = coreIndexKCore(g);
	end = clock();
	time = end - start;
	cout << "run time: " << time << endl;
	result << graph_name << "," << graph_alias << "," << time << "," << km << endl;
	result.close();
}

// i dataset index
void core2size3d(int i) {
	int alphacutoffbegin = 2;
	int betacutoffbegin = 2;
	int alphacutoffend = 147;
	int betacutoffend = 147;
	string result_location = OUTPUT_LOC + "core2size3d_" + DATASET[i] + ".txt";
	ofstream result(result_location);
	string graph_location = INPUT_LOC + DATASET[i];
	BiGraph g(graph_location);
	time_t start, end, time;
	cout << "start " << ALIAS[i] << " index construction!" << endl;
	start = clock();
	coreIndexKCore(g);
	end = clock();
	time = end - start;
	cout << "run time: " << time << endl;
	result << "alpha" << "," << "beta" << "," << "U" << "," << "V" << "," << "U+V" << endl;
	bool stop_search = true;
	for (int a = alphacutoffbegin; a <= g.v1_max_degree; a++) {
		if (a > alphacutoffend) break;
		stop_search = false;
		for (int b = betacutoffbegin; b <= g.v2_max_degree; b++) {
			if (b > betacutoffend) break;
			if (stop_search) {
				result << a << "," << b << "," << 0 << "," << 0 << "," << 0 << endl;
			}
			int left_count = 0;
			int right_count = 0;
			for (vid_t u = 0; u < g.num_v1; u++) {
				if (g.degree_v1[u] < a) continue;
				if (g.left_index[u][a] >= b) left_count++;
			}
			for (vid_t v = 0; v < g.num_v2; v++) {
				if (g.degree_v2[v] < b) continue;
				if (g.right_index[v][b] >= a) right_count++;
			}
			result << a << "," << b << "," << left_count << "," << right_count << "," << left_count + right_count << endl;
			if (left_count == 0) stop_search = true;
		}
	}
	result.close();
	cout << DATASET[i] << "," << ALIAS[i] << "," << time << endl;
}

void scalability_sample() {
	for (int i = 0; i < scala_dataset.size(); i++) {
		string output_dir = "../scalability/";
		string graph_location = INPUT_LOC + dataset2name[scala_dataset[i]];
		BiGraph g(graph_location);
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<int> dis(0, 4);
		// node sample
		for (int p = 0; p < 5; p++) {
			// generate random nodes
			vector<int> left_node_sample_map; left_node_sample_map.resize(g.num_v1); fill_n(left_node_sample_map.begin(), left_node_sample_map.size(), -1);
			vector<int> right_node_sample_map; right_node_sample_map.resize(g.num_v2); fill_n(right_node_sample_map.begin(), right_node_sample_map.size(), -1);
			int left_sample_num = 0; int right_sample_num = 0;
			for (vid_t u = 0; u < g.num_v1; u++) {
				int r = dis(gen);
				if (r <= p) {
					left_node_sample_map[u] = left_sample_num;
					left_sample_num++;
				}
			}
			for (vid_t v = 0; v < g.num_v2; v++) {
				int r = dis(gen);
				if (r <= p) {
					right_node_sample_map[v] = right_sample_num;
					right_sample_num++;
				}
			}
			BiGraph sg;
			sg.init(left_sample_num, right_sample_num);
			for (vid_t u = 0; u < g.num_v1; u++) {
				if (left_node_sample_map[u] != -1) {
					vector<vid_t> neigh = g.neighbor_v1[u];
					for (int t = 0; t < g.degree_v1[u]; t++) {
						vid_t v = neigh[t];
						if (right_node_sample_map[v] != -1) {
							sg.addEdge(left_node_sample_map[u], right_node_sample_map[v]);
						}
					}
				}
			}
			stringstream output_location; output_location << output_dir << dataset2name[scala_dataset[i]] << "node" << (p * 20 + 20);
			sg.dir = output_location.str();
			sg.print();
		}
		// edge sample
		for (int p = 0; p < 5; p++) {
			BiGraph sg;
			sg.init(g.num_v1, g.num_v2);
			for (vid_t u = 0; u < g.num_v1; u++) {
				vector<vid_t> neigh = g.neighbor_v1[u];
				for (int t = 0; t < g.degree_v1[u]; t++) {
					vid_t v = neigh[t];
					int r = dis(gen);
					if (r <= p) {
						sg.addEdge(u, v);
					}
				}
			}
			stringstream output_location; output_location << output_dir << dataset2name[scala_dataset[i]] << "edge" << (p * 20 + 20);
			sg.dir = output_location.str();
			sg.print();
		}
	}	
}

void scalabilityTest_basic(int ds) {
	//string noderesult_location = "../experiments/" + DATASET[ds] + "basic_node_scale.txt";
	string edgeresult_location = "../experiments/" + DATASET[ds] + "basic_edge_scale.txt";
	//ofstream noderesult(noderesult_location);
	ofstream edgeresult(edgeresult_location);
	// node test
	/*noderesult << "sampling_rate" << "," << "time" << endl;
	for (int p = 0; p < 4; p++) {
		stringstream graph_location; graph_location << "../scalability/" << DATASET[ds] << "node" << (p * 20 + 20);
		BiGraph g(graph_location.str());
		cout << "start " << ALIAS[ds] << " node scale " << (p * 20 + 20) << " index_basic construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexBasic(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		noderesult << double(p * 20 + 20) / 100 << "," << time << endl;
	}
	noderesult.close();*/
	// edge test
	edgeresult << "sampling_rate" << "," << "time" << endl;
	for (int p = 0; p < 4; p++) {
		stringstream graph_location; graph_location << "../scalability/" << DATASET[ds] << "edge" << (p * 20 + 20);
		BiGraph g(graph_location.str());
		cout << "start " << ALIAS[ds] << " edge scale " << (p * 20 + 20) << " index_basic construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexBasic(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		edgeresult << double(p * 20 + 20) / 100 << "," << time << endl;
	}
	edgeresult.close();
}

void scalabilityTest_swap(int ds) {
	string noderesult_location = "../experiments/" + DATASET[ds] + "swap_node_scale.txt";
	string edgeresult_location = "../experiments/" + DATASET[ds] + "swap_edge_scale.txt";
	ofstream noderesult(noderesult_location);
	ofstream edgeresult(edgeresult_location);
	// node test
	noderesult << "sampling_rate" << "," << "time" << endl;
	for (int p = 0; p < 4; p++) {
		stringstream graph_location; graph_location << "../scalability/" << DATASET[ds] << "node" << (p * 20 + 20);
		BiGraph g(graph_location.str());
		cout << "start " << ALIAS[ds] << " node scale " << (p * 20 + 20) << " index_swap construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexSwap(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		noderesult << double(p * 20 + 20) / 100 << "," << time << endl;
	}
	noderesult.close();
	// edge test
	edgeresult << "sampling_rate" << "," << "time" << endl;
	for (int p = 0; p < 4; p++) {
		stringstream graph_location; graph_location << "../scalability/" << DATASET[ds] << "edge" << (p * 20 + 20);
		BiGraph g(graph_location.str());
		cout << "start " << ALIAS[ds] << " edge scale " << (p * 20 + 20) << " index_swap construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexSwap(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		edgeresult << double(p * 20 + 20) / 100 << "," << time << endl;
	}
	edgeresult.close();
}

void scalabilityTest_kcore(int ds) {
	string noderesult_location = "../experiments/" + DATASET[ds] + "kcore_node_scale.txt";
	string edgeresult_location = "../experiments/" + DATASET[ds] + "kcore_edge_scale.txt";
	ofstream noderesult(noderesult_location);
	ofstream edgeresult(edgeresult_location);
	// node test
	noderesult << "sampling_rate" << "," << "time" << endl;
	for (int p = 0; p < 4; p++) {
		stringstream graph_location; graph_location << "../scalability/" << DATASET[ds] << "node" << (p * 20 + 20);
		BiGraph g(graph_location.str());
		cout << "start " << ALIAS[ds] << " node scale " << (p * 20 + 20) << " index_kcore construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexKCore(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		noderesult << double(p * 20 + 20) / 100 << "," << time << endl;
	}
	noderesult.close();
	// edge test
	edgeresult << "sampling_rate" << "," << "time" << endl;
	for (int p = 0; p < 4; p++) {
		stringstream graph_location; graph_location << "../scalability/" << DATASET[ds] << "edge" << (p * 20 + 20);
		BiGraph g(graph_location.str());
		cout << "start " << ALIAS[ds] << " edge scale " << (p * 20 + 20) << " index_kcore construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexKCore(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		edgeresult << double(p * 20 + 20) / 100 << "," << time << endl;
	}
	edgeresult.close();
}

void alpha_copy_peel_vs_alpha_peel() {
	string result_location = "../experiments/apVSacp.txt";
	ofstream result(result_location);
	result << "dataset" << "," << "alias" << "," << "alphaPeel_DELNUM" << "," << "aptime" << "," << "alphaCopyPeel_DELNUM" << "," << "acptime" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		string graph_location = INPUT_LOC + DATASET[i];
		BiGraph g(graph_location);
		time_t start, end, time;
		cout << "start " << ALIAS[i] << " index_kcore construction!" << endl;
		start = clock();
		DELETECOUNT = 0;
		coreIndexKCore_with_alpha_peel_and_cross_updating_del_count(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		result << DATASET[i] << "," << ALIAS[i] << "," << DELETECOUNT << "," << time << ",";
		g.left_index.clear(); g.right_index.clear();
		cout << "start " << ALIAS[i] << " index_kcore construction!" << endl;
		start = clock();
		DELETECOUNT = 0;
		coreIndexKCore_del_count(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		result << DELETECOUNT << "," << time << endl;
	}
	result.close();
}

void coresize() {
	for (int i = 0; i < DATASET.size(); i++) {
		string graph_location = INPUT_LOC + DATASET[i];
		BiGraph g(graph_location);
		time_t start, end, time;
		cout << "start " << ALIAS[i] << " index_kcore construction!" << endl;
		start = clock();
		int km = coreIndexKCore(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		string result_location = "../experiments/" + DATASET[i] + "coresize.txt";
		ofstream result(result_location);
		result << "alpha" << "," << "beta" << "," << "U" << "," << "V" << "," << "U+V" << endl;
		for (int a = 1; a <= km; a++) {
			for (int b = 1; b <= km; b++) {
				int left_count = 0;
				int right_count = 0;
				for (vid_t u = 0; u < g.num_v1; u++) {
					if (g.degree_v1[u] < a) continue;
					if (g.left_index[u][a] >= b) left_count++;
				}
				for (vid_t v = 0; v < g.num_v2; v++) {
					if (g.degree_v2[v] < b) continue;
					if (g.right_index[v][b] >= a) right_count++;
				}
				result << a << "," << b << "," << left_count << "," << right_count << "," << left_count + right_count << endl;
			}
		}
		result.close();
	}
}

void coresize_cutoff(int ds, int cut_begin, int cut_end) {
	string input_location = "../experiments/" + DATASET[ds] + "coresize.txt";
	ifstream input; input.open(input_location);
	string output_location = "../experiments/" + DATASET[ds] + "coresize_cutoff.txt";
	ofstream output(output_location);
	string line; // header
	input >> line;
	output << line << endl;
	while (input >> line) {
		stringstream ss(line);
		int alpha, beta;
		ss >> alpha; ss.ignore(); ss >> beta;
		if (cut_begin <= alpha && cut_end >= alpha && cut_begin <= beta && cut_end >= beta) {
			output << line << endl;
		}
	}
	input.close();
	output.close();
}

void peeltimes() {
	string output_location = OUTPUT_LOC + "peeltimes.txt";
	ofstream output(output_location);
	output << "dataset" << "," << "alias" << "," << "basic" << "," << "swap" << "," << "kcore" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		string input_location = INPUT_LOC + DATASET[i];
		BiGraph g(input_location);
		if (g.v1_max_degree < g.v2_max_degree) g.v1_max_degree = g.v2_max_degree;
		output << DATASET[i] << "," << ALIAS[i] << "," << g.v1_max_degree << "," << (ceil(sqrt(g.num_edges)) - 1)<< "," << endl;
	}
	output.close();
}

// node 1e7 with avg 2 to 10
void generate_ER() {
	for (int avg = 2; avg <= 10; avg += 2) {
		cout << "generating ER avg: " << avg << endl;
		generate_Erdos_Renyi_random_graph(1e7, 1e7, avg);
	}
}

// node 1e7 with avg 2 to 10
void generate_BA() {
	for (int avg = 2; avg <= 10; avg += 2) {
		cout << "generating BA avg: " << avg << endl;
		generate_Barabasi_Albert_random_graph(1e7, 1e7, avg);
	}
}

void vary_average_deg_basic() {
	//string ERsult_location = "../experiments/ER_avg_deg_basic.txt";
	string BAsult_location = "../experiments/BA_avg_deg_basic.txt";
	//ofstream ERresult(ERsult_location);
	ofstream BAresult(BAsult_location);
	// ER test
	/*ERresult << "avg_degree" << "," << "time" << endl;
	for (int avg = 2; avg <= 10; avg += 2) {
		stringstream graph_location; graph_location << "../randomGraph/ER_" << "l_" << 10000000 << "_r_" << 10000000 << "_avg_" << avg;
		BiGraph g(graph_location.str());
		cout << "start Erdos Renyi average degree index_basic construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexBasic(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		ERresult << avg << "," << time << endl;
	}
	ERresult.close();*/
	// BA test
	BAresult << "avg_degree" << "," << "time" << endl;
	for (int avg = 2; avg <= 10; avg += 2) {
		stringstream graph_location; graph_location << "../randomGraph/BA_" << "l_" << 10000000 << "_r_" << 10000000 << "_avg_" << avg;
		BiGraph g(graph_location.str());
		cout << "start Barabasi Albert average degree index_basic construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexBasic(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		BAresult << avg << "," << time << endl;
	}
	BAresult.close();
}

void vary_average_deg_swap() {
	//string ERsult_location = "../experiments/ER_avg_deg_swap.txt";
	string BAsult_location = "../experiments/BA_avg_deg_swap.txt";
	//ofstream ERresult(ERsult_location);
	ofstream BAresult(BAsult_location);
	// ER test
	/*ERresult << "avg_degree" << "," << "time" << endl;
	for (int avg = 2; avg <= 10; avg += 2) {
		stringstream graph_location; graph_location << "../randomGraph/ER_" << "l_" << 10000000 << "_r_" << 10000000 << "_avg_" << avg;
		BiGraph g(graph_location.str());
		cout << "start Erdos Renyi average degree index_swap construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexSwap(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		ERresult << avg << "," << time << endl;
	}
	ERresult.close();*/
	// BA test
	BAresult << "avg_degree" << "," << "time" << endl;
	for (int avg = 2; avg <= 10; avg += 2) {
		stringstream graph_location; graph_location << "../randomGraph/BA_" << "l_" << 10000000 << "_r_" << 10000000 << "_avg_" << avg;
		BiGraph g(graph_location.str());
		cout << "start Barabasi Albert average degree index_swap construction!" << endl;
		time_t start, end, time;
		start = clock();
		coreIndexSwap(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		BAresult << avg << "," << time << endl;
	}
	BAresult.close();
}

void vary_average_deg_kcore() {
	//string ERsult_location = "../experiments/ER_avg_deg_kcore.txt";
	string BAsult_location = "../experiments/BA_avg_deg_kcore.txt";
	//ofstream ERresult(ERsult_location);
	ofstream BAresult(BAsult_location);
	// ER test
	/*ERresult << "avg_degree" << "," << "time" << "," << "km" << endl;
	for (int avg = 2; avg <= 10; avg += 2) {
		stringstream graph_location; graph_location << "../randomGraph/ER_" << "l_" << 10000000 << "_r_" << 10000000 << "_avg_" << avg;
		BiGraph g(graph_location.str());
		cout << "start Erdos Renyi average degree index_kcore construction!" << endl;
		time_t start, end, time;
		start = clock();
		int km = coreIndexKCore(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		ERresult << avg << "," << time << "," << km << endl;
	}
	ERresult.close();*/
	// BA test
	BAresult << "avg_degree" << "," << "time" << "," << "km" << endl;
	for (int avg = 2; avg <= 10; avg += 2) {
		stringstream graph_location; graph_location << "../randomGraph/BA_" << "l_" << 10000000 << "_r_" << 10000000 << "_avg_" << avg;
		BiGraph g(graph_location.str());
		cout << "start Barabasi Albert average degree index_kcore construction!" << endl;
		time_t start, end, time;
		start = clock();
		int km = coreIndexKCore(g);
		end = clock();
		time = end - start;
		cout << "run time: " << time << endl;
		BAresult << avg << "," << time << "," << km << endl;
	}
	BAresult.close();
}

void bicore_index_vs_peeling() {
	vector<int> cases = { 12, 14 };
	for (int i = 0; i < cases.size(); i++) {	
		string input_location = INPUT_LOC + DATASET[cases[i]];
		BiGraph g(input_location);
		time_t start, end, time;
		cout << "start " << ALIAS[cases[i]] << " index_kcore construction!" << endl;
		start = clock();
		int km = coreIndexKCore(g);
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
		vector<int> standard = {10, 100, 400};
		for (int j = 0; j < standard.size(); j++) {
			int alpha = standard[j];
			stringstream output_location; output_location << OUTPUT_LOC << DATASET[cases[i]] << "bciVSpeel_alpha" << alpha << ".txt";
			ofstream output(output_location.str());
			output <<  "beta" << "," << "peeling" << "," << "naiveindex" << "," << "bicoreIndex" << endl;
			for (int beta = 40; beta <= 400; beta += 40) {
				output << beta << ",";
				vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2; vector<bool> left3; vector<bool> right3;
				left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2); left3.resize(g.num_v1); right3.resize(g.num_v2);
				start = clock();
				findabcore_peeling(alpha, beta, g, left1, right1);
				end = clock();
				time = end - start;
				output << time << ",";
				start = clock();
				findabcore_core_index(alpha, beta, g, left2, right2);
				end = clock();
				time = end - start;
				output << time << ",";
				start = clock();
				retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left3, right3, alpha, beta);
				end = clock();
				time = end - start;
				output << time << endl;
				if (!test_abcore_equal(left1, right1, left3, right3)) {
					cout << "error" << endl;
					return;
				}
			}
		}
		for (int j = 0; j < standard.size(); j++) {
			int beta = standard[j];
			stringstream output_location; output_location << OUTPUT_LOC << DATASET[cases[i]] << "bciVSpeel_beta" << beta << ".txt";
			ofstream output(output_location.str());
			output << "alpha" << "," << "peeling" << "," << "naiveindex" << "," << "bicoreIndex" << endl;
			for (int alpha = 40; alpha <= 400; alpha += 40) {
				output << alpha << ",";
				vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2; vector<bool> left3; vector<bool> right3;
				left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2); left3.resize(g.num_v1); right3.resize(g.num_v2);
				start = clock();
				findabcore_peeling(alpha, beta, g, left1, right1);
				end = clock();
				time = end - start;
				output << time << ",";
				start = clock();
				findabcore_core_index(alpha, beta, g, left2, right2);
				end = clock();
				time = end - start;
				output << time << ",";
				start = clock();
				retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left3, right3, alpha, beta);
				end = clock();
				time = end - start;
				output << time << endl;
				if (!test_abcore_equal(left1, right1, left3, right3)) {
					cout << "error" << endl;
					return;
				}
			}
		}
	}
}

// consider orkuts and wiki
void core_index_vs_peeling() {
	vector<int> cases = { 12, 14 };
	for (int i = 0; i < cases.size(); i++) {
		string output_location = OUTPUT_LOC + DATASET[cases[i]] + "ciVSpeel.txt";
		ofstream output(output_location);
		output << "alpha" <<  "," << "beta" << "," << "peeling" << "," << "bicoreIndex" << endl;
		string input_location = INPUT_LOC + DATASET[cases[i]];
		BiGraph g(input_location);
		time_t start, end, time;
		cout << "start " << ALIAS[cases[i]] << " index_kcore construction!" << endl;
		start = clock();
		int km = coreIndexKCore(g);
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
		for (int alpha = 2; alpha <= 10; alpha++) {
			for (int beta = 2; beta <= 10; beta++) {
				output << alpha << "," << beta << ",";
				vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2;
				left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2);
				start = clock();
				findabcore_peeling(alpha, beta, g, left1, right1);
				end = clock();
				time = end - start;
				output << time << ",";
				start = clock();
				retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left2, right2, alpha, beta);
				end = clock();
				time = end - start;
				output << time << endl;
				if (!test_abcore_equal(left1, right1, left2, right2)) {
					cout << "error" << endl;
					return;
				}
			}
		}
		output.close();
		// test for larger alpha and beta
		if (cases[i] == 14) {			
			string output_location = OUTPUT_LOC + DATASET[14] + "ciVSpeel_larger.txt";
			ofstream output(output_location);
			output << "alpha" << "," << "beta" << "," << "peeling" << "," << "bicoreIndex" << endl;
			time_t start, end, time;
			for (int alpha = 401; alpha <= 410; alpha++) {
				for (int beta = 401; beta <= 410; beta++) {
					output << alpha << "," << beta << ",";
					vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2;
					left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2);
					start = clock();
					findabcore_peeling(alpha, beta, g, left1, right1);
					end = clock();
					time = end - start;
					output << time << ",";
					start = clock();
					retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left2, right2, alpha, beta);
					end = clock();
					time = end - start;
					output << time << endl;
					if (!test_abcore_equal(left1, right1, left2, right2)) {
						cout << "error" << endl;
						return;
					}
				}
			}
			output.close();
		}
	}
}

void graph_stas() {
	string result_location = OUTPUT_LOC + "graph_stas.txt";
	ofstream result(result_location);
	result << "dataset" << "," << "alias" << "," << "U_dm" << "," << "V_dm" << "," << "U_dm+V_dm" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		string graph_location = INPUT_LOC + DATASET[i];
		BiGraph g(graph_location);
		result << DATASET[i] << "," << ALIAS[i] << "," << g.v1_max_degree << "," << g.v2_max_degree << "," << g.v1_max_degree + g.v2_max_degree << endl;
	}
	result.close();
}

void memory_size_and_retrieve_test() {
	//vector<int> alpha_fix = {5,5,5,5,18,36,12,7,38,16,8,135,360,90,360};
	//vector<int> beta_fix = {5,5,5,5,18,36,12,7,38,16,8,135,360,90,360};
	//vector<int> start_pos = {1,1,1,1,2,4,1,1,4,2,1,15,40,10,40};
	//vector<int> gap = {1,1,1,1,2,4,1,1,4,2,1,15,40,10,40};
	vector<int> alpha_fix = { 140,140,180,500,180 };
	vector<int> beta_fix = { 140,140,180,500,180 };
	vector<int> start_pos = { 15,15,20,50,20 };
	vector<int> gap = { 15,15,20,50,20 };
	//vector<int> alpha_fix = { 10,10,10,10,10,10,10,10,10,10,10,10,10,10 };
	//vector<int> beta_fix = { 10,10,10,10,10,10,10,10,10,10,10,10,10,10 };
	//vector<int> start_pos = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
	//vector<int> gap = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
	stringstream memory_size; memory_size << OUTPUT_LOC << "graph_and_index_size.txt";
	ofstream memory_summary(memory_size.str());
	memory_summary << "dataset" << "," << "graph" << "," << "index" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		string input_location = INPUT_LOC + DATASET[i];
		BiGraph g(input_location);
		time_t start, end, time;
		cout << "start " << ALIAS[i] << " index_kcore construction!" << endl;
		start = clock();
		int km = coreIndexKCore(g);
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
		// compute size
		long long int graph_size = 0; long long int index_size = 0;
		index_size += g.num_edges * 2; index_size *= sizeof(int);
		for (int cc = 0; cc < bicore_index_u.size(); cc++) {
			index_size += (bicore_index_u[cc].size() * 8);
		}
		for (int cc = 0; cc < bicore_index_v.size(); cc++) {
			index_size += (bicore_index_v[cc].size() * 8);
		}
		graph_size += g.num_v1 + g.num_v2 + g.num_edges * 2; graph_size *= sizeof(int);
		graph_size += (g.num_v1 + g.num_v2) * 4;
		memory_summary << ALIAS[i] << "," << graph_size << "," << index_size << endl;


		for (int alpha = alpha_fix[i]; alpha < 1000; alpha += 1000) {
			stringstream output_location; output_location << OUTPUT_LOC << DATASET[i] << "bciVSpeel_alpha" << alpha << ".txt";
			ofstream output(output_location.str());
			output << "beta" << "," << "peeling" << "," << "naiveindex" << "," << "bicoreIndex" << endl;
			for (int beta = start_pos[i]; beta <= start_pos[i] + gap[i] * 9; beta += gap[i]) {
				output << beta << ",";
				vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2; vector<bool> left3; vector<bool> right3;
				left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2); left3.resize(g.num_v1); right3.resize(g.num_v2);
				start = clock();
				findabcore_peeling(alpha, beta, g, left1, right1);
				end = clock();
				time = end - start;
				output << time << ",";
				start = clock();
				findabcore_core_index(alpha, beta, g, left2, right2);
				end = clock();
				time = end - start;
				output << time << ",";
				start = clock();
				retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left3, right3, alpha, beta);
				end = clock();
				time = end - start;
				output << time << endl;
				if (!test_abcore_equal(left1, right1, left3, right3)) {
					cout << "error" << endl;
					return;
				}
			}
			output.close();
		}
		for (int beta = beta_fix[i]; beta < 1000; beta += 1000) {
			stringstream output_location; output_location << OUTPUT_LOC << DATASET[i] << "bciVSpeel_beta" << beta << ".txt";
			ofstream output(output_location.str());
			output << "alpha" << "," << "peeling" << "," << "naiveindex" << "," << "bicoreIndex" << endl;
			for (int alpha = start_pos[i]; alpha <= start_pos[i] + gap[i] * 9; alpha += gap[i]) {
				output << alpha << ",";
				vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2; vector<bool> left3; vector<bool> right3;
				left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2); left3.resize(g.num_v1); right3.resize(g.num_v2);
				start = clock();
				findabcore_peeling(alpha, beta, g, left1, right1);
				end = clock();
				time = end - start;
				output << time << ",";
				start = clock();
				findabcore_core_index(alpha, beta, g, left2, right2);
				end = clock();
				time = end - start;
				output << time << ",";
				start = clock();
				retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left3, right3, alpha, beta);
				end = clock();
				time = end - start;
				output << time << endl;
				if (!test_abcore_equal(left1, right1, left3, right3)) {
					cout << "error" << endl;
					return;
				}
			}
			output.close();
		}
	}
	memory_summary.close();
}

void memory_size_and_all_dataset() {
	stringstream memory_size; memory_size << OUTPUT_LOC << "graph_and_index_size.txt";
	ofstream memory_summary(memory_size.str());
	memory_summary << "dataset" << "," << "graph" << "," << "index" << "," << "index_space_saver" << "," << "index_space_saver_dual_pointer" << endl;
	stringstream querytime; querytime << OUTPUT_LOC << "querytime_all_a10b10.txt";
	ofstream querytime_all(querytime.str());
	querytime_all << "dataset" << "," << "Peel" << "," << "Base" << "," << "Bicore" << "," << "Bicore_space_saver" << "," << "Bicore_space_saver_dual_pointer" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		string input_location = INPUT_LOC + DATASET[i];
		BiGraph g(input_location);
		time_t start, end, time;
		cout << "start " << ALIAS[i] << " index_kcore construction!" << endl;
		start = clock();
		int km = coreIndexKCore(g);
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

		vector<vector<bicore_index_block*>> bicore_index_u_s; vector<vector<bicore_index_block*>> bicore_index_v_s;
		cout << "start bicore index space saver construction!" << endl;
		start = clock();
		//build_bicore_index_space_saver(g, bicore_index_u_s, bicore_index_v_s);
		end = clock();
		time = end - start;
		cout << "buidling time: " << time << endl;

		vector<vector<bicore_index_block_dual_pointer*>> bicore_index_u_s_d; vector<vector<bicore_index_block_dual_pointer*>> bicore_index_v_s_d;
		cout << "start bicore index space saver dual pointer construction!" << endl;
		start = clock();
		//build_bicore_index_space_saver_dual_pointer(g, bicore_index_u_s_d, bicore_index_v_s_d);
		end = clock();
		time = end - start;
		cout << "buidling time: " << time << endl;
		// compute size
		long long int graph_size = 0; long long int index_size = 0; long long int index_size_space_saver = 0; long long int index_size_space_saver_dual_pointer = 0;
		index_size += g.num_edges * 2; index_size *= sizeof(int);
		for (int cc = 0; cc < bicore_index_u.size(); cc++) {
			index_size += (bicore_index_u[cc].size() * 8);
		}
		for (int cc = 0; cc < bicore_index_v.size(); cc++) {
			index_size += (bicore_index_v[cc].size() * 8);
		}
		/*for (int cc = 0; cc < bicore_index_u_s.size(); cc++) {
			index_size_space_saver += (bicore_index_u_s[cc].size() * 8);
			index_size_space_saver_dual_pointer += (bicore_index_u_s[cc].size() * 2 * 8);
			for (int ccc = 1; ccc < bicore_index_u_s[cc].size(); ccc++) {
				index_size_space_saver += bicore_index_u_s[cc][ccc]->nodeset.size() * sizeof(int);
				index_size_space_saver_dual_pointer += bicore_index_u_s[cc][ccc]->nodeset.size() * sizeof(int);
			}
		}
		for (int cc = 0; cc < bicore_index_v_s.size(); cc++) {
			index_size_space_saver += (bicore_index_v_s[cc].size() * 8);
			index_size_space_saver_dual_pointer += (bicore_index_v_s[cc].size() * 2 * 8);
			for (int ccc = 1; ccc < bicore_index_v_s[cc].size(); ccc++) {
				index_size_space_saver += bicore_index_v_s[cc][ccc]->nodeset.size() * sizeof(int);
				index_size_space_saver_dual_pointer += bicore_index_v_s[cc][ccc]->nodeset.size() * sizeof(int);
			}
		}*/
		graph_size += g.num_v1 + g.num_v2 + g.num_edges * 2; graph_size *= sizeof(int);
		graph_size += (g.num_v1 + g.num_v2) * 4;
		memory_summary << ALIAS[i] << "," << graph_size << "," << index_size << "," << index_size_space_saver << "," << index_size_space_saver_dual_pointer << endl;
		querytime_all << ALIAS[i] << ",";
		// query
		int alpha = 10; int beta = 10;
		vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2; vector<bool> left3; vector<bool> right3;
		vector<bool> left4; vector<bool> right4; vector<bool> left5; vector<bool> right5;
		left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2); left3.resize(g.num_v1); right3.resize(g.num_v2);
		left4.resize(g.num_v1); right4.resize(g.num_v2); left5.resize(g.num_v1); right5.resize(g.num_v2);
		start = clock();
		findabcore_peeling(alpha, beta, g, left1, right1);
		end = clock();
		time = end - start;
		querytime_all << time << ",";
		start = clock();
		findabcore_core_index(alpha, beta, g, left2, right2);
		end = clock();
		time = end - start;
		querytime_all << time << ",";
		start = clock();
		retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left3, right3, alpha, beta);
		end = clock();
		time = end - start;
		querytime_all << time << ",";

		start = clock();
		retrieve_via_bicore_index_space_saver(g, bicore_index_u_s, bicore_index_v_s, left4, right4, alpha, beta);
		end = clock();
		time = end - start;
		querytime_all << time << ",";
		start = clock();
		retrieve_via_bicore_index_space_saver_dual_pointer(g, bicore_index_u_s_d, bicore_index_v_s_d, left5, right5, alpha, beta);
		end = clock();
		time = end - start;
		querytime_all << time << endl;
		if (!test_abcore_equal(left1, right1, left3, right3)) {
			cout << "error" << endl;
			return;
		}
	}
	memory_summary.close();
	querytime_all.close();
}

void index_detail_print(string GRAPH_NAME) {
	string graph_location = "../realGraph/" + GRAPH_NAME;
	BiGraph g(graph_location);
	BiGraph copy1 = g;
	BiGraph copy2 = g;
	cout << "BasicCon:" << endl;
	coreIndexBasic(g);
	cout << "FPCon:" << endl;
	coreIndexSwap(g);
	cout << "ADCon:" << endl;
	coreIndexKCore(g);
}

///////////////////////////////
// dynamic part
//////////////////////////////



/*void graph_stas() {
	string result_location = OUTPUT_LOC + "graph_stas.txt";
	ofstream result(result_location);
	result << "dataset" << "," << "alias" << "," << "U_dm" << "," << "V_dm" << "," << "U_dm+V_dm" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		string graph_location = INPUT_LOC + DATASET[i];
		BiGraph g(graph_location);
		result << DATASET[i] << "," << ALIAS[i] << "," << g.v1_max_degree << "," << g.v2_max_degree << "," << g.v1_max_degree + g.v2_max_degree << endl;
	}
	result.close();
}*/

void generate_edge_streams() {
	// generate edge streams for all dataset
	for (int i = 0; i < DATASET.size(); i++) {
		stringstream input_location;
		input_location << INPUT_LOC << DATASET[i];
		stringstream output_location;
		output_location << "../dynamic/" << DATASET[i] << "addition_" << DYNAMIC_TEST_NUM << ".txt";
		generate_random_addition_stream_for_real_graph(input_location.str(), output_location.str(), DYNAMIC_TEST_NUM);
		output_location.str("");
		output_location << "../dynamic/" << DATASET[i] << "deletion_" << DYNAMIC_TEST_NUM << ".txt";
		generate_random_deletion_stream_for_real_graph(input_location.str(), output_location.str(), DYNAMIC_TEST_NUM);
	}
	// generate edge streams for scalability
	for (int i = 0; i < scala_dataset.size(); i++) {
		for (int p = 20; p <= 100; p += 20) {
			stringstream input_location; stringstream output_location;
			// node scalability
			input_location << "../scalability/" << dataset2name[scala_dataset[i]] << "node" << p;
			output_location << "../dynamic/" << dataset2name[scala_dataset[i]] << "node" << p << "_addition_" << DYNAMIC_TEST_NUM << ".txt";
			generate_random_addition_stream_for_real_graph(input_location.str(), output_location.str(), DYNAMIC_TEST_NUM);
			output_location.str("");
			output_location << "../dynamic/" << dataset2name[scala_dataset[i]] << "node" << p << "_deletion_" << DYNAMIC_TEST_NUM << ".txt";
			generate_random_deletion_stream_for_real_graph(input_location.str(), output_location.str(), DYNAMIC_TEST_NUM);
			input_location.str(""); output_location.str("");
			// edge scalability
			input_location << "../scalability/" << dataset2name[scala_dataset[i]] << "edge" << p;
			output_location << "../dynamic/" << dataset2name[scala_dataset[i]] << "edge" << p << "_addition_" << DYNAMIC_TEST_NUM << ".txt";
			generate_random_addition_stream_for_real_graph(input_location.str(), output_location.str(), DYNAMIC_TEST_NUM);
			output_location.str("");
			output_location << "../dynamic/" << dataset2name[scala_dataset[i]] << "edge" << p << "_deletion_" << DYNAMIC_TEST_NUM << ".txt";
			generate_random_deletion_stream_for_real_graph(input_location.str(), output_location.str(), DYNAMIC_TEST_NUM);
		}
	}
}

void dynamic_test_all_addition() {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/dynamic_all_addition.txt";
	ofstream result; result.open(output_location.str());
	result << "dataset" << "," << "alias" << "," << "time" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		cout << "start " << ALIAS[i] << " dynmaic addition test" << endl;
		stringstream graph_location, stream_location;
		graph_location << "../realGraph/" << DATASET[i];
		stream_location << "../dynamic/" << DATASET[i] << "addition_" << DYNAMIC_TEST_NUM << ".txt";
		time_t start, end, time;
		BiGraph g(graph_location.str());
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		ifstream edge_stream;
		edge_stream.open(stream_location.str());
		int u, v, type;
		if (DYNAMIC_WITH_LIMIT) {
			start = clock();
			while (edge_stream >> u >> v >> type) {
				if (type == 1) {
					update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, true);
				}
				else {
					update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, false);
				}
			}
			end = clock();
		}
		else {
			start = clock();
			while (edge_stream >> u >> v >> type) {
				if (type == 1) {
					update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, true);
				}
				else {
					update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, false);
				}
			}
			end = clock();
		}		
		time = end - start; // this is the total time cost of DYNAMIC_TEST_NUM operations
		cout << "time cost: " << time << endl;
		result << DATASET[i] << "," << ALIAS[i] << "," << time << endl;
	}
}

void dynamic_test_all_deletion() {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/dynamic_all_deletion.txt";
	ofstream result; result.open(output_location.str());
	result << "dataset" << "," << "alias" << "," << "time" << endl;
	for (int i = 0; i < DATASET.size(); i++) {
		cout << "start " << ALIAS[i] << " dynmaic deletion test" << endl;
		stringstream graph_location, stream_location;
		graph_location << "../realGraph/" << DATASET[i];
		stream_location << "../dynamic/" << DATASET[i] << "deletion_" << DYNAMIC_TEST_NUM << ".txt";
		time_t start, end, time;
		BiGraph g(graph_location.str());
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		ifstream edge_stream;
		edge_stream.open(stream_location.str());
		int u, v, type;
		if (DYNAMIC_WITH_LIMIT) {
			start = clock();
			while (edge_stream >> u >> v >> type) {
				if (type == 1) {
					update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, true);
				}
				else {
					update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, false);
				}
			}
			end = clock();
		}
		else {
			start = clock();
			while (edge_stream >> u >> v >> type) {
				if (type == 1) {
					update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, true);
				}
				else {
					update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, false);
				}
			}
			end = clock();
		}
		time = end - start; // this is the total time cost of DYNAMIC_TEST_NUM operations
		cout << "time cost: " << time << endl;
		result << DATASET[i] << "," << ALIAS[i] << "," << time << endl;
	}
}

void dynamic_test_scala_vary_node_addition() {
	for (int i = 0; i < scala_dataset.size(); i++) {
		stringstream output_location;
		output_location << "../experiments/dynamic_exp/dynamic_scala_node_" << dataset2name[scala_dataset[i]] << "addition.txt";
		ofstream result; result.open(output_location.str());
		result << "sampling_rate" << "," << "time" << endl;
		for (int p = 20; p <= 100; p += 20) {
			cout << "start " << scala_dataset[i] << " dynmaic addition scalability node test" << endl;
			stringstream graph_location, stream_location;
			// node scalability
			graph_location << "../scalability/" << dataset2name[scala_dataset[i]] << "node" << p;
			stream_location << "../dynamic/" << dataset2name[scala_dataset[i]] << "node" << p << "_addition_" << DYNAMIC_TEST_NUM << ".txt";
			time_t start, end, time;
			BiGraph g(graph_location.str());
			coreIndexKCore(g);
			vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
			build_bicore_index(g, bicore_index_u, bicore_index_v);
			ifstream edge_stream;
			edge_stream.open(stream_location.str());
			int u, v, type;
			if (DYNAMIC_WITH_LIMIT) {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			else {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			time = end - start; // this is the total time cost of DYNAMIC_TEST_NUM operations
			cout << "time cost: " << time << endl;
			result << double(p) / 100 << "," << time << endl;
		}
	}
}

void dynamic_test_scala_vary_node_deletion() {
	for (int i = 0; i < scala_dataset.size(); i++) {
		stringstream output_location;
		output_location << "../experiments/dynamic_exp/dynamic_scala_node_" << dataset2name[scala_dataset[i]] << "deletion.txt";
		ofstream result; result.open(output_location.str());
		result << "sampling_rate" << "," << "time" << endl;
		for (int p = 20; p <= 100; p += 20) {
			cout << "start " << scala_dataset[i] << " dynmaic deletion scalability node test" << endl;
			stringstream graph_location, stream_location;
			// node scalability
			graph_location << "../scalability/" << dataset2name[scala_dataset[i]] << "node" << p;
			stream_location << "../dynamic/" << dataset2name[scala_dataset[i]] << "node" << p << "_deletion_" << DYNAMIC_TEST_NUM << ".txt";
			time_t start, end, time;
			BiGraph g(graph_location.str());
			coreIndexKCore(g);
			vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
			build_bicore_index(g, bicore_index_u, bicore_index_v);
			ifstream edge_stream;
			edge_stream.open(stream_location.str());
			int u, v, type;
			if (DYNAMIC_WITH_LIMIT) {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			else {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			time = end - start; // this is the total time cost of DYNAMIC_TEST_NUM operations
			cout << "time cost: " << time << endl;
			result << double(p) / 100 << "," << time << endl;
		}
	}
}

void dynamic_test_scala_vary_edge_addition() {
	for (int i = 0; i < scala_dataset.size(); i++) {
		stringstream output_location;
		output_location << "../experiments/dynamic_exp/dynamic_scala_edge_" << dataset2name[scala_dataset[i]] << "addition.txt";
		ofstream result; result.open(output_location.str());
		result << "sampling_rate" << "," << "time" << endl;
		for (int p = 20; p <= 100; p += 20) {
			cout << "start " << scala_dataset[i] << " dynmaic addition scalability edge test" << endl;
			stringstream graph_location, stream_location;
			// node scalability
			graph_location << "../scalability/" << dataset2name[scala_dataset[i]] << "edge" << p;
			stream_location << "../dynamic/" << dataset2name[scala_dataset[i]] << "edge" << p << "_addition_" << DYNAMIC_TEST_NUM << ".txt";
			time_t start, end, time;
			BiGraph g(graph_location.str());
			coreIndexKCore(g);
			vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
			build_bicore_index(g, bicore_index_u, bicore_index_v);
			ifstream edge_stream;
			edge_stream.open(stream_location.str());
			int u, v, type;
			if (DYNAMIC_WITH_LIMIT) {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			else {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			time = end - start; // this is the total time cost of DYNAMIC_TEST_NUM operations
			cout << "time cost: " << time << endl;
			result << double(p) / 100 << "," << time << endl;
		}
	}
}

void dynamic_test_scala_vary_edge_deletion() {
	for (int i = 0; i < scala_dataset.size(); i++) {
		stringstream output_location;
		output_location << "../experiments/dynamic_exp/dynamic_scala_edge_" << dataset2name[scala_dataset[i]] << "deletion.txt";
		ofstream result; result.open(output_location.str());
		result << "sampling_rate" << "," << "time" << endl;
		for (int p = 20; p <= 100; p += 20) {
			cout << "start " << scala_dataset[i] << " dynmaic deletion scalability edge test" << endl;
			stringstream graph_location, stream_location;
			// node scalability
			graph_location << "../scalability/" << dataset2name[scala_dataset[i]] << "edge" << p;
			stream_location << "../dynamic/" << dataset2name[scala_dataset[i]] << "edge" << p << "_deletion_" << DYNAMIC_TEST_NUM << ".txt";
			time_t start, end, time;
			BiGraph g(graph_location.str());
			coreIndexKCore(g);
			vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
			build_bicore_index(g, bicore_index_u, bicore_index_v);
			ifstream edge_stream;
			edge_stream.open(stream_location.str());
			int u, v, type;
			if (DYNAMIC_WITH_LIMIT) {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			else {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			time = end - start; // this is the total time cost of DYNAMIC_TEST_NUM operations
			cout << "time cost: " << time << endl;
			result << double(p) / 100 << "," << time << endl;
		}
	}
}

void dynamic_degree_out_put(string graph_location, string edge_stream_location) {
	BiGraph g(graph_location);
	int u, v, type;
	ifstream edge_stream;
	edge_stream.open(edge_stream_location);
	cout << "degree_u" << " " << "degree_v" << endl;
	while (edge_stream >> u >> v >> type) {
		cout << g.degree_v1[u] << " " << g.degree_v2[v] << endl;
	}
	edge_stream.close();
}

void dynamic_test_scala_vary_node_deletion_single_test() {
	vector<int> scala_dataset = { 12 };
	for (int i = 0; i < scala_dataset.size(); i++) {
		stringstream output_location;
		output_location << "../experiments/dynamic_exp/dynamic_scala_node_" << DATASET[scala_dataset[i]] << "deletion.txt";
		ofstream result; result.open(output_location.str());
		result << "sampling_rate" << "," << "time" << endl;
		for (int p = 20; p <= 20; p += 20) {
			cout << "start " << ALIAS[scala_dataset[i]] << " dynmaic deletion scalability node test" << endl;
			stringstream graph_location, stream_location;
			// node scalability
			graph_location << "../scalability/" << DATASET[scala_dataset[i]] << "node" << p;
			stream_location << "../dynamic/" << DATASET[scala_dataset[i]] << "node" << p << "_deletion_" << DYNAMIC_TEST_NUM << ".txt";
			time_t start, end, time;
			BiGraph g(graph_location.str());
			coreIndexKCore(g);
			vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
			build_bicore_index(g, bicore_index_u, bicore_index_v);
			ifstream edge_stream;
			edge_stream.open(stream_location.str());
			int u, v, type;
			if (DYNAMIC_WITH_LIMIT) {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_with_limit(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			else {
				start = clock();
				while (edge_stream >> u >> v >> type) {
					if (type == 1) {
						update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, true);
					}
					else {
						update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, false);
					}
				}
				end = clock();
			}
			time = end - start; // this is the total time cost of DYNAMIC_TEST_NUM operations
			cout << "time cost: " << time << endl;
			result << double(p) / 100 << "," << time << endl;
		}
	}
}

void dynamic_test_insertion_single_test_for_optimize_purpose(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/insertion_" << ds << "_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "dataset" << "," << "alias" << "," << "time(avg)" << endl;
	cout << "start " << ds << " " << iu << " dynmaic insertion test" << endl;
	stringstream graph_location, stream_location;
	graph_location << "../realGraph/" << dataset2name[ds];
	BiGraph g(graph_location.str());

	stream_location << "../dynamic/" << dataset2name[ds] << "deletion_" << DYNAMIC_TEST_NUM << ".txt";
	ifstream edge_stream_removal;
	edge_stream_removal.open(stream_location.str());
	int u, v, type;
	while (edge_stream_removal >> u >> v >> type)
	{
		if (u < g.num_v1 && v < g.num_v2) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
			int ss = tmp_neigh_.size();
			bool deleted = false;
			for (int i = 0; i < ss; i++) {
				if (tmp_neigh_[i] == v) {
					deleted = true;
					g.deleteEdge(u, v);
					break;
				}
			}
			if (!deleted) {
				cout << "illegal deletion 1" << endl;
				return;
			}
		}
		else {
			cout << "illegal deletion 2" << endl;
			return;
		}
	}
	edge_stream_removal.close();


	auto start = chrono::system_clock::now();
	coreIndexKCore(g);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << elapsed_seconds.count() << endl;
	result << "index construction time: " << elapsed_seconds.count() << endl;
	vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
	//build_bicore_index(g, bicore_index_u, bicore_index_v);
	ifstream edge_stream;
	edge_stream.open(stream_location.str());
	double time = 0;
	
	while (edge_stream >> u >> v >> type) {
		type = 1;
		if (iu == Index_update::withlimit) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_base_opt) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_dfs) {
			time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withoutlimit) {
			time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
		else if (iu == Index_update::withlimit_dfs_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
	}
	cout << "time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
	result << dataset2name[ds] << "," << ds << "," << time / DYNAMIC_TEST_NUM << endl;
}

void dynamic_test_insertion_single_test_random_edge_version(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/insertion_" << ds << "_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "dataset" << "," << "alias" << "," << "time(avg)" << endl;
	cout << "start " << ds << " " << iu << " dynmaic insertion test" << endl;
	stringstream graph_location, stream_location;
	graph_location << "../realGraph/" << dataset2name[ds];
	stream_location << "../dynamic/" << dataset2name[ds] << "addition_" << DYNAMIC_TEST_NUM << ".txt";
	BiGraph g(graph_location.str());
	auto start = chrono::system_clock::now();
	coreIndexKCore(g);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << elapsed_seconds.count() << endl;
	result << "index construction time: " << elapsed_seconds.count() << endl;
	vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
	//build_bicore_index(g, bicore_index_u, bicore_index_v);
	ifstream edge_stream;
	edge_stream.open(stream_location.str());
	double time = 0;
	int u, v, type;
	while (edge_stream >> u >> v >> type) {
		if (iu == Index_update::withlimit) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_base_opt) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_dfs) {
			time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withoutlimit) {
			time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
		else if (iu == Index_update::withlimit_dfs_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
	}
	cout << "time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
	result << dataset2name[ds] << "," << ds << "," << time / DYNAMIC_TEST_NUM << endl;
}

double dynamic_test_insertion_single_test(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/insertion_" << ds << "_" << iu << "_" << threads << ".txt";
	ofstream result; result.open(output_location.str());
	result << "dataset" << "," << "alias" << "," << "time(avg)" << endl;
	cout << "start " << ds << " " << iu << " dynmaic insertion test" << endl;
	stringstream graph_location, stream_location;
	graph_location << "../realGraph/" << dataset2name[ds];
	BiGraph g(graph_location.str());

	stream_location << "../dynamic/" << dataset2name[ds] << "deletion_" << DYNAMIC_TEST_NUM << ".txt";
	ifstream edge_stream_removal;
	edge_stream_removal.open(stream_location.str());
	int u, v, type;
	while (edge_stream_removal >> u >> v >> type)
	{
		if (u < g.num_v1 && v < g.num_v2) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
			int ss = tmp_neigh_.size();
			bool deleted = false;
			for (int i = 0; i < ss; i++) {
				if (tmp_neigh_[i] == v) {
					deleted = true;
					g.deleteEdge(u, v);
					break;
				}
			}
			if (!deleted) {
				cout << "illegal deletion 1" << endl;
				return 0;
			}
		}
		else {
			cout << "illegal deletion 2" << endl;
			return 0;
		}
	}
	edge_stream_removal.close();


	auto start = chrono::system_clock::now();
	coreIndexKCore(g);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << elapsed_seconds.count() << endl;
	result << "index construction time: " << elapsed_seconds.count() << endl;
	vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
	//build_bicore_index(g, bicore_index_u, bicore_index_v);
	ifstream edge_stream;
	edge_stream.open(stream_location.str());
	double time = 0;

	while (edge_stream >> u >> v >> type) {
		type = 1;
		if (iu == Index_update::withlimit) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_base_opt) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_dfs) {
			time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withoutlimit) {
			time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
		else if (iu == Index_update::withlimit_dfs_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
	}
	cout << "time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
	result << dataset2name[ds] << "," << ds << "," << time / DYNAMIC_TEST_NUM << endl;
	return time / DYNAMIC_TEST_NUM;
}


void dynamic_test_removal_single_test_for_optimize_purpose(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/removal_" << ds << "_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "dataset" << "," << "alias" << "," << "time(avg)" << endl;
	cout << "start " << ds << " " << iu << " dynmaic removal test" << endl;
	stringstream graph_location, stream_location;
	graph_location << "../realGraph/" << dataset2name[ds];
	BiGraph g(graph_location.str());

	stream_location << "../dynamic/" << dataset2name[ds] << "addition_" << DYNAMIC_TEST_NUM << ".txt";
	ifstream edge_stream_insert;
	edge_stream_insert.open(stream_location.str());
	int u, v, type;
	while (edge_stream_insert >> u >> v >> type) {
		if (u < g.num_v1 && v < g.num_v2) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
			int ss = tmp_neigh_.size();
			for (int i = 0; i < ss; i++) {
				if (tmp_neigh_[i] == v) {
					cout << "illegal insertion 1" << endl;
					return;
				}
			}
			g.addEdge(u, v);
		}
		else {
			cout << "illegal insertion 2" << endl;
			return;
		}
	}
	edge_stream_insert.close();

	
	auto start = chrono::system_clock::now();
	coreIndexKCore(g);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << elapsed_seconds.count() << endl;
	result << "index construction time: " << elapsed_seconds.count() << endl;
	vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
	//build_bicore_index(g, bicore_index_u, bicore_index_v);
	ifstream edge_stream;
	edge_stream.open(stream_location.str());
	double time = 0;
	while (edge_stream >> u >> v >> type) {
		type = 0;
		if (iu == Index_update::withlimit) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_base_opt) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_dfs) {
			time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withoutlimit) {
			time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
		else if (iu == Index_update::withlimit_dfs_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
	}
	cout << "time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
	result << dataset2name[ds] << "," << ds << "," << time / DYNAMIC_TEST_NUM << endl;
}

double dynamic_test_removal_single_test(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/removal_" << ds << "_" << iu << "_" << threads << ".txt";
	ofstream result; result.open(output_location.str());
	result << "dataset" << "," << "alias" << "," << "time(avg)" << endl;
	cout << "start " << ds << " " << iu <<  " dynmaic removal test" << endl;
	stringstream graph_location, stream_location;
	graph_location << "../realGraph/" << dataset2name[ds];
	stream_location << "../dynamic/" << dataset2name[ds] << "deletion_" << DYNAMIC_TEST_NUM << ".txt";
	BiGraph g(graph_location.str());
	auto start = chrono::system_clock::now();
	coreIndexKCore(g);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << elapsed_seconds.count() << endl;
	result << "index construction time: " << elapsed_seconds.count() << endl;
	vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
	//build_bicore_index(g, bicore_index_u, bicore_index_v);
	ifstream edge_stream;
	edge_stream.open(stream_location.str());
	double time = 0;
	int u, v, type;
	while (edge_stream >> u >> v >> type) {
		if (iu == Index_update::withlimit) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_base_opt) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_dfs) {
			time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withoutlimit) {
			time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
		}
		else if (iu == Index_update::withlimit_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
		else if (iu == Index_update::withlimit_dfs_parallel) {
			time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
		}
	}
	cout << "time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
	result << dataset2name[ds] << "," << ds << "," << time / DYNAMIC_TEST_NUM << endl;
	return time / DYNAMIC_TEST_NUM;
}

void dynamic_test_batch_single_test(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/batch_" << ds << "_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "dataset" << "," << "alias" << "," << "time(avg)" << endl;
	cout << "start " << ds << " " << iu << " dynmaic batch test" << endl;
	stringstream graph_location, stream_location;
	graph_location << "../realGraph/" << dataset2name[ds];
	// should be change to batch
	stream_location << "../dynamic/" << dataset2name[ds] << "deletion_" << DYNAMIC_TEST_NUM << ".txt";
	BiGraph g(graph_location.str());
	auto start = chrono::system_clock::now();
	coreIndexKCore(g);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << elapsed_seconds.count() << endl;
	result << "index construction time: " << elapsed_seconds.count() << endl;
	vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
	//build_bicore_index(g, bicore_index_u, bicore_index_v);
	ifstream edge_stream;
	edge_stream.open(stream_location.str());
	double time = 0;
	int u, v, type;
	batch_set I, R;
	while (edge_stream >> u >> v >> type) {
		if (type == 0) {
			auto e = make_pair(u, v);
			if (I.find(e) != I.end()) {
				I.erase(e);
			}
			else {
				R.insert(e);
			}
		}
		if (type == 1) {
			auto e = make_pair(u, v);
			if (R.find(e) != R.end()) {
				R.erase(e);
			}
			else {
				I.insert(e);
			}
		}
	}
	if (iu == Index_update::batch) {
		time = Dyn_rebuild::update_bicore_index_batch(g, bicore_index_u, bicore_index_v, I, R);
	}
	else {
		time = Dyn_rebuild::update_bicore_index_batch_parallel(g, bicore_index_u, bicore_index_v, I, R, threads);
	}
	cout << "time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
	result << dataset2name[ds] << "," << ds << "," << time / DYNAMIC_TEST_NUM << endl;
}

void dynamic_test_all_insertion() {
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::WC, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::WC, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::WC, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::WC, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::FG, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::FG, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::FG, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::FG, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::EP, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::EP, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::EP, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::EP, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::DE, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::DE, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::DE, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::DE, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::RE, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::RE, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::RE, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::RE, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::TR, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::TR, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::TR, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::TR, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::DUI, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::DUI, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::DUI, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::DUI, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::LJ, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::LJ, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 10;
	//dynamic_test_insertion_single_test(Datasets::LJ, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::LJ, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::WT, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::WT, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 10;
	//dynamic_test_insertion_single_test(Datasets::WT, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::WT, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_insertion_single_test(Datasets::OG, Index_update::withlimit, 12);
	dynamic_test_insertion_single_test(Datasets::OG, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 10;
	//dynamic_test_insertion_single_test(Datasets::OG, Index_update::withlimit_parallel, 12);
	dynamic_test_insertion_single_test(Datasets::OG, Index_update::withlimit_dfs_parallel, 12);
}

void dynamic_test_all_removal() {
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::WC, Index_update::withlimit, 12);
	dynamic_test_removal_single_test(Datasets::WC, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::WC, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::WC, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::FG, Index_update::withlimit, 12);
	dynamic_test_removal_single_test(Datasets::FG, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::FG, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::FG, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::EP, Index_update::withlimit, 12);
	dynamic_test_removal_single_test(Datasets::EP, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::EP, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::EP, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::DE, Index_update::withlimit, 12);
	dynamic_test_removal_single_test(Datasets::DE, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::DE, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::DE, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::RE, Index_update::withlimit, 12);
	dynamic_test_removal_single_test(Datasets::RE, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::RE, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::RE, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::TR, Index_update::withlimit, 12);
	dynamic_test_removal_single_test(Datasets::TR, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::TR, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::TR, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	//dynamic_test_removal_single_test(Datasets::DUI, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::DUI, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 10;
	//dynamic_test_removal_single_test(Datasets::LJ, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::LJ, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 10;
	//dynamic_test_removal_single_test(Datasets::WT, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::WT, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 10;
	//dynamic_test_removal_single_test(Datasets::OG, Index_update::withlimit_parallel, 12);
	dynamic_test_removal_single_test(Datasets::OG, Index_update::withlimit_dfs_parallel, 12);

	//DYNAMIC_TEST_NUM = 10;
	//dynamic_test_removal_single_test(Datasets::DUI, Index_update::withlimit, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_removal_single_test(Datasets::DUI, Index_update::withlimit_dfs, 12);
	//
	//
	//DYNAMIC_TEST_NUM = 10;
	//dynamic_test_removal_single_test(Datasets::LJ, Index_update::withlimit, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_removal_single_test(Datasets::LJ, Index_update::withlimit_dfs, 12);
	//

	//DYNAMIC_TEST_NUM = 10;
	//dynamic_test_removal_single_test(Datasets::WT, Index_update::withlimit, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_removal_single_test(Datasets::WT, Index_update::withlimit_dfs, 12);
	//

	//DYNAMIC_TEST_NUM = 10;
	//dynamic_test_removal_single_test(Datasets::OG, Index_update::withlimit, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_removal_single_test(Datasets::OG, Index_update::withlimit_dfs, 12);
}

void dynamic_test_all_batch() {
	DYNAMIC_TEST_NUM = 1000;
	dynamic_test_batch_single_test(Datasets::WC, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::FG, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::EP, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::DE, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::RE, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::TR, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::DUI, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::LJ, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::WT, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::OG, Index_update::batch, 12);

	dynamic_test_batch_single_test(Datasets::WC, Index_update::batch_parallel, 12);

	dynamic_test_batch_single_test(Datasets::FG, Index_update::batch_parallel, 12);

	dynamic_test_batch_single_test(Datasets::EP, Index_update::batch_parallel, 12);

	dynamic_test_batch_single_test(Datasets::DE, Index_update::batch_parallel, 12);

	dynamic_test_batch_single_test(Datasets::RE, Index_update::batch_parallel, 12);

	dynamic_test_batch_single_test(Datasets::TR, Index_update::batch_parallel, 12);

	dynamic_test_batch_single_test(Datasets::DUI, Index_update::batch_parallel, 12);

	dynamic_test_batch_single_test(Datasets::LJ, Index_update::batch_parallel, 12);

	dynamic_test_batch_single_test(Datasets::WT, Index_update::batch_parallel, 12);

	dynamic_test_batch_single_test(Datasets::OG, Index_update::batch_parallel, 12);
}

void dynamic_test_scala_vary_node_removal_single_test(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/scalability/dynamic_scala_node_" << ds << "_removal_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "sampling_rate" << "," << "time(avg)" << endl;
	if (iu == Index_update::withlimit_dfs_parallel || iu == Index_update::withlimit_parallel) {
		result << "numcores: " << threads << endl;
	}
	cout << "start " << ds << iu << " dynmaic removal scalability node test" << endl;
	for (int p = 20; p <= 80; p += 20) {
		stringstream graph_location, stream_location;
		// node scalability
		graph_location << "../scalability/" << dataset2name[ds] << "node" << p;
		stream_location << "../dynamic/" << dataset2name[ds] << "node" << p << "_deletion_" << DYNAMIC_TEST_NUM << ".txt";
		BiGraph g(graph_location.str());
		auto start = chrono::system_clock::now();
		coreIndexKCore(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		cout << elapsed_seconds.count() << endl;
		result << "index construction time: " << elapsed_seconds.count() << endl;
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		ifstream edge_stream;
		edge_stream.open(stream_location.str());
		double time = 0;
		int u, v, type;
		while (edge_stream >> u >> v >> type) {
			if (iu == Index_update::withlimit) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_base_opt) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_dfs) {
				time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withoutlimit) {
				time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
			else if (iu == Index_update::withlimit_dfs_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
		}
		cout << "sampling rate: " << double(p) / 100 << " time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
		result << double(p) / 100 << "," << time / DYNAMIC_TEST_NUM << endl;
	}
}

void dynamic_test_scala_vary_node_insertion_single_test_random_edge_version(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/scalability/dynamic_scala_node_" << ds << "_insertion_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "sampling_rate" << "," << "time(avg)" << endl;
	if (iu == Index_update::withlimit_dfs_parallel || iu == Index_update::withlimit_parallel) {
		result << "numcores: " << threads << endl;
	}
	cout << "start " << ds << iu << " dynmaic insertion scalability node test" << endl;
	for (int p = 20; p <= 80; p += 20) {
		stringstream graph_location, stream_location;
		// node scalability
		graph_location << "../scalability/" << dataset2name[ds] << "node" << p;
		stream_location << "../dynamic/" << dataset2name[ds] << "node" << p << "_addition_" << DYNAMIC_TEST_NUM << ".txt";
		BiGraph g(graph_location.str());
		auto start = chrono::system_clock::now();
		coreIndexKCore(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		cout << elapsed_seconds.count() << endl;
		result << "index construction time: " << elapsed_seconds.count() << endl;
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		ifstream edge_stream;
		edge_stream.open(stream_location.str());
		double time = 0;
		int u, v, type;
		while (edge_stream >> u >> v >> type) {
			if (iu == Index_update::withlimit) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_base_opt) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_dfs) {
				time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withoutlimit) {
				time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
			else if (iu == Index_update::withlimit_dfs_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
		}
		cout << "sampling rate: " << double(p) / 100 << " time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
		result << double(p) / 100 << "," << time / DYNAMIC_TEST_NUM << endl;
	}
}

void dynamic_test_scala_vary_node_insertion_single_test(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/scalability/dynamic_scala_node_" << ds << "_insertion_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "sampling_rate" << "," << "time(avg)" << endl;
	if (iu == Index_update::withlimit_dfs_parallel || iu == Index_update::withlimit_parallel) {
		result << "numcores: " << threads << endl;
	}
	cout << "start " << ds << iu << " dynmaic insertion scalability node test" << endl;
	for (int p = 20; p <= 80; p += 20) {
		stringstream graph_location, stream_location;
		// node scalability
		graph_location << "../scalability/" << dataset2name[ds] << "node" << p;
		stream_location << "../dynamic/" << dataset2name[ds] << "node" << p << "_deletion_" << DYNAMIC_TEST_NUM << ".txt";
		BiGraph g(graph_location.str());

		ifstream edge_stream_removal;
		edge_stream_removal.open(stream_location.str());
		int u, v, type;
		while (edge_stream_removal >> u >> v >> type)
		{
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				bool deleted = false;
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						deleted = true;
						g.deleteEdge(u, v);
						break;
					}
				}
				if (!deleted) {
					cout << "illegal deletion 1" << endl;
					return;
				}
			}
			else {
				cout << "illegal deletion 2" << endl;
				return;
			}
		}
		edge_stream_removal.close();

		auto start = chrono::system_clock::now();
		coreIndexKCore(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		cout << elapsed_seconds.count() << endl;
		result << "index construction time: " << elapsed_seconds.count() << endl;
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		ifstream edge_stream;
		edge_stream.open(stream_location.str());
		double time = 0;
		while (edge_stream >> u >> v >> type) {
			type = 1;
			if (iu == Index_update::withlimit) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_base_opt) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_dfs) {
				time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withoutlimit) {
				time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
			else if (iu == Index_update::withlimit_dfs_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
		}
		cout << "sampling rate: " << double(p) / 100 << " time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
		result << double(p) / 100 << "," << time / DYNAMIC_TEST_NUM << endl;
	}
}

void dynamic_test_scala_vary_edge_removal_single_test(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/scalability/dynamic_scala_edge_" << ds << "_removal_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "sampling_rate" << "," << "time(avg)" << endl;
	if (iu == Index_update::withlimit_dfs_parallel || iu == Index_update::withlimit_parallel) {
		result << "numcores: " << threads << endl;
	}
	cout << "start " << ds << iu << " dynmaic removal scalability edge test" << endl;
	for (int p = 20; p <= 80; p += 20) {
		stringstream graph_location, stream_location;
		// node scalability
		graph_location << "../scalability/" << dataset2name[ds] << "edge" << p;
		stream_location << "../dynamic/" << dataset2name[ds] << "edge" << p << "_deletion_" << DYNAMIC_TEST_NUM << ".txt";
		BiGraph g(graph_location.str());
		auto start = chrono::system_clock::now();
		coreIndexKCore(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		cout << elapsed_seconds.count() << endl;
		result << "index construction time: " << elapsed_seconds.count() << endl;
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		ifstream edge_stream;
		edge_stream.open(stream_location.str());
		double time = 0;
		int u, v, type;
		while (edge_stream >> u >> v >> type) {
			if (iu == Index_update::withlimit) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_base_opt) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_dfs) {
				time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withoutlimit) {
				time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
			else if (iu == Index_update::withlimit_dfs_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
		}
		cout << "sampling rate: " << double(p) / 100 << " time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
		result << double(p) / 100 << "," << time / DYNAMIC_TEST_NUM << endl;
	}
}

void dynamic_test_scala_vary_edge_insertion_single_test_random_edge_version(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/scalability/dynamic_scala_edge_" << ds << "_insertion_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "sampling_rate" << "," << "time(avg)" << endl;
	if (iu == Index_update::withlimit_dfs_parallel || iu == Index_update::withlimit_parallel) {
		result << "numcores: " << threads << endl;
	}
	cout << "start " << ds << iu << " dynmaic insertion scalability edge test" << endl;
	for (int p = 20; p <= 80; p += 20) {
		stringstream graph_location, stream_location;
		// node scalability
		graph_location << "../scalability/" << dataset2name[ds] << "edge" << p;
		stream_location << "../dynamic/" << dataset2name[ds] << "edge" << p << "_addition_" << DYNAMIC_TEST_NUM << ".txt";
		BiGraph g(graph_location.str());
		auto start = chrono::system_clock::now();
		coreIndexKCore(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		cout << elapsed_seconds.count() << endl;
		result << "index construction time: " << elapsed_seconds.count() << endl;
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		ifstream edge_stream;
		edge_stream.open(stream_location.str());
		double time = 0;
		int u, v, type;
		while (edge_stream >> u >> v >> type) {
			if (iu == Index_update::withlimit) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_base_opt) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_dfs) {
				time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withoutlimit) {
				time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
			else if (iu == Index_update::withlimit_dfs_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
		}
		cout << "sampling rate: " << double(p) / 100 << " time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
		result << double(p) / 100 << "," << time / DYNAMIC_TEST_NUM << endl;
	}
}

void dynamic_test_scala_vary_edge_insertion_single_test(Datasets ds, Index_update iu, int threads) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/scalability/dynamic_scala_edge_" << ds << "_insertion_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	result << "sampling_rate" << "," << "time(avg)" << endl;
	if (iu == Index_update::withlimit_dfs_parallel || iu == Index_update::withlimit_parallel) {
		result << "numcores: " << threads << endl;
	}
	cout << "start " << ds << iu << " dynmaic insertion scalability edge test" << endl;
	for (int p = 20; p <= 80; p += 20) {
		stringstream graph_location, stream_location;
		// node scalability
		graph_location << "../scalability/" << dataset2name[ds] << "edge" << p;
		stream_location << "../dynamic/" << dataset2name[ds] << "edge" << p << "_deletion_" << DYNAMIC_TEST_NUM << ".txt";
		BiGraph g(graph_location.str());

		ifstream edge_stream_removal;
		edge_stream_removal.open(stream_location.str());
		int u, v, type;
		while (edge_stream_removal >> u >> v >> type)
		{
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				bool deleted = false;
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						deleted = true;
						g.deleteEdge(u, v);
						break;
					}
				}
				if (!deleted) {
					cout << "illegal deletion 1" << endl;
					return;
				}
			}
			else {
				cout << "illegal deletion 2" << endl;
				return;
			}
		}
		edge_stream_removal.close();

		auto start = chrono::system_clock::now();
		coreIndexKCore(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		cout << elapsed_seconds.count() << endl;
		result << "index construction time: " << elapsed_seconds.count() << endl;
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		ifstream edge_stream;
		edge_stream.open(stream_location.str());
		double time = 0;
		while (edge_stream >> u >> v >> type) {
			type = 1;
			if (iu == Index_update::withlimit) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_base_opt) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_base_opt(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_dfs) {
				time += Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withoutlimit) {
				time += Dyn_rebuild::update_bicore_index_without_limit_swap(g, bicore_index_u, bicore_index_v, u, v, type);
			}
			else if (iu == Index_update::withlimit_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
			else if (iu == Index_update::withlimit_dfs_parallel) {
				time += Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, u, v, type, threads);
			}
		}
		cout << "sampling rate: " << double(p) / 100 << " time cost (avg): " << time / DYNAMIC_TEST_NUM << endl;
		result << double(p) / 100 << "," << time / DYNAMIC_TEST_NUM << endl;
	}
}

void dynamic_test_scala_vary_node_insertion() {
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::WC, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::WC, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::WC, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::WC, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::FG, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::FG, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::FG, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::FG, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::TR, Index_update::withlimit, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::TR, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::TR, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::TR, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::DUI, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::DUI, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::DUI, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::DUI, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::LJ, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::LJ, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 10;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::LJ, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::LJ, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::WT, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::WT, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 10;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::WT, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::WT, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::OG, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::OG, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 10;
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::OG, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_insertion_single_test(Datasets::OG, Index_update::withlimit_dfs_parallel, 12);
}

void dynamic_test_scala_vary_node_removal() {
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::WC, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::WC, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::WC, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::WC, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::FG, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::FG, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::FG, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::FG, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::TR, Index_update::withlimit, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::TR, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::TR, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::TR, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::DUI, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::DUI, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::DUI, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::DUI, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::LJ, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::LJ, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 10;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::LJ, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::LJ, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::WT, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::WT, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 10;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::WT, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::WT, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::OG, Index_update::withlimit, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::OG, Index_update::withlimit_dfs, 12);
	DYNAMIC_TEST_NUM = 10;
	dynamic_test_scala_vary_node_removal_single_test(Datasets::OG, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_node_removal_single_test(Datasets::OG, Index_update::withlimit_dfs_parallel, 12);
}

void dynamic_test_scala_vary_edge_insertion() {
	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::WC, Index_update::withlimit, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::WC, Index_update::withlimit_dfs, 12);
	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::WC, Index_update::withlimit_parallel, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::WC, Index_update::withlimit_dfs_parallel, 12);

	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::FG, Index_update::withlimit, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::FG, Index_update::withlimit_dfs, 12);
	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::FG, Index_update::withlimit_parallel, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::FG, Index_update::withlimit_dfs_parallel, 12);

	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::TR, Index_update::withlimit_parallel, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::TR, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_insertion_single_test(Datasets::TR, Index_update::withlimit, 12);
	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_insertion_single_test(Datasets::TR, Index_update::withlimit_dfs, 12);

	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::DUI, Index_update::withlimit, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::DUI, Index_update::withlimit_dfs, 12);
	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::DUI, Index_update::withlimit_parallel, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::DUI, Index_update::withlimit_dfs_parallel, 12);

	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::LJ, Index_update::withlimit, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::LJ, Index_update::withlimit_dfs, 12);
	//DYNAMIC_TEST_NUM = 10;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::LJ, Index_update::withlimit_parallel, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::LJ, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_insertion_single_test(Datasets::WT, Index_update::withlimit, 12);
	dynamic_test_scala_vary_edge_insertion_single_test(Datasets::WT, Index_update::withlimit_dfs, 12);
	//DYNAMIC_TEST_NUM = 10;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::WT, Index_update::withlimit_parallel, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::WT, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_insertion_single_test(Datasets::OG, Index_update::withlimit, 12);
	dynamic_test_scala_vary_edge_insertion_single_test(Datasets::OG, Index_update::withlimit_dfs, 12);

	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::OG, Index_update::withlimit_parallel, 12);
	//dynamic_test_scala_vary_edge_insertion_single_test(Datasets::OG, Index_update::withlimit_dfs_parallel, 12);


	
}

void dynamic_test_scala_vary_edge_removal() {

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::WC, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::WC, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::FG, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::FG, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::TR, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::TR, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::DUI, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::DUI, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 10;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::LJ, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::LJ, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 10;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::WT, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::WT, Index_update::withlimit_dfs_parallel, 12);

	DYNAMIC_TEST_NUM = 10;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::OG, Index_update::withlimit_parallel, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::OG, Index_update::withlimit_dfs_parallel, 12);

	/*DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::WC, Index_update::withlimit, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::WC, Index_update::withlimit_dfs, 12);
	

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::FG, Index_update::withlimit, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::FG, Index_update::withlimit_dfs, 12);*/
	

	//DYNAMIC_TEST_NUM = 10;
	//dynamic_test_scala_vary_edge_removal_single_test(Datasets::TR, Index_update::withlimit, 12);
	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_removal_single_test(Datasets::TR, Index_update::withlimit_dfs, 12);


	/*DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::DUI, Index_update::withlimit, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::DUI, Index_update::withlimit_dfs, 12);
	

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::LJ, Index_update::withlimit, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::LJ, Index_update::withlimit_dfs, 12);
	

	DYNAMIC_TEST_NUM = 100;
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::WT, Index_update::withlimit, 12);
	dynamic_test_scala_vary_edge_removal_single_test(Datasets::WT, Index_update::withlimit_dfs, 12);*/
	

	//DYNAMIC_TEST_NUM = 10;
	//dynamic_test_scala_vary_edge_removal_single_test(Datasets::OG, Index_update::withlimit, 12);
	//DYNAMIC_TEST_NUM = 100;
	//dynamic_test_scala_vary_edge_removal_single_test(Datasets::OG, Index_update::withlimit_dfs, 12);
	
}

void index_construction_scala_vary_node_single_test(Datasets ds, int threads) {
	stringstream output_location;
	output_location << "../experiments/multithread_scala_node_" << ds << ".txt";
	ofstream result; result.open(output_location.str());
	// node test
	result << "sampling_rate" << "," << "time" << endl;
	for (int p = 0; p < 4; p++) {
		stringstream graph_location; graph_location << "../scalability/" << dataset2name[ds] << "node" << (p * 20 + 20);
		BiGraph g(graph_location.str());
		cout << "start " << ds << " node scale " << (p * 20 + 20) << " multithread index_kcore construction!" << endl;
		auto start = chrono::system_clock::now();
		multithread_index_construction(g, threads);
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		cout << elapsed_seconds.count() << endl;
		result << double(p * 20 + 20) / 100 << "," << elapsed_seconds.count() << endl;
	}
	stringstream graph_location; graph_location << "../realGraph/" << dataset2name[ds];
	BiGraph g(graph_location.str());
	cout << "start " << ds << " multithread index_kcore construction!" << endl;
	auto start = chrono::system_clock::now();
	multithread_index_construction(g, threads);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << elapsed_seconds.count() << endl;
	result << "1.0" << "," << elapsed_seconds.count() << endl;
	result.close();
}

void index_construction_scala_vary_edge_single_test(Datasets ds, int threads) {
	stringstream output_location;
	output_location << "../experiments/multithread_scala_edge_" << ds << ".txt";
	ofstream result; result.open(output_location.str());
	// node test
	result << "sampling_rate" << "," << "time" << endl;
	for (int p = 0; p < 4; p++) {
		stringstream graph_location; graph_location << "../scalability/" << dataset2name[ds] << "edge" << (p * 20 + 20);
		BiGraph g(graph_location.str());
		cout << "start " << ds << " edge scale " << (p * 20 + 20) << " multithread index_kcore construction!" << endl;
		auto start = chrono::system_clock::now();
		multithread_index_construction(g, threads);
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		cout << elapsed_seconds.count() << endl;
		result << double(p * 20 + 20) / 100 << "," << elapsed_seconds.count() << endl;
	}
	stringstream graph_location; graph_location << "../realGraph/" << dataset2name[ds];
	BiGraph g(graph_location.str());
	cout << "start " << ds << " multithread index_kcore construction!" << endl;
	auto start = chrono::system_clock::now();
	multithread_index_construction(g, threads);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << elapsed_seconds.count() << endl;
	result << "1.0" << "," << elapsed_seconds.count() << endl;
	result.close();
}

void index_construction_scalability() {
	index_construction_scala_vary_node_single_test(Datasets::TR, 12);
	index_construction_scala_vary_node_single_test(Datasets::WT, 12);
	index_construction_scala_vary_node_single_test(Datasets::OG, 12);

	index_construction_scala_vary_edge_single_test(Datasets::TR, 12);
	index_construction_scala_vary_edge_single_test(Datasets::WT, 12);
	index_construction_scala_vary_edge_single_test(Datasets::OG, 12);
}

void parallel_vary_core_insertion_single_test(Datasets ds, Index_update iu) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/parallel_varycore_insertion_" << ds << "_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	cout << "begin parallel varycore insertion test " << ds << " " << iu << endl;
	double time = 0;
	for (int i = 2; i <= 20; i += 2) {
		//dynamic_test_removal_single_test(Datasets::OG, Index_update::withlimit_parallel, 12);
		if (iu == Index_update::withlimit_dfs_parallel) {
			time = dynamic_test_insertion_single_test(ds, iu, i);
		}
		else if (iu == Index_update::withlimit_parallel) {
			time = dynamic_test_insertion_single_test(ds, iu, i);
		}
		cout << i << " " << time << endl;
		result << i << "	" << time << endl;
	}
	result.close();
}

void parallel_vary_core_removal_single_test(Datasets ds, Index_update iu) {
	stringstream output_location;
	output_location << "../experiments/dynamic_exp/parallel_varycore_removal_" << ds << "_" << iu << ".txt";
	ofstream result; result.open(output_location.str());
	cout << "begin parallel varycore removal test " << ds << " " << iu << endl;
	double time = 0;
	for (int i = 2; i <= 20; i += 2) {
		//dynamic_test_removal_single_test(Datasets::OG, Index_update::withlimit_parallel, 12);
		if (iu == Index_update::withlimit_dfs_parallel) {
			time = dynamic_test_removal_single_test(ds, iu, i);
		}
		else if (iu == Index_update::withlimit_parallel) {
			time = dynamic_test_removal_single_test(ds, iu, i);
		}
		cout << i << " " << time << endl;
		result << i << "	" << time << endl;
	}
	result.close();
}

void parallel_vary_core_removal() {
	DYNAMIC_TEST_NUM = 10;
	parallel_vary_core_removal_single_test(Datasets::TR, Index_update::withlimit_dfs_parallel);

	DYNAMIC_TEST_NUM = 10;
	parallel_vary_core_removal_single_test(Datasets::WT, Index_update::withlimit_dfs_parallel);

	DYNAMIC_TEST_NUM = 10;
	parallel_vary_core_removal_single_test(Datasets::OG, Index_update::withlimit_dfs_parallel);
}

void parallel_vary_core_insertion() {
	DYNAMIC_TEST_NUM = 10;
	parallel_vary_core_insertion_single_test(Datasets::TR, Index_update::withlimit_dfs_parallel);

	DYNAMIC_TEST_NUM = 10;
	parallel_vary_core_insertion_single_test(Datasets::WT, Index_update::withlimit_dfs_parallel);

	DYNAMIC_TEST_NUM = 10;
	parallel_vary_core_insertion_single_test(Datasets::OG, Index_update::withlimit_dfs_parallel);
}


void ext_exp(string graph_location, int threads, string output_location) {
	ofstream result_file(output_location);
	BiGraph g(graph_location);
	cout << "graph loaded" << endl;
	BiGraph tmp = g;
	auto start = chrono::system_clock::now();
	int delta = coreIndexKCore(g);
	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	result_file << "delta: " << delta << endl;
	cout << "nonmultithread running time: " << elapsed_seconds.count() << endl;
	result_file << "nonmultithread running time: " << elapsed_seconds.count() << endl;
	start = chrono::system_clock::now();
	multithread_index_construction(tmp, threads);
	end = chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "multithread running time: " << elapsed_seconds.count() << endl;
	result_file << "multithread running time: " << elapsed_seconds.count() << endl;
	if (isIndexEqual(g, tmp)) {
		cout << "correct" << endl;
		result_file << "correct" << endl;
	}
	else {
		cout << "error" << endl;
		result_file << "error" << endl;
	}
	vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
	cout << "start bicore index construction!" << endl;
	build_bicore_index(g, bicore_index_u, bicore_index_v);
	time_t _start, _end, _time;
	vector<int> standard = { 400 };
	for (int j = 0; j < standard.size(); j++) {
		int alpha = standard[j];
		result_file << "--------------fix alpha = " << alpha << "-------------" << endl;
		result_file << "beta" << "," << "peeling" << "," << "naiveindex" << "," << "bicoreIndex" << endl;
		for (int beta = 40; beta <= 400; beta += 40) {
			result_file << beta << ",";
			vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2; vector<bool> left3; vector<bool> right3;
			left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2); left3.resize(g.num_v1); right3.resize(g.num_v2);
			_start = clock();
			for (int ijk = 0; ijk < 100; ijk++) {
				findabcore_peeling(alpha, beta, g, left1, right1);
			}
			_end = clock();
			_time = _end - _start;
			result_file << _time << ",";
			_start = clock();
			for (int ijk = 0; ijk < 100; ijk++) {
				findabcore_core_index(alpha, beta, g, left2, right2);
			}
			_end = clock();
			_time = _end - _start;
			result_file << _time << ",";
			_start = clock();
			for (int ijk = 0; ijk < 100; ijk++) {
				retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left3, right3, alpha, beta);
			}
			_end = clock();
			_time = _end - _start;
			result_file << _time << endl;
			if (!test_abcore_equal(left1, right1, left3, right3)) {
				cout << "error" << endl;
				result_file << "error" << endl;
			}
		}
	}
	for (int j = 0; j < standard.size(); j++) {
		int beta = standard[j];
		result_file << "--------------fix beta = " << beta << "-------------" << endl;
		result_file << "alpha" << "," << "peeling" << "," << "naiveindex" << "," << "bicoreIndex" << endl;
		for (int alpha = 40; alpha <= 400; alpha += 40) {
			result_file << alpha << ",";
			vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2; vector<bool> left3; vector<bool> right3;
			left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2); left3.resize(g.num_v1); right3.resize(g.num_v2);
			_start = clock();
			for (int ijk = 0; ijk < 100; ijk++) {
				findabcore_peeling(alpha, beta, g, left1, right1);
			}
			_end = clock();
			_time = _end - _start;
			result_file << _time << ",";
			_start = clock();
			for (int ijk = 0; ijk < 100; ijk++) {
				findabcore_core_index(alpha, beta, g, left2, right2);
			}
			_end = clock();
			_time = _end - _start;
			result_file << _time << ",";
			_start = clock();
			for (int ijk = 0; ijk < 100; ijk++) {
				retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left3, right3, alpha, beta);
			}
			_end = clock();
			_time = _end - _start;
			result_file << _time << endl;
			if (!test_abcore_equal(left1, right1, left3, right3)) {
				cout << "error" << endl;
				result_file << "error" << endl;
			}
		}
	}
	result_file << "--------------fix beta = 10 and alpha = 10-------------" << endl;
	result_file << "peeling" << "," << "naiveindex" << "," << "bicoreIndex" << endl;
	vector<bool> left1; vector<bool> right1; vector<bool> left2; vector<bool> right2; vector<bool> left3; vector<bool> right3;
	left1.resize(g.num_v1); right1.resize(g.num_v2); left2.resize(g.num_v1); right2.resize(g.num_v2); left3.resize(g.num_v1); right3.resize(g.num_v2);
	_start = clock();
	for (int ijk = 0; ijk < 100; ijk++) {
		findabcore_peeling(10, 10, g, left1, right1);
	}
	_end = clock();
	_time = _end - _start;
	result_file << _time << ",";
	_start = clock();
	for (int ijk = 0; ijk < 100; ijk++) {
		findabcore_core_index(10, 10, g, left2, right2);
	}
	_end = clock();
	_time = _end - _start;
	result_file << _time << ",";
	_start = clock();
	for (int ijk = 0; ijk < 100; ijk++) {
		retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left3, right3, 10, 10);
	}
	_end = clock();
	_time = _end - _start;
	result_file << _time << endl;
	if (!test_abcore_equal(left1, right1, left3, right3)) {
		cout << "error" << endl;
		result_file << "error" << endl;
	}
	result_file.close();
}