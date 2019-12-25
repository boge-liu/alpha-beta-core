#include "random.h"
#include "static.h"
#include "dynamic.h"
#include "gephi.h"
#include "preprocessing.h"
#include "incremental_test.h"
#include <ctime>
#include "experiment.h"
#include "paper.h"
#include "mischievous.h"
#include "Multithread.h"
using namespace std;

vector<string> DATASET = { "wiki-en-cat_", "flickr-groupmemberships_",
"epinions-rating_", "edit-dewiki_", "reuters_", "gottron-trec_",
"delicious-ui_", "livejournal-groupmemberships_", "trackers_", "orkut-groupmemberships_"
};
vector<string> ALIAS = { "WC", "FG",
"EP", "DE", "RE", "TR",
"DUI", "LJ", "WT", "OG"
};



int main(int argc, char **argv) {
	if (argc == 1) {
		cout << "error in number of arguments" << endl;
	}
	string exec_type = argv[1];
	if (exec_type == "-BasicDecom") {
		cout << "start BasicDecom for " << argv[2] << endl;
		BiGraph g(argv[2]);
		auto start = chrono::system_clock::now();
		coreIndexBasic(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "run time: " << time.count() << endl;
	}
	else if (exec_type == "-ComShrDecom") {
		cout << "start ComShrDecom for " << argv[2] << endl;
		BiGraph g(argv[2]);
		auto start = chrono::system_clock::now();
		coreIndexKCore(g);
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "run time: " << time.count() << endl;
	}
	else if (exec_type == "-ParallelDecom") {
		cout << "start ParallelDecom for " << argv[2] << endl;
		BiGraph g(argv[2]);
		multithread_index_construction(g, stoi(argv[3]));
	}
	else if (exec_type == "-Query") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		vector<bool> left; vector<bool> right;
		// all the vertices in query result are set as true
		left.resize(g.num_v1, false); right.resize(g.num_v2, false);
		auto start = chrono::system_clock::now();
		retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left, right, stoi(argv[3]), stoi(argv[4]));
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "query time: " << time.count() << endl;
	}
	else if (exec_type == "-BiCore-Index-Ins") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), true);
		cout << "BiCore-Index-Ins running time: " << time << endl;
	}
	else if (exec_type == "-BiCore-Index-Rem") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_with_limit_swap(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), false);
		cout << "BiCore-Index-Rem running time: " << time << endl;
	}
	else if (exec_type == "-BiCore-Index-Ins*") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), true);
		cout << "BiCore-Index-Ins* running time: " << time << endl;
	}
	else if (exec_type == "-BiCore-Index-Rem*") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_swap_with_dfs(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), false);
		cout << "BiCore-Index-Rem* running time: " << time << endl;
	}
	else if (exec_type == "-ParallelIns") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), true, stoi(argv[5]));
		cout << "ParallelIns running time: " << time << endl;
	}
	else if (exec_type == "-ParallelRem") {
		BiGraph g(argv[2]);
		coreIndexKCore(g);
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
		build_bicore_index(g, bicore_index_u, bicore_index_v);
		auto time = Dyn_rebuild::update_bicore_index_with_limit_swap_dfs_parallel(g, bicore_index_u, bicore_index_v, stoi(argv[3]), stoi(argv[4]), false, stoi(argv[5]));
		cout << "ParallelRem running time: " << time << endl;
	}
	else {
		cout << "illegal arguments" << endl;
	}
	return 0;
}