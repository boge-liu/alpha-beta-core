#include"static.h"
#include <algorithm>
#include <ctime>
#include <iostream>
#include <set>
#include <list>
using namespace std;

void binsort(vector<list<vid_t>>& left_bin, vector<list<vid_t>>& right_bin,
	vector<list<vid_t>::iterator>& left_list_pointer, vector<list<vid_t>::iterator>& right_list_pointer, BiGraph& tmp) {
	left_bin.clear(); right_bin.clear(); left_list_pointer.clear(); right_list_pointer.clear();
	left_bin.resize(tmp.v1_max_degree+1);
	left_list_pointer.resize(tmp.getV1Num());
	for (int i = 0; i < tmp.getV1Num(); i++) {
		left_bin[tmp.getV1Degree(i)].push_back(i);
		left_list_pointer[i] = left_bin[tmp.getV1Degree(i)].end();
		left_list_pointer[i]--;
	}
	right_list_pointer.resize(tmp.getV2Num());
	right_bin.resize(tmp.v2_max_degree + 1);
	for (int i = 0; i < tmp.getV2Num(); i++) {
		right_bin[tmp.getV2Degree(i)].push_back(i);
		right_list_pointer[i] = right_bin[tmp.getV2Degree(i)].end();
		right_list_pointer[i]--;
	}
}

void reorder_left(vector<list<vid_t>>& bin, vector<list<vid_t>::iterator>& list_pointer, int new_degree, int old_degree, vid_t vertex) {
	bin[old_degree].erase(list_pointer[vertex]);
	bin[new_degree].push_back(vertex);
	list_pointer[vertex] = bin[new_degree].end();
	list_pointer[vertex]--;
}

void reorder_right(vector<list<vid_t>>& bin, vector<list<vid_t>::iterator>& list_pointer, int new_degree, int old_degree, vid_t vertex) {
	bin[old_degree].erase(list_pointer[vertex]);
	bin[new_degree].push_back(vertex);
	list_pointer[vertex] = bin[new_degree].end();
	list_pointer[vertex]--;
}

void update_x_k_with_k_x(BiGraph& g, int left_k, int k_x, vid_t v) {
	for (int right_k = 1; right_k < g.right_index[v].size(); right_k++) {
		if (right_k <= k_x) {
			if (g.right_index[v][right_k] < left_k) {
				g.right_index[v][right_k] = left_k;
			}
		}
		else {
			break;
		}
	}
}

void peel_with_fixed_left_k(int left_k, BiGraph& tmp, BiGraph& g) {
	vector<list<vid_t>> left_bin;
	vector<list<vid_t>::iterator> left_list_pointer;
	vector<list<vid_t>> right_bin;
	vector<list<vid_t>::iterator> right_list_pointer;
	binsort(left_bin, right_bin, left_list_pointer, right_list_pointer, tmp);
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		bool left_unchange = false;
		bool right_unchange = false;
		while (!left_unchange || !right_unchange) {
			left_unchange = true;
			right_unchange = true;
			// peel left
			for (int i = 1; i < left_k; i++) {
				for (auto j = left_bin[i].begin(); j != left_bin[i].end();) {
					vid_t u = *j;
					vector<vid_t> neigh = tmp.getV1Neighbors(u);
					int old_degree = tmp.getV1Degree(u); // should equal to i
					for (int k = 0; k < old_degree; k++) {
						vid_t v = neigh[k];
						int old_deg = tmp.getV2Degree(v); // should equal to tmp.getV2Degree(v) + 1
						tmp.deleteEdge(u, v);
						reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), old_deg, v);
						// degree = 0 will not be considered in next round
						if (tmp.getV2Degree(v) == 0 && right_k - 1 > 0) {
							update_x_k_with_k_x(g, left_k, right_k - 1, v);
						}
					}
					++j; // j becomes unavailable after reordering
					reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), old_degree, u);
					if (right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
					}
					left_unchange = false;
				}
			}
			// peel right
			for (int i = 1; i < right_k; i++) {
				for (auto j = right_bin[i].begin(); j != right_bin[i].end(); ) {
					vid_t v = *j;
					vector<vid_t> neigh = tmp.getV2Neighbors(v);
					int old_degree = tmp.getV2Degree(v);
					for (int k = 0; k < old_degree; k++) {
						vid_t u = neigh[k];
						int old_deg = tmp.getV1Degree(u);
						tmp.deleteEdge(u, v);
						reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), old_deg, u);
						// degree = 0 will not be considered in next round
						if (tmp.getV1Degree(u) == 0 && right_k - 1 > 0) {
							g.left_index[u][left_k] = right_k - 1;
						}
					}
					++j;
					reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), old_degree, v);
					if (right_k - 1 > 0) {
						update_x_k_with_k_x(g, left_k, right_k - 1, v);
					}					
					right_unchange = false;
				}
			}
		}
	}
}

void inverse(BiGraph& g) {
	swap(g.degree_v1,g.degree_v2);
	swap(g.left_index, g.right_index);
	swap(g.num_v1, g.num_v2);
	swap(g.neighbor_v1, g.neighbor_v2);
	swap(g.neighborHash_v1, g.neighborHash_v2);
	swap(g.v1_max_degree, g.v2_max_degree);
	swap(g.left_delete, g.right_delete);
}

void static_construct_index_with_peeling_algorithm_basic(BiGraph& g) {
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
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i)+1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i)+1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	for (int i = 1; i <= g.v1_max_degree; i++) {
		BiGraph tmp = g;
		peel_with_fixed_left_k(i, tmp, g);
	}
	// calculate right
	/*inverse(g);
	for (int i = 1; i <= g.v1_max_degree; i++) {
		BiGraph tmp = g;
		peel_with_fixed_left_k(i, tmp, g);
	}
	inverse(g);*/
}

// with left fixed
void calculate_left_side_upper_bound(BiGraph& g, vid_t u, int left_k, vector<int>& right_tmp_k){
	vector<vid_t> neigh = g.getV1Neighbors(u);
	int degree = g.getV1Degree(u);	
	// bin sort
	vector<int> bins;
	bins.resize(g.v2_max_degree + 1);
	fill_n(bins.begin(), bins.size(), 0);
	for (int i = 0; i < degree; i++) {
		bins[right_tmp_k[neigh[i]]]++;
	}
	for (int i = bins.size() - 1; i > 0; i--) {
		if (bins[i] >= left_k) {
			g.left_index[u][left_k] = i;
			return;
		}
		bins[i - 1] += bins[i];		
	}
}

// with left fixed
void calculate_right_side_upper_bound(BiGraph& g, vid_t v, int left_k, vector<int>& right_tmp_k) {
	vector<vid_t> neigh = g.getV2Neighbors(v);
	int degree = g.getV2Degree(v);
	// bin sort
	vector<int> bins;
	bins.resize(g.v2_max_degree + 1);
	fill_n(bins.begin(), bins.size(), 0);
	for (int i = 0; i < degree; i++) {
		if (g.getV1Degree(neigh[i]) < left_k) continue;
		bins[g.left_index[neigh[i]][left_k]]++;
	}
	for (int i = bins.size() - 1; i > 0; i--) {
		if (bins[i] >= i) {
			right_tmp_k[v] = i;
			return;
		}
		bins[i - 1] += bins[i];
	}
}

// zhu yi left_k yue jie de wen ti
void converge_with_fixed_left_k(int left_k, BiGraph& g) {
	vector<int> right_tmp_k;
	right_tmp_k.resize(g.getV2Num());
	for (int i = 0; i < right_tmp_k.size(); i++) {
		right_tmp_k[i] = g.getV2Degree(i);
	}
	while (true) {
		// left side
		bool unchange = true;		
		for (vid_t i = 0; i < g.getV1Num(); i++) {
			if (g.getV1Degree(i) < left_k) continue;
			int old_value = g.left_index[i][left_k];
			calculate_left_side_upper_bound(g, i, left_k, right_tmp_k);
			int new_value = g.left_index[i][left_k];
			if (old_value != new_value) unchange = false;
		}
		if (unchange) break;
		// right_side
		unchange = true;
		for (vid_t i = 0; i < g.getV2Num(); i++) {
			//if (g.getV1Degree(i) < left_k) continue;
			int old_value = right_tmp_k[i];
			calculate_right_side_upper_bound(g, i, left_k, right_tmp_k);
			int new_value = right_tmp_k[i];
			if (old_value != new_value) unchange = false;
		}
		if (unchange) break;
	}
}

void static_construct_index_with_convergence_algorithm_basic(BiGraph& g) {
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
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	for (int i = 1; i <= g.v1_max_degree; i++) {
		converge_with_fixed_left_k(i, g);
	}
	inverse(g);
	for (int i = 1; i <= g.v1_max_degree; i++) {
		converge_with_fixed_left_k(i, g);
	}
	inverse(g);
}

void peel_with_fixed_left_k_acc(int left_k, BiGraph& tmp, BiGraph& g) {
	BiGraph newtmp;
	vector<list<vid_t>> left_bin;
	vector<list<vid_t>::iterator> left_list_pointer;
	vector<list<vid_t>> right_bin;
	vector<list<vid_t>::iterator> right_list_pointer;
	binsort(left_bin, right_bin, left_list_pointer, right_list_pointer, tmp);
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		bool left_unchange = false;
		bool right_unchange = false;
		while (!left_unchange || !right_unchange) {
			left_unchange = true;
			right_unchange = true;
			// peel left
			for (int i = 1; i < left_k; i++) {
				for (auto j = left_bin[i].begin(); j != left_bin[i].end();) {
					vid_t u = *j;
					vector<vid_t> neigh = tmp.getV1Neighbors(u);
					int old_degree = tmp.getV1Degree(u); // should equal to i
					for (int k = 0; k < old_degree; k++) {
						vid_t v = neigh[k];
						int old_deg = tmp.getV2Degree(v); // should equal to tmp.getV2Degree(v) + 1
						tmp.deleteEdge(u, v);
						reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), old_deg, v);
						// degree = 0 will not be considered in next round
						if (tmp.getV2Degree(v) == 0 && right_k - 1 > 0) {
							update_x_k_with_k_x(g, left_k, right_k - 1, v);
						}
					}
					++j; // j becomes unavailable after reordering
					reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), old_degree, u);
					if (right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
					}
					left_unchange = false;
				}
			}
			// peel right
			for (int i = 1; i < right_k; i++) {
				for (auto j = right_bin[i].begin(); j != right_bin[i].end(); ) {
					vid_t v = *j;
					vector<vid_t> neigh = tmp.getV2Neighbors(v);
					int old_degree = tmp.getV2Degree(v);
					for (int k = 0; k < old_degree; k++) {
						vid_t u = neigh[k];
						int old_deg = tmp.getV1Degree(u);
						tmp.deleteEdge(u, v);
						reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), old_deg, u);
						// degree = 0 will not be considered in next round
						if (tmp.getV1Degree(u) == 0 && right_k - 1 > 0) {
							g.left_index[u][left_k] = right_k - 1;
						}
					}
					++j;
					reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), old_degree, v);
					if (right_k - 1 > 0) {
						update_x_k_with_k_x(g, left_k, right_k - 1, v);
					}
					right_unchange = false;
				}
			}
		}
		if (right_k == 1) {
			newtmp = tmp;
		}
	}
	tmp = newtmp;
}

void static_construct_index_with_peeling_algorithm_acc(BiGraph& g) {
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
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	BiGraph tmp = g;
	for (int i = 1; i <= g.v1_max_degree; i++) {		
		peel_with_fixed_left_k_acc(i, tmp, g);
	}
	// calculate right
	/*inverse(g);
	tmp = g;
	for (int i = 1; i <= g.v1_max_degree; i++) {		
		peel_with_fixed_left_k_acc(i, tmp, g);
	}
	inverse(g);*/
}

// accelerate convergence
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////
vector<int> left_binns;
vector<int> right_binns;
void compute_left_side_k_x_upper_bound_and_score_acc(BiGraph& g, vid_t u, int left_k, vector<int>& left_tmp_k_x, vector<int>& left_score, const vector<int>& right_tmp_k_x) {
	//vector<vid_t> neigh = g.getV1Neighbors(u);
	int degree = g.getV1Degree(u);
	//left_binns.resize(left_tmp_k_x[u] + 1);
	int boundary = left_tmp_k_x[u] + 1;
	fill_n(left_binns.begin(), boundary, 0);
	for (int i = 0; i < degree; i++) {
		int right_k = right_tmp_k_x[g.neighbor_v1[u][i]];
		if (right_k >= boundary) left_binns[boundary - 1]++;
		else left_binns[right_k]++;
	}
	for (int right_k = boundary - 1; right_k > 0; right_k--) {
		if (left_binns[right_k] >= left_k) {
			left_score[u] = left_binns[right_k];
			left_tmp_k_x[u] = right_k;
			return;
		}
		left_binns[right_k - 1] += left_binns[right_k];
	}
	// retire
	left_tmp_k_x[u] = 0;
	left_score[u] = g.v1_max_degree + g.v2_max_degree + 100;
}

void compute_right_side_k_x_upper_bound_and_score_acc(BiGraph& g, vid_t v, int left_k, vector<int>& right_tmp_k_x, vector<int>& right_score, const vector<int>& left_tmp_k_x) {
	//vector<vid_t> neigh = g.getV2Neighbors(v);
	int degree = g.getV2Degree(v);
	//bins.resize(degree + 1);
	int boundary = degree + 1;
	fill_n(right_binns.begin(), boundary, 0);
	for (int i = 0; i < degree; i++) {
		int right_k = left_tmp_k_x[g.neighbor_v2[v][i]];
		if (right_k >= boundary) right_binns[boundary - 1]++;
		else right_binns[right_k]++;
	}
	for (int right_k = boundary - 1; right_k > 0; right_k--) {
		if (right_binns[right_k] >= right_k) {
			right_score[v] = right_binns[right_k];
			right_tmp_k_x[v] = right_k;
			return;
		}
		right_binns[right_k - 1] += right_binns[right_k];
	}
	// retire
	right_tmp_k_x[v] = 0;
	right_score[v] = g.v1_max_degree + g.v2_max_degree + 100;
}

void converge_with_fixed_left_k_acc(int left_k, BiGraph& g, vector<int>& left_tmp_k_x, vector<int>& right_tmp_k_x, vector<int>& left_score, vector<int>& right_score) {
	// first time initialize
	if (left_k == 1) {
		right_score.clear(); right_score.resize(g.getV2Num());
		right_tmp_k_x.clear(); right_tmp_k_x.resize(g.getV2Num());
		for (vid_t v = 0; v < right_tmp_k_x.size(); v++) {
			right_tmp_k_x[v] = g.getV2Degree(v);
			right_score[v] = g.getV2Degree(v);
		}
		left_tmp_k_x.clear(); left_tmp_k_x.resize(g.getV1Num()); fill_n(left_tmp_k_x.begin(), left_tmp_k_x.size(), -1);
		left_score.clear(); left_score.resize(g.getV1Num()); fill_n(left_score.begin(), left_score.size(), -1);
		for (vid_t u = 0; u < g.getV1Num(); u++) {
			//vector<vid_t> neigh = g.getV1Neighbors(u);
			int degree = g.getV1Degree(u);
			int k_x = 0; int score = 0;
			for (int i = 0; i < degree; i++) {
				vid_t v = g.neighbor_v1[u][i];
				int deg = g.getV2Degree(v);
				if (deg == k_x) {
					score++;
				}
				else if (deg > k_x) {
					k_x = deg;
					score = 1;
				}
			}
			left_tmp_k_x[u] = k_x;
			left_score[u] = score;
		}
	}
	// initialize recompute nodes
	vector<vid_t> left_recompute;
	vector<vid_t> right_recompute;
	for (vid_t u = 0; u < g.getV1Num(); u++) {
		if (left_score[u] >= left_k) continue;
		else left_recompute.push_back(u);
	}
	/*for (vid_t v = 0; v < g.getV2Num(); v++) {
		if (right_score[v] >= right_tmp_k_x[v]) continue;
		else right_recompute.insert(v);
	}*/
	while (!left_recompute.empty() || !right_recompute.empty()) {
		// left side
		for (auto it = left_recompute.begin(); it != left_recompute.end(); it++) {
			vid_t u = *it;
			int old_value = left_tmp_k_x[u];
			compute_left_side_k_x_upper_bound_and_score_acc(g, u, left_k, left_tmp_k_x, left_score, right_tmp_k_x);
			int new_value = left_tmp_k_x[u];
			//vector<vid_t> neigh = g.getV1Neighbors(u);
			int degree = g.getV1Degree(u);
			for (int i = 0; i < degree; i++) {
				vid_t v = g.neighbor_v1[u][i];
				if (right_tmp_k_x[v] <= old_value && right_tmp_k_x[v] > new_value) {
					right_score[v]--;				
				}
			}
		}
		for (vid_t v = 0; v < g.getV2Num(); v++) {
			if (right_score[v] < right_tmp_k_x[v]) {
				right_recompute.push_back(v);
			}
		}
		left_recompute.clear();
		// right side 
		for (auto it = right_recompute.begin(); it != right_recompute.end(); it++) {
			vid_t v = *it;
			int old_value = right_tmp_k_x[v];			
			compute_right_side_k_x_upper_bound_and_score_acc(g, v, left_k, right_tmp_k_x, right_score, left_tmp_k_x);
			int new_value = right_tmp_k_x[v];
			//vector<vid_t> neigh = g.getV2Neighbors(v);
			int degree = g.getV2Degree(v);
			for (int i = 0; i < degree; i++) {
				vid_t u = g.neighbor_v2[v][i];
				if (left_tmp_k_x[u] <= old_value && left_tmp_k_x[u] > new_value) {
					left_score[u]--;				
				}
			}
		}
		for (vid_t u = 0; u < g.getV1Num(); u++) {
			if (left_score[u] < left_k) {
				left_recompute.push_back(u);
			}
		}
		right_recompute.clear();
	}
	// update k_x and x_k
	for (vid_t u = 0; u < left_tmp_k_x.size(); u++) {
		if (g.getV1Degree(u) >= left_k) {
			g.left_index[u][left_k] = left_tmp_k_x[u];
		}
	}
	for (vid_t v = 0; v < right_tmp_k_x.size(); v++) {
		update_x_k_with_k_x(g, left_k, right_tmp_k_x[v], v);
	}
}

void static_construct_index_with_convergence_algorithm_acc(BiGraph& g) {
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
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	vector<int> right_tmp_k_x;
	vector<int> left_tmp_k_x;
	vector<int> right_score;
	vector<int> left_score;
	left_binns.resize(g.v2_max_degree + 1);
	right_binns.resize(g.v2_max_degree + 1);
	for (int i = 1; i <= g.v1_max_degree; i++) {
		converge_with_fixed_left_k_acc(i, g, left_tmp_k_x, right_tmp_k_x, left_score, right_score);
	}
}

// long yuan's peel
//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////

void update_k_x_with_x_k(BiGraph& g, int right_k, int x_k, vid_t u) {
	for (int left_k = 1; left_k < g.left_index[u].size(); left_k++) {
		if (left_k <= x_k) {
			if (g.left_index[u][left_k] < right_k) {
				g.left_index[u][left_k] = right_k;
			}
		}
		else {
			break;
		}
	}
}

void peel_with_fixed_right_k_acc(int right_k, BiGraph& tmp, BiGraph& g) {
	BiGraph newtmp;
	vector<list<vid_t>> left_bin;
	vector<list<vid_t>::iterator> left_list_pointer;
	vector<list<vid_t>> right_bin;
	vector<list<vid_t>::iterator> right_list_pointer;
	binsort(left_bin, right_bin, left_list_pointer, right_list_pointer, tmp);
	for (int left_k = 1; left_k <= g.v1_max_degree + 1; left_k++) {
		bool left_unchange = false;
		bool right_unchange = false;
		while (!left_unchange || !right_unchange) {
			left_unchange = true;
			right_unchange = true;
			// peel right
			for (int i = 1; i < right_k; i++) {
				for (auto j = right_bin[i].begin(); j != right_bin[i].end();) {
					vid_t v = *j;
					vector<vid_t> neigh = tmp.getV2Neighbors(v);
					int old_degree = tmp.getV2Degree(v); // should equal to i
					for (int k = 0; k < old_degree; k++) {
						vid_t u = neigh[k];
						int old_deg = tmp.getV1Degree(u); // should equal to tmp.getV2Degree(v) + 1
						tmp.deleteEdge(u, v);
						reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), old_deg, u);
						// degree = 0 will not be considered in next round
						if (tmp.getV1Degree(u) == 0 && left_k - 1 > 0) {
							update_k_x_with_x_k(g, right_k, left_k - 1, u);
						}
					}
					++j; // j becomes unavailable after reordering
					reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), old_degree, v);
					if (left_k - 1 > 0) {
						g.right_index[v][right_k] = left_k - 1;
					}
					right_unchange = false;
				}
			}
			// peel left
			for (int i = 1; i < left_k; i++) {
				for (auto j = left_bin[i].begin(); j != left_bin[i].end(); ) {
					vid_t u = *j;
					vector<vid_t> neigh = tmp.getV1Neighbors(u);
					int old_degree = tmp.getV1Degree(u);
					for (int k = 0; k < old_degree; k++) {
						vid_t v = neigh[k];
						int old_deg = tmp.getV2Degree(v);
						tmp.deleteEdge(u, v);
						reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), old_deg, v);
						// degree = 0 will not be considered in next round
						if (tmp.getV2Degree(v) == 0 && left_k - 1 > 0) {
							g.right_index[v][right_k] = left_k - 1;
						}
					}
					++j;
					reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), old_degree, u);
					if (left_k - 1 > 0) {
						update_k_x_with_x_k(g, right_k, left_k - 1, u);
					}
					left_unchange = false;
				}
			}
		}
		if (left_k == 1) {
			newtmp = tmp;
		}
	}
	tmp = newtmp;
}

void static_construct_index_with_peeling_algorithm_long_yuan(BiGraph& g) {
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
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	BiGraph tmp = g;
	int right_max = 0;
	for (int left_k = 1; left_k <= g.v1_max_degree; left_k++) {
		peel_with_fixed_left_k_acc(left_k, tmp, g);
		right_max = 0;
		for (vid_t u = 0; u < g.getV1Num(); u++) {
			if (g.left_index[u].size() <= left_k) continue;
			int right_k = g.left_index[u][left_k];
			if (right_max < right_k) right_max = right_k;
		}
		if (right_max <= left_k) break;
	}
	
	inverse(g);
	tmp = g;
	for (int left_k = 1; left_k <= right_max; left_k++) {
		peel_with_fixed_left_k_acc(left_k, tmp, g);
	}
	inverse(g);
}


// no edge deletion
//////////////////////////////////////
//////////////////////////////////////
//////////////////////////////////////

void peel_with_fixed_left_k_acc_no_edge_deletion(int left_k, BiGraph& tmp, BiGraph& g) {
	BiGraph newtmp;
	vector<list<vid_t>> left_bin;
	vector<list<vid_t>::iterator> left_list_pointer;
	vector<list<vid_t>> right_bin;
	vector<list<vid_t>::iterator> right_list_pointer;
	binsort(left_bin, right_bin, left_list_pointer, right_list_pointer, tmp);
	// avoid too many reorder times
	vector<int> right_side_old_degree;
	for (vid_t v = 0; v < tmp.getV2Num(); v++) {
		right_side_old_degree.push_back(tmp.getV2Degree(v));
	}
	vector<int> left_side_old_degree;
	for (vid_t u = 0; u < tmp.getV1Num(); u++) {
		left_side_old_degree.push_back(tmp.getV1Degree(u));
	}
	for (int right_k = 1; right_k <= tmp.v2_max_degree + 1; right_k++) {
		bool left_unchange = false;
		bool right_unchange = false;
		while (!left_unchange || !right_unchange) {
			left_unchange = true;
			right_unchange = true;
			// peel left
			vector<vid_t> reorder_right_vertices;
			for (int i = 1; i < left_k; i++) {
				for (auto j = left_bin[i].begin(); j != left_bin[i].end();) {
					vid_t u = *j;
					//vector<vid_t> neigh = tmp.getV1Neighbors(u);
					int old_degree = tmp.getV1Degree(u); // should equal to i
					for (int k = 0; k < tmp.neighbor_v1[u].size(); k++) {
						vid_t v = tmp.neighbor_v1[u][k];
						if (tmp.right_delete[v]) continue;
						int old_deg = tmp.getV2Degree(v); // should equal to tmp.getV2Degree(v) + 1
						//tmp.deleteEdge(u, v);
						tmp.degree_v2[v]--;
						tmp.degree_v1[u]--;
						reorder_right_vertices.push_back(v);
						//reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), old_deg, v);
						// degree = 0 will not be considered in next round
						if (tmp.getV2Degree(v) == 0 && right_k - 1 > 0) {
							update_x_k_with_k_x(g, left_k, right_k - 1, v);
						}
					}
					tmp.left_delete[u] = true;
					if (tmp.getV1Degree(u) != 0) {
						int kkkkkk = 0;
					}
					++j; // j becomes unavailable after reordering
					reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), old_degree, u);
					if (right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
					}
					left_unchange = false;
				}
			}
			for (int i = 0; i < reorder_right_vertices.size(); i++) {
				vid_t v = reorder_right_vertices[i];
				if (tmp.right_delete[v]) continue;
				if (right_side_old_degree[v] != tmp.getV2Degree(v)) {
					right_side_old_degree[v] = tmp.getV2Degree(v);
					reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), right_side_old_degree[v], v);
					if (tmp.getV2Degree(v) == 0) {
						tmp.right_delete[v] = true;
					}
				}
			}
			// peel right
			vector<vid_t> reorder_left_vertices;
			for (int i = 1; i < right_k; i++) {
				for (auto j = right_bin[i].begin(); j != right_bin[i].end(); ) {
					vid_t v = *j;
					//vector<vid_t> neigh = tmp.getV2Neighbors(v);
					int old_degree = tmp.getV2Degree(v);
					for (int k = 0; k < tmp.neighbor_v2[v].size(); k++) {
						vid_t u = tmp.neighbor_v2[v][k];
						if (tmp.left_delete[u]) continue;
						int old_deg = tmp.getV1Degree(u);
						//tmp.deleteEdge(u, v);
						tmp.degree_v2[v]--;
						tmp.degree_v1[u]--;
						reorder_left_vertices.push_back(u);
						//reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), old_deg, u);
						// degree = 0 will not be considered in next round
						if (tmp.getV1Degree(u) == 0 && right_k - 1 > 0) {
							g.left_index[u][left_k] = right_k - 1;
						}
					}
					tmp.right_delete[v] = true;
					if (tmp.getV2Degree(v) != 0) {
						int kkkkkk = 0;
					}
					++j;
					reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), old_degree, v);
					if (right_k - 1 > 0) {
						update_x_k_with_k_x(g, left_k, right_k - 1, v);
					}
					right_unchange = false;
				}
			}
			for (int i = 0; i < reorder_left_vertices.size(); i++) {
				vid_t u = reorder_left_vertices[i];
				if (tmp.left_delete[u]) continue;
				if (left_side_old_degree[u] != tmp.getV1Degree(u)) {
					left_side_old_degree[u] = tmp.getV1Degree(u);
					reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), left_side_old_degree[u], u);
					if (tmp.getV1Degree(u) == 0) {
						tmp.left_delete[u] = true;
					}
				}
			}
		}
		if (right_k == 1) {
			newtmp = tmp;
		}
	}
	tmp = newtmp;
	tmp.v1_max_degree = 0;
	tmp.v2_max_degree = 0;
	for (vid_t u = 0; u < tmp.degree_v1.size(); u++) {
		if (tmp.v1_max_degree < tmp.degree_v1[u]) tmp.v1_max_degree = tmp.degree_v1[u];
	}
	for (vid_t v = 0; v < tmp.degree_v2.size(); v++) {
		if (tmp.v2_max_degree < tmp.degree_v2[v]) tmp.v2_max_degree = tmp.degree_v2[v];
	}
}

void peel_with_fixed_left_k_acc_no_edge_deletion_no_tmp(int left_k, BiGraph& g) {
	vector<bool> left_deletion_next_round;
	vector<bool> right_deletion_next_round;
	vector<int> left_degree_next_round;
	vector<int> right_degree_next_round;
	vector<list<vid_t>> left_bin;
	vector<list<vid_t>::iterator> left_list_pointer;
	vector<list<vid_t>> right_bin;
	vector<list<vid_t>::iterator> right_list_pointer;
	binsort(left_bin, right_bin, left_list_pointer, right_list_pointer, g);
	// avoid too many reorder times
	vector<int> right_side_old_degree;
	for (vid_t v = 0; v < g.getV2Num(); v++) {
		right_side_old_degree.push_back(g.getV2Degree(v));
	}
	vector<int> left_side_old_degree;
	for (vid_t u = 0; u < g.getV1Num(); u++) {
		left_side_old_degree.push_back(g.getV1Degree(u));
	}
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		bool left_unchange = false;
		bool right_unchange = false;
		while (!left_unchange || !right_unchange) {
			left_unchange = true;
			right_unchange = true;
			// peel left
			vector<vid_t> reorder_right_vertices;
			for (int i = 1; i < left_k; i++) {
				for (auto j = left_bin[i].begin(); j != left_bin[i].end();) {
					vid_t u = *j;
					//vector<vid_t> neigh = tmp.getV1Neighbors(u);
					int old_degree = g.getV1Degree(u); // should equal to i
					for (int k = 0; k < g.neighbor_v1[u].size(); k++) {
						vid_t v = g.neighbor_v1[u][k];
						if (g.right_delete[v]) continue;
						int old_deg = g.getV2Degree(v); // should equal to tmp.getV2Degree(v) + 1
														//tmp.deleteEdge(u, v);
						g.degree_v2[v]--;
						g.degree_v1[u]--;
						reorder_right_vertices.push_back(v);
						//reorder_right(right_bin, right_list_pointer, tmp.getV2Degree(v), old_deg, v);
						// degree = 0 will not be considered in next round
						if (g.getV2Degree(v) == 0 && right_k - 1 > 0) {
							update_x_k_with_k_x(g, left_k, right_k - 1, v);
						}
					}
					g.left_delete[u] = true;
					if (g.getV1Degree(u) != 0) {
						int kkkkkk = 0;
					}
					++j; // j becomes unavailable after reordering
					reorder_left(left_bin, left_list_pointer, g.getV1Degree(u), old_degree, u);
					if (right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
					}
					left_unchange = false;
				}
			}
			for (int i = 0; i < reorder_right_vertices.size(); i++) {
				vid_t v = reorder_right_vertices[i];
				if (g.right_delete[v]) continue;
				if (right_side_old_degree[v] != g.getV2Degree(v)) {
					right_side_old_degree[v] = g.getV2Degree(v);
					reorder_right(right_bin, right_list_pointer, g.getV2Degree(v), right_side_old_degree[v], v);
					if (g.getV2Degree(v) == 0) {
						g.right_delete[v] = true;
					}
				}
			}
			// peel right
			vector<vid_t> reorder_left_vertices;
			for (int i = 1; i < right_k; i++) {
				for (auto j = right_bin[i].begin(); j != right_bin[i].end(); ) {
					vid_t v = *j;
					//vector<vid_t> neigh = tmp.getV2Neighbors(v);
					int old_degree = g.getV2Degree(v);
					for (int k = 0; k < g.neighbor_v2[v].size(); k++) {
						vid_t u = g.neighbor_v2[v][k];
						if (g.left_delete[u]) continue;
						int old_deg = g.getV1Degree(u);
						//tmp.deleteEdge(u, v);
						g.degree_v2[v]--;
						g.degree_v1[u]--;
						reorder_left_vertices.push_back(u);
						//reorder_left(left_bin, left_list_pointer, tmp.getV1Degree(u), old_deg, u);
						// degree = 0 will not be considered in next round
						if (g.getV1Degree(u) == 0 && right_k - 1 > 0) {
							g.left_index[u][left_k] = right_k - 1;
						}
					}
					g.right_delete[v] = true;
					if (g.getV2Degree(v) != 0) {
						int kkkkkk = 0;
					}
					++j;
					reorder_right(right_bin, right_list_pointer, g.getV2Degree(v), old_degree, v);
					if (right_k - 1 > 0) {
						update_x_k_with_k_x(g, left_k, right_k - 1, v);
					}
					right_unchange = false;
				}
			}
			for (int i = 0; i < reorder_left_vertices.size(); i++) {
				vid_t u = reorder_left_vertices[i];
				if (g.left_delete[u]) continue;
				if (left_side_old_degree[u] != g.getV1Degree(u)) {
					left_side_old_degree[u] = g.getV1Degree(u);
					reorder_left(left_bin, left_list_pointer, g.getV1Degree(u), left_side_old_degree[u], u);
					if (g.getV1Degree(u) == 0) {
						g.left_delete[u] = true;
					}
				}
			}
		}
		if (right_k == 1) {
			left_degree_next_round = g.degree_v1;
			right_degree_next_round = g.degree_v2;
			left_deletion_next_round = g.left_delete;
			right_deletion_next_round = g.right_delete;
		}
	}
	g.degree_v1 = left_degree_next_round;
	g.degree_v2 = right_degree_next_round;
	g.left_delete = left_deletion_next_round;
	g.right_delete = right_deletion_next_round;
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
}

void peel_with_fixed_left_k_acc_no_edge_deletion_no_tmp_no_bin_sort(int left_k, BiGraph& g) {
	vector<bool> left_deletion_next_round;
	vector<bool> right_deletion_next_round;
	vector<int> left_degree_next_round;
	vector<int> right_degree_next_round;
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;
	for (vid_t u = 0; u < g.getV1Num(); u++) {
		if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
			left_vertices_to_be_peeled.push_back(u);
		}
	}
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		bool ter = true;
		for (vid_t v = 0; v < g.getV2Num(); v++) {
			if (!g.right_delete[v]) ter = false;
			else continue;
			if (g.degree_v2[v] < right_k) {
				right_vertices_to_be_peeled.push_back(v);
			}
		}
		if (ter) break;
		while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
			// peel left
			for (auto j = left_vertices_to_be_peeled.begin(); j != left_vertices_to_be_peeled.end(); j++) {
				vid_t u = *j;
				if (g.left_delete[u]) continue;
				for (int k = 0; k < g.neighbor_v1[u].size(); k++) {
					vid_t v = g.neighbor_v1[u][k];
					if (g.right_delete[v]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.getV2Degree(v) == 0 && right_k - 1 > 0) {
						update_x_k_with_k_x(g, left_k, right_k - 1, v);
						g.right_delete[v] = true;
					}
					if (g.degree_v2[v] < right_k && !g.right_delete[v]) {
						right_vertices_to_be_peeled.push_back(v);
					}				
				}
				g.left_delete[u] = true;
				if (right_k - 1 > 0) {
					g.left_index[u][left_k] = right_k - 1;
				}
			}
			left_vertices_to_be_peeled.clear();
			// peel right
			for (auto j = right_vertices_to_be_peeled.begin(); j != right_vertices_to_be_peeled.end(); j++) {
				vid_t v = *j;
				if (g.right_delete[v]) continue;
				for (int k = 0; k < g.neighbor_v2[v].size(); k++) {
					vid_t u = g.neighbor_v2[v][k];
					if (g.left_delete[u]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.getV1Degree(u) == 0 && right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
						g.left_delete[u] = true;
					}
					if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.right_delete[v] = true;
				if (right_k - 1 > 0) {
					update_x_k_with_k_x(g, left_k, right_k - 1, v);
				}
			}
			right_vertices_to_be_peeled.clear();
		}
		if (right_k == 1) {
			left_degree_next_round = g.degree_v1;
			right_degree_next_round = g.degree_v2;
			left_deletion_next_round = g.left_delete;
			right_deletion_next_round = g.right_delete;
		}
	}
	g.degree_v1 = left_degree_next_round;
	g.degree_v2 = right_degree_next_round;
	g.left_delete = left_deletion_next_round;
	g.right_delete = right_deletion_next_round;
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
}

bool compare_degree_array(const pair<int, int>& p1, const pair<int, int>& p2) {
	return p1.second < p2.second;
}

void peel_with_fixed_left_k_acc_no_edge_deletion_dyn_array(int left_k, BiGraph& g) {
	vector<bool> left_deletion_next_round;
	vector<bool> right_deletion_next_round;
	vector<int> left_degree_next_round;
	vector<int> right_degree_next_round;
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;
	// definition of dyn array
	vector<int> right_k_x_index; right_k_x_index.resize(g.v2_max_degree + 2);
	fill_n(right_k_x_index.begin(), right_k_x_index.size(), -1); right_k_x_index[0] = 0; right_k_x_index[g.v2_max_degree + 1] = g.num_v2;
	vector<int> right_vertices_index; right_vertices_index.resize(g.num_v2);
	vector<pair<int, int>> degree_array; degree_array.resize(g.num_v2);
	// end of definition
	for (vid_t u = 0; u < g.getV1Num(); u++) {
		if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
			left_vertices_to_be_peeled.push_back(u);
		}
	}
	// initialize dyn array
	for (vid_t v = 0; v < g.num_v2; v++) {
		degree_array[v] = make_pair(v, g.degree_v2[v]);
	}
	sort(degree_array.begin(), degree_array.end(), compare_degree_array);
	for (int i = 0; i < degree_array.size(); i++) {
		if (right_k_x_index[degree_array[i].second] == -1) {
			right_k_x_index[degree_array[i].second] = i;
		}
		right_vertices_index[degree_array[i].first] = i;
	}
	for (int i = right_k_x_index.size() - 1; i >= 0; i--) {
		if (right_k_x_index[i] == -1) {
			right_k_x_index[i] = right_k_x_index[i + 1];
		}
	}
	// end of initialization
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		if (right_k_x_index[right_k - 1] == g.num_v2) break;
		for (int i = right_k_x_index[right_k - 1]; i < right_k_x_index[right_k]; i++) {
			right_vertices_to_be_peeled.push_back(degree_array[i].first);
		}
		while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
			// peel left
			for (auto j = left_vertices_to_be_peeled.begin(); j != left_vertices_to_be_peeled.end(); j++) {
				vid_t u = *j;
				if (g.left_delete[u]) continue;
				for (int k = 0; k < g.neighbor_v1[u].size(); k++) {
					vid_t v = g.neighbor_v1[u][k];
					if (g.right_delete[v]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.getV2Degree(v) == 0 && right_k - 1 > 0) {
						update_x_k_with_k_x(g, left_k, right_k - 1, v);
						g.right_delete[v] = true;
					}
					if (g.degree_v2[v] < right_k && !g.right_delete[v]) {
						right_vertices_to_be_peeled.push_back(v);
					}
					if (g.degree_v2[v] >= right_k - 1) {
						int olddegree = g.degree_v2[v] + 1;
						vid_t stack = degree_array[right_k_x_index[olddegree]].first;
						swap(degree_array[right_k_x_index[olddegree]], degree_array[right_vertices_index[v]]);
						swap(right_vertices_index[stack], right_vertices_index[v]);
						right_k_x_index[olddegree]++;
					}
				}
				g.left_delete[u] = true;
				if (right_k - 1 > 0) {
					g.left_index[u][left_k] = right_k - 1;
				}
			}
			left_vertices_to_be_peeled.clear();
			// peel right
			for (auto j = right_vertices_to_be_peeled.begin(); j != right_vertices_to_be_peeled.end(); j++) {
				vid_t v = *j;
				if (g.right_delete[v]) continue;
				for (int k = 0; k < g.neighbor_v2[v].size(); k++) {
					vid_t u = g.neighbor_v2[v][k];
					if (g.left_delete[u]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.getV1Degree(u) == 0 && right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
						g.left_delete[u] = true;
					}
					if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.right_delete[v] = true;
				if (right_k - 1 > 0) {
					update_x_k_with_k_x(g, left_k, right_k - 1, v);
				}
			}
			right_vertices_to_be_peeled.clear();
		}
		if (right_k == 1) {
			left_degree_next_round = g.degree_v1;
			right_degree_next_round = g.degree_v2;
			left_deletion_next_round = g.left_delete;
			right_deletion_next_round = g.right_delete;
		}
	}
	g.degree_v1 = left_degree_next_round;
	g.degree_v2 = right_degree_next_round;
	g.left_delete = left_deletion_next_round;
	g.right_delete = right_deletion_next_round;
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
}

void static_construct_index_with_peeling_algorithm_long_yuan_no_edge_deletion(BiGraph& g) {
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
	int right_max = 0;
	for (int left_k = 1; left_k <= g.v1_max_degree; left_k++) {
		peel_with_fixed_left_k_acc_no_edge_deletion_dyn_array(left_k, g);
		right_max = 0;
		for (vid_t u = 0; u < g.getV1Num(); u++) {
			if (g.left_index[u].size() <= left_k) continue;
			int right_k = g.left_index[u][left_k];
			if (right_max < right_k) right_max = right_k;
		}
		if (right_max <= left_k) break;
	}
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = 0;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = 0;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}

	cout << "right_max: " << right_max << endl;
	inverse(g);
	for (int left_k = 1; left_k <= right_max; left_k++) {
		peel_with_fixed_left_k_acc_no_edge_deletion_dyn_array(left_k, g);
	}
	inverse(g);
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = 0;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = 0;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
}

void static_construct_index_with_peeling_algorithm_acc_no_edge_deletion(BiGraph& g) {
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
	BiGraph tmp = g;
	for (int i = 1; i <= g.v1_max_degree; i++) {
		peel_with_fixed_left_k_acc_no_edge_deletion(i, tmp, g);
	}
}

//////////////long yuan's converange algorithm
//////////////////////////////////////
/////////////////////////////////////
////////////////////////////////////


void static_construct_index_with_convergence_algorithm_long_yuan(BiGraph& g) {
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
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	vector<int> right_tmp_k_x;
	vector<int> left_tmp_k_x;
	vector<int> right_score;
	vector<int> left_score;
	left_binns.resize(g.v2_max_degree + 1);
	right_binns.resize(g.v2_max_degree + 1);
	int right_max = 0;
	for (int left_k = 1; left_k <= g.v1_max_degree; left_k++) {
		converge_with_fixed_left_k_acc(left_k, g, left_tmp_k_x, right_tmp_k_x, left_score, right_score);
		right_max = 0;
		for (vid_t u = 0; u < g.getV1Num(); u++) {
			if (g.left_index[u].size() <= left_k) continue;
			int right_k = g.left_index[u][left_k];
			if (right_max < right_k) right_max = right_k;
		}
		if (right_max <= left_k) break;
	}
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = 0;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = 0;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}

	cout << "right_max: " << right_max << endl;
	inverse(g);
	right_tmp_k_x.clear(); left_tmp_k_x.clear(); right_score.clear(); left_score.clear(); left_binns.resize(g.v2_max_degree + 1); right_binns.resize(g.v2_max_degree + 1);
	for (int left_k = 1; left_k <= right_max; left_k++) {
		converge_with_fixed_left_k_acc(left_k, g, left_tmp_k_x, right_tmp_k_x, left_score, right_score);
	}
	inverse(g);
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = 0;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = 0;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
}