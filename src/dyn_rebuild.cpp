#include "dyn_rebuild.h"
#include <cmath>
using namespace std;
// use dynamic_test to test correctness
namespace Dyn_rebuild
{
	mutex g_mutex;
	volatile int alpha_mutex(0);
	volatile int beta_mutex(0);

	void dyn_crossUpdate_addition(BiGraph& g, int alpha, int k_x, vid_t v) {
		for (int beta = k_x; beta > 0; beta--) {
			int oldalpha = g.right_index[v][beta];
			if (oldalpha < alpha) {
				g.right_index[v][beta] = alpha;
			}
			else {
				break;
			}
		}
	}
	void dyn_crossUpdate_deletion(BiGraph& g, int alpha, int k_x, vid_t v) {
		int newalpha = alpha - 1;
		int truedegree = g.neighbor_v2[v].size();
		for (int i = k_x; i <= truedegree; i++) {
			int oldalpha = g.right_index[v][i];
			if (oldalpha > newalpha) {
				g.right_index[v][i] = newalpha;
			}
			else {
				break;
			}
		}
	}
	int calculate_ualpha_with_fixed_beta(BiGraph& g, vid_t u, int beta) {
		int ss = g.left_index[u].size();
		for (int alpha = 1; alpha < ss; alpha++) {
			if (g.left_index[u][alpha] < beta) {
				return alpha - 1;
			}
		}
		return ss - 1;
	}
	int calculate_vbeta_with_fixed_alpha(BiGraph& g, vid_t v, int alpha) {
		int ss = g.right_index[v].size();
		for (int beta = 1; beta < ss; beta++) {
			if (g.right_index[v][beta] < alpha) {
				return beta - 1;
			}
		}
		return ss - 1;
	}
	// dynamic update with swap
	pair<int, int> calculate_Delta(vector<vector<bicore_index_block*>>& bicore_index_u) {
		for (int alpha = 1; alpha < bicore_index_u.size(); alpha++) {
			if (bicore_index_u[alpha].size() - 1 <= alpha) {
				return make_pair(alpha, bicore_index_u[alpha].size() - 1);
			}
		}
		return make_pair(bicore_index_u.size(), 0);
	}

	pair<int, int> calculate_Delta(BiGraph& g) {
		int delta = 0;
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (g.left_index[u].size() > delta + 1) {
				for (int alpha = delta + 1; alpha < g.left_index[u].size(); alpha++) {
					if (g.left_index[u][alpha] >= alpha) {
						delta = alpha;
					}
					else {
						break;
					}
				}
			}
		}
		if (delta == 0) {
			return make_pair(1, 0);
		}
		else {
			return make_pair(delta, delta);
		}
	}

	int calculate_bound_for_left_node(BiGraph& g, vid_t u, int alpha) {
		// deletion operation when removing level
		if (alpha > g.degree_v1[u]) {
			return 0;
		}
		vector<int> k_index;
		vector<vid_t>& neigh = g.neighbor_v1[u];
		int ss = neigh.size();
		for (int i = 0; i < ss; i++) {
			int vbeta = calculate_vbeta_with_fixed_alpha(g, neigh[i], alpha);
			/*if (vbeta < 1) {
			cout << "error in calculate_bound_for_left_node" << endl;
			return -1;
			}*/
			k_index.push_back(vbeta);
		}
		sort(k_index.begin(), k_index.end());
		return k_index[k_index.size() - alpha];
	}

	int calculate_bound_for_right_node(BiGraph& g, vid_t v, int beta) {
		// deletion operation when removing level
		if (beta > g.degree_v2[v]) {
			return 0;
		}
		vector<int> k_index;
		vector<vid_t>& neigh = g.neighbor_v2[v];
		int ss = neigh.size();
		for (int i = 0; i < ss; i++) {
			int ualpha = calculate_ualpha_with_fixed_beta(g, neigh[i], beta);
			k_index.push_back(ualpha);
		}
		sort(k_index.begin(), k_index.end());
		return k_index[k_index.size() - beta];
	}

	void compute_a_b_core(BiGraph& g, int alpha, int beta) {
		// intialize g.left_delete and g.right_deletion
		//fill_n(g.left_delete.begin(), g.left_delete.size(), false);
		//fill_n(g.right_delete.begin(), g.right_delete.size(), false);
		vector<vid_t> left_vertices_to_be_peeled;
		vector<vid_t> right_vertices_to_be_peeled;
		bool stop = true;
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (!g.left_delete[u]) {
				if (g.degree_v1[u] < alpha) {
					left_vertices_to_be_peeled.push_back(u);
				}
			}
		}
		for (vid_t v = 0; v < g.num_v2; v++) {
			if (!g.right_delete[v]) {
				if (g.degree_v2[v] < beta) {
					right_vertices_to_be_peeled.push_back(v);
				}
			}
		}
		while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
			// peel left
			int oo_ = left_vertices_to_be_peeled.size();
			for (int j = 0; j < oo_; j++) {
				vid_t u = left_vertices_to_be_peeled[j];
				if (g.left_delete[u]) continue;
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int k = 0; k < ss; k++) {
					vid_t v = tmp_neigh_[k];
					if (g.right_delete[v]) continue;
					int dd_ = --g.degree_v2[v];
					if (dd_ == 0) {
						// core part
						g.right_delete[v] = true;
					}
					if (dd_ == beta - 1) {
						right_vertices_to_be_peeled.push_back(v);
					}
				}
				g.degree_v1[u] = 0;
				g.left_delete[u] = true;
			}
			left_vertices_to_be_peeled.clear();
			// peel right
			oo_ = right_vertices_to_be_peeled.size();
			for (int j = 0; j < oo_; j++) {
				vid_t v = right_vertices_to_be_peeled[j];
				if (g.right_delete[v]) continue;
				vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v];
				int ss = tmp_neigh_.size();
				for (int k = 0; k < ss; k++) {
					vid_t u = tmp_neigh_[k];
					if (g.left_delete[u]) continue;
					int dd_ = --g.degree_v1[u];
					if (dd_ == 0) {
						g.left_delete[u] = true;
					}
					if (dd_ == alpha - 1) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.degree_v2[v] = 0;
				g.right_delete[v] = true;
			}
			right_vertices_to_be_peeled.clear();
		}
	}
	void update_index_with_fixed_left_k_addition_with_limit_swap(BiGraph& g, int alpha, int tau_alpha, vid_t u, vid_t v) {
		if (tau_alpha > g.left_index[u][alpha]) {
			g.left_index[u][alpha] = tau_alpha;
		}
		int beta = tau_alpha;
		// compute alpha-beta+1-core
		compute_a_b_core(g, alpha, beta + 1);
		// update core value
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (!g.left_delete[u]) {
				int oldbeta = g.left_index[u][alpha];
				if (oldbeta < beta + 1) {
					g.left_index[u][alpha] = beta + 1;
				}
			}
		}
		for (vid_t v = 0; v < g.num_v2; v++) {
			if (!g.right_delete[v]) {
				dyn_crossUpdate_addition(g, alpha, beta + 1, v);
			}
		}
	}

	int upgrade_vbeta_insertion(BiGraph& g, vid_t v, int alpha, int beta) {
		int a1 = g.right_index[v][beta];
		if (beta == 0) {
			a1 = alpha;
		}
		int a2 = g.right_index[v][beta + 1];
		if (a1 < alpha) {
			return -1;
		}
		else if (a1 >= alpha && a2 < alpha) {
			return beta;
		}
		else {
			return 1e9;
		}
	}

	int upgrade_vbeta_removal(BiGraph& g, vid_t v, int alpha, int beta) {
		int a1 = g.right_index[v][beta];
		int a2 = -1;
		if (beta + 1 < g.right_index[v].size()) {
			a2 = g.right_index[v][beta + 1];
		}
		if (a1 < alpha) {
			return -1;
		}
		else if (a1 >= alpha && a2 < alpha) {
			return beta;
		}
		else {
			return 1e9;
		}
	}
	// may incur errors !!!
	int compute_vbeta_regarding_specific_alpha(BiGraph& g, vid_t v, int alpha) {
		int vbeta = -1;
		bool changed = false;
		for (int b = 1; b <= g.right_index[v].size() - 1; b++) {
			if (g.right_index[v][b] < alpha) {
				vbeta = b - 1;
				changed = true;
				break;
			}
		}
		if (!changed) {
			vbeta = g.right_index[v].size() - 1;
		}
		if (vbeta < 0) {
			//cout << "error: vbeta" << endl;
		}
		return vbeta;
	}
	int compute_vbeta_regarding_specific_alpha__(BiGraph& g, vid_t v, int alpha) {
		int vbeta = -1;
		bool changed = false;
		for (int b = 1; b <= g.degree_v2[v]; b++) {
			if (g.right_index[v][b] < alpha) {
				vbeta = b - 1;
				changed = true;
				break;
			}
		}
		if (!changed) {
			vbeta = g.degree_v2[v];
		}
		if (vbeta < 0) {
			cout << "error: vbeta" << endl;
		}
		return vbeta;
	}

	void casacade_remove_candidate_nodes(vid_t node, int alpha, int beta, bool left_side,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
		unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right,
		unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_right) {
		if (left_side) {
			candidate_nodes_left.erase(node);
			for (auto v : visited_nodes_2_info_left[node].neighbors_in_candidates) {
				if (candidate_nodes_right.find(v) != candidate_nodes_right.end()) {
					visited_nodes_2_info_right[v].current_supports--;
					if (visited_nodes_2_info_right[v].current_supports <= beta) {
						casacade_remove_candidate_nodes(v, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
		}
		else {
			candidate_nodes_right.erase(node);
			for (auto u : visited_nodes_2_info_right[node].neighbors_in_candidates) {
				if (candidate_nodes_left.find(u) != candidate_nodes_left.end()) {
					visited_nodes_2_info_left[u].current_supports--;
					if (visited_nodes_2_info_left[u].current_supports < alpha) {
						casacade_remove_candidate_nodes(u, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
		}
	}

	//void dyn_update_insertion_dfs(BiGraph& g, vid_t node, int alpha, int beta, bool left_side,
	//	unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
	//	unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right,
	//	unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_right) {
	//	// node + left_side * (g.num_v1 + g.num_v2) to distinguish between left and right
	//	vector<vid_t> nodes_to_be_visited;
	//	if (left_side) {
	//		candidate_nodes_left.insert(node);
	//		visited_nodes_left.insert(node);
	//		vector<vid_t>& neigh = g.neighbor_v1[node];
	//		for (auto v : neigh) {
	//			// could be a candidate only if v's deg > beta
	//			// aim at finding whether v belongs to alpha-beta+1-core
	//			if (g.degree_v2[v] > beta) {
	//				// could be optimized???
	//				int vbeta = upgrade_vbeta_insertion(g, v, alpha, beta);
	//				if (vbeta > beta) {
	//					visited_nodes_2_info_left[node].current_supports++;
	//				}
	//				else if (vbeta == beta) {
	//					if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
	//						visited_nodes_2_info_left[node].current_supports++;
	//						nodes_to_be_visited.push_back(v);
	//					}
	//					else if (visited_nodes_right.find(v) != visited_nodes_right.end() && candidate_nodes_right.find(v) != candidate_nodes_right.end()) {
	//						visited_nodes_2_info_left[node].current_supports++;
	//						visited_nodes_2_info_left[node].neighbors_in_candidates.push_back(v);
	//						visited_nodes_2_info_right[v].neighbors_in_candidates.push_back(node);
	//					}
	//					else {
	//						//do nothing
	//					}
	//				}
	//				else {
	//					// do nothing
	//				}
	//			}
	//		}
	//		if (visited_nodes_2_info_left[node].current_supports >= alpha) {
	//			for (auto v : nodes_to_be_visited) {
	//				// check whether need to continue
	//				if (candidate_nodes_left.find(node) == candidate_nodes_left.end()) {
	//					break;
	//				}
	//				//vid_t v = visited_nodes_2_info_left[node].neighbors_in_candidates[i];
	//				if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
	//					dyn_update_insertion_dfs(g, v, alpha, beta, !left_side,
	//						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
	//						visited_nodes_2_info_left, visited_nodes_2_info_right);
	//				}
	//			}
	//		}
	//		else {
	//			casacade_remove_candidate_nodes(node, alpha, beta, left_side,
	//				visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
	//				visited_nodes_2_info_left, visited_nodes_2_info_right);
	//		}
	//	}
	//	else {
	//		candidate_nodes_right.insert(node);
	//		visited_nodes_right.insert(node);
	//		vector<vid_t>& neigh = g.neighbor_v2[node];
	//		for (auto u : neigh) {
	//			if (g.degree_v1[u] >= alpha) {
	//				int ubeta = g.left_index[u][alpha];
	//				if (ubeta > beta) {
	//					visited_nodes_2_info_right[node].current_supports++;
	//				}
	//				else if (ubeta == beta) {
	//					if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
	//						visited_nodes_2_info_right[node].current_supports++;
	//						nodes_to_be_visited.push_back(u);
	//					}
	//					else if (visited_nodes_left.find(u) != visited_nodes_left.end() && candidate_nodes_left.find(u) != candidate_nodes_left.end()) {
	//						visited_nodes_2_info_right[node].current_supports++;
	//						visited_nodes_2_info_right[node].neighbors_in_candidates.push_back(u);
	//						visited_nodes_2_info_left[u].neighbors_in_candidates.push_back(node);
	//					}
	//					else {
	//						//do nothing
	//					}
	//				}
	//				else {
	//					// do nothing
	//				}
	//			}
	//		}
	//		if (visited_nodes_2_info_right[node].current_supports > beta) {
	//			for (auto u : nodes_to_be_visited) {
	//				// check whether need to continue
	//				if (candidate_nodes_right.find(node) == candidate_nodes_right.end()) {
	//					break;
	//				}
	//				//vid_t u = visited_nodes_2_info_right[node].neighbors_in_candidates[i];
	//				if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
	//					dyn_update_insertion_dfs(g, u, alpha, beta, !left_side,
	//						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
	//						visited_nodes_2_info_left, visited_nodes_2_info_right);
	//				}
	//			}
	//		}
	//		else {
	//			casacade_remove_candidate_nodes(node, alpha, beta, left_side,
	//				visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
	//				visited_nodes_2_info_left, visited_nodes_2_info_right);
	//		}
	//	}
	//}

	void dyn_update_insertion_dfs(BiGraph& g, vid_t node, int alpha, int beta, bool left_side,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
		unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right,
		unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_right) {
		// node + left_side * (g.num_v1 + g.num_v2) to distinguish between left and right
		vector<pair<pair<vid_t,vid_t>,bool>> stack; // pre cur left_side
		stack.push_back(make_pair(make_pair(-1,node), left_side));
		bool start_root = true;
		while (!stack.empty()) {
			auto back = stack.back();
			stack.pop_back();
			//vector<vid_t> nodes_to_be_visited;
			if (back.second) {
				vid_t pre = back.first.first;
				vid_t node = back.first.second;
				if (visited_nodes_left.find(node) != visited_nodes_left.end()) {
					continue;
				}
				if (candidate_nodes_right.find(pre) == candidate_nodes_right.end() && !start_root) {
					continue;
				}
				if (start_root) {
					start_root = false;
				}
				candidate_nodes_left.insert(node);
				visited_nodes_left.insert(node);
				vector<vid_t>& neigh = g.neighbor_v1[node];
				for (auto v : neigh) {
					// could be a candidate only if v's deg > beta
					// aim at finding whether v belongs to alpha-beta+1-core
					if (g.degree_v2[v] > beta) {
						// could be optimized???
						int vbeta = upgrade_vbeta_insertion(g, v, alpha, beta);
						if (vbeta > beta) {
							visited_nodes_2_info_left[node].current_supports++;
						}
						else if (vbeta == beta) {
							if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
								visited_nodes_2_info_left[node].current_supports++;
								stack.push_back(make_pair(make_pair(node,v),false));
							}
							else if (visited_nodes_right.find(v) != visited_nodes_right.end() && candidate_nodes_right.find(v) != candidate_nodes_right.end()) {
								visited_nodes_2_info_left[node].current_supports++;
								visited_nodes_2_info_left[node].neighbors_in_candidates.push_back(v);
								visited_nodes_2_info_right[v].neighbors_in_candidates.push_back(node);
							}
							else {
								//do nothing
							}
						}
						else {
							// do nothing
						}
					}
				}
				if (visited_nodes_2_info_left[node].current_supports >= alpha) {
					// do nothing
				}
				else {
					casacade_remove_candidate_nodes(node, alpha, beta, true,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
				}
			}
			else {
				vid_t pre = back.first.first;
				vid_t node = back.first.second;
				if (visited_nodes_right.find(node) != visited_nodes_right.end()) {
					continue;
				}
				if (candidate_nodes_left.find(pre) == candidate_nodes_left.end() && !start_root) {
					continue;
				}
				if (start_root) {
					start_root = false;
				}
				candidate_nodes_right.insert(node);
				visited_nodes_right.insert(node);
				vector<vid_t>& neigh = g.neighbor_v2[node];
				for (auto u : neigh) {
					if (g.degree_v1[u] >= alpha) {
						int ubeta = g.left_index[u][alpha];
						if (ubeta > beta) {
							visited_nodes_2_info_right[node].current_supports++;
						}
						else if (ubeta == beta) {
							if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
								visited_nodes_2_info_right[node].current_supports++;
								stack.push_back(make_pair(make_pair(node, u), true));
							}
							else if (visited_nodes_left.find(u) != visited_nodes_left.end() && candidate_nodes_left.find(u) != candidate_nodes_left.end()) {
								visited_nodes_2_info_right[node].current_supports++;
								visited_nodes_2_info_right[node].neighbors_in_candidates.push_back(u);
								visited_nodes_2_info_left[u].neighbors_in_candidates.push_back(node);
							}
							else {
								//do nothing
							}
						}
						else {
							// do nothing
						}
					}
				}
				if (visited_nodes_2_info_right[node].current_supports > beta) {
					
				}
				else {
					casacade_remove_candidate_nodes(node, alpha, beta, false,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
				}
			}
		}
	}

	void update_index_with_fixed_left_k_addition_swap_dfs(BiGraph& g, int alpha, int tau_alpha, vid_t u, vid_t v) {
		if (tau_alpha > g.left_index[u][alpha]) {
			g.left_index[u][alpha] = tau_alpha;
		}
		int beta = tau_alpha;
		unordered_set<vid_t> visited_nodes_left; unordered_set<vid_t> visited_nodes_right;
		unordered_set<vid_t> candidate_nodes_left; unordered_set<vid_t> candidate_nodes_right;
		unordered_map<vid_t, traversal_info_insertion> visited_nodes_2_info_left; unordered_map<vid_t, traversal_info_insertion> visited_nodes_2_info_right;
		// starting from u
		if (g.left_index[u][alpha] == beta) {
			dyn_update_insertion_dfs(g, u, alpha, beta, true, visited_nodes_left, visited_nodes_right,
				candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
		}
		// starting from v
		else {
			dyn_update_insertion_dfs(g, v, alpha, beta, false, visited_nodes_left, visited_nodes_right,
				candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
		}
		for (auto u : candidate_nodes_left) {
			g.left_index[u][alpha] = beta + 1;
		}
		for (auto v : candidate_nodes_right) {
			dyn_crossUpdate_addition(g, alpha, beta + 1, v);
		}
	}

	void casacade_insert_candidate_nodes(BiGraph& g, vid_t node, int alpha, int beta, bool left_side,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
		unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right,
		unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_right) {
		if (left_side) {
			candidate_nodes_left.insert(node);
			for (auto v : visited_nodes_2_info_left[node].neighbors_not_in_candidates) {
				if (candidate_nodes_right.find(v) == candidate_nodes_right.end()) {
					visited_nodes_2_info_right[v].current_supports--;
					if (visited_nodes_2_info_right[v].current_supports < beta) {
						casacade_insert_candidate_nodes(g, v, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
			for (auto v : visited_nodes_2_info_left[node].nodes_to_be_visited) {
				if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
					dyn_update_removal_dfs(g, v, alpha, beta, !left_side,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
				}
			}
		}
		else {
			candidate_nodes_right.insert(node);
			for (auto u : visited_nodes_2_info_right[node].neighbors_not_in_candidates) {
				if (candidate_nodes_left.find(u) == candidate_nodes_left.end()) {
					visited_nodes_2_info_left[u].current_supports--;
					if (visited_nodes_2_info_left[u].current_supports < alpha) {
						casacade_insert_candidate_nodes(g, u, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
			for (auto u : visited_nodes_2_info_right[node].nodes_to_be_visited) {
				if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
					dyn_update_removal_dfs(g, u, alpha, beta, !left_side,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
				}
			}
		}
	}

	void dyn_update_removal_dfs(BiGraph& g, vid_t node, int alpha, int beta, bool left_side,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
		unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right,
		unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_right) {
		if (!left_side && node == 4) {
			int kkkkkk = 0;// for debug
		}
		if (left_side) {
			visited_nodes_left.insert(node);
			vector<vid_t>& neigh = g.neighbor_v1[node];
			for (auto v : neigh) {
				// could be a candidate only if v's deg >= beta
				// aim at finding whether v still belongs to alpha-beta-core
				if (g.degree_v2[v] >= beta) {
					// could be optimized???
					int vbeta = upgrade_vbeta_removal(g, v, alpha, beta);
					if (vbeta > beta) {
						visited_nodes_2_info_left[node].current_supports++;
					}
					else if (vbeta == beta) {
						if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
							visited_nodes_2_info_left[node].current_supports++;
							visited_nodes_2_info_left[node].nodes_to_be_visited.push_back(v);
						}
						else if (visited_nodes_right.find(v) != visited_nodes_right.end() && candidate_nodes_right.find(v) == candidate_nodes_right.end()) {
							visited_nodes_2_info_left[node].current_supports++;
							visited_nodes_2_info_left[node].neighbors_not_in_candidates.push_back(v);
							visited_nodes_2_info_right[v].neighbors_not_in_candidates.push_back(node);
						}
						else {
							//do nothing
						}
					}
					else {
						// do nothing
					}
				}
			}
			// expand current node
			if (visited_nodes_2_info_left[node].current_supports < alpha) {
				casacade_insert_candidate_nodes(g, node, alpha, beta, left_side,
					visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
					visited_nodes_2_info_left, visited_nodes_2_info_right);
				for (auto v : visited_nodes_2_info_left[node].nodes_to_be_visited) {
					if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
						dyn_update_removal_dfs(g, v, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
		}
		else {
			visited_nodes_right.insert(node);
			vector<vid_t>& neigh = g.neighbor_v2[node];
			for (auto u : neigh) {
				if (g.degree_v1[u] >= alpha) {
					int ubeta = g.left_index[u][alpha];
					if (ubeta > beta) {
						visited_nodes_2_info_right[node].current_supports++;
					}
					else if (ubeta == beta) {
						if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
							visited_nodes_2_info_right[node].current_supports++;
							visited_nodes_2_info_right[node].nodes_to_be_visited.push_back(u);
						}
						// can be simplified
						else if (visited_nodes_left.find(u) != visited_nodes_left.end() && candidate_nodes_left.find(u) == candidate_nodes_left.end()) {
							visited_nodes_2_info_right[node].current_supports++;
							// this node is doomed! must be the end point of insertion edge otherwise it can not be visited!
							if (g.right_index[node][beta] < alpha || g.degree_v2[node] < beta) continue;
							visited_nodes_2_info_right[node].neighbors_not_in_candidates.push_back(u);
							visited_nodes_2_info_left[u].neighbors_not_in_candidates.push_back(node);
						}
						else {
							//do nothing
						}
					}
					else {
						// do nothing
					}
				}
			}
			if (visited_nodes_2_info_right[node].current_supports < beta) {
				casacade_insert_candidate_nodes(g, node, alpha, beta, left_side,
					visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
					visited_nodes_2_info_left, visited_nodes_2_info_right);
				for (auto u : visited_nodes_2_info_right[node].nodes_to_be_visited) {
					if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
						dyn_update_removal_dfs(g, u, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
		}
	}
	void update_index_with_fixed_left_k_deletion_swap_dfs(BiGraph& g, int alpha, int tau_alpha, int bound, vid_t u, vid_t v) {
		int beta = tau_alpha;
		unordered_set<vid_t> visited_nodes_left; unordered_set<vid_t> visited_nodes_right;
		unordered_set<vid_t> candidate_nodes_left; unordered_set<vid_t> candidate_nodes_right;
		unordered_map<vid_t, traversal_info_removal> visited_nodes_2_info_left; unordered_map<vid_t, traversal_info_removal> visited_nodes_2_info_right;
		// starting from u
		if (g.left_index[u][alpha] == beta) {
			dyn_update_removal_dfs(g, u, alpha, beta, true, visited_nodes_left, visited_nodes_right,
				candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
		}
		// compute vbeta before removing edge
		int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
		// starting from v, vbeta may become smaller than beta during updating process
		// avoid duplicate visiting v if v is already in candidates
		if (vbeta <= beta && candidate_nodes_right.find(v) == candidate_nodes_right.end()) {
			dyn_update_removal_dfs(g, v, alpha, beta, false, visited_nodes_left, visited_nodes_right,
				candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
		}
		for (auto u : candidate_nodes_left) {
			g.left_index[u][alpha] = beta - 1;
		}
		for (auto v : candidate_nodes_right) {
			dyn_crossUpdate_deletion(g, alpha, beta, v);
		}
		// contrary to insertion
		if (g.left_index[u][alpha] > bound) {
			g.left_index[u][alpha] = bound;
		}
	}
	double update_bicore_index_without_limit_swap(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition)
	{
		return 0;
	}
	void update_index_with_fixed_left_k_deletion_with_limit_swap(BiGraph& g, int alpha, int tau_alpha, int start_bound, vid_t u, vid_t v) {
		int beta = tau_alpha;

		compute_a_b_core(g, alpha, beta);

		for (vid_t u = 0; u < g.num_v1; u++) {
			if (g.left_delete[u] && g.left_index[u].size() >= alpha + 1) {
				int oldbeta = g.left_index[u][alpha];
				if (oldbeta == beta) {
					g.left_index[u][alpha] = beta - 1;
				}
			}
		}
		for (vid_t v = 0; v < g.num_v2; v++) {
			if (g.right_delete[v] && g.right_index[v].size() >= beta + 1) {
				dyn_crossUpdate_deletion(g, alpha, beta, v);
			}
		}
		// contrary to insertion
		int bound = start_bound;
		if (g.left_index[u][alpha] > bound) {
			g.left_index[u][alpha] = bound;
		}
	}
	double update_bicore_index_with_limit_swap(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition)
	{
		auto start = chrono::system_clock::now();
		// calculate swap threshold
		pair<int, int> r = calculate_Delta(g);
		// check whether operation is legal and insert/remove the edge
		if (addition) {
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						cout << "illegal insertion 1" << endl;
						return 0;
					}
				}
				g.addEdge(u, v);
			}
			else {
				cout << "illegal insertion 2" << endl;
				return 0;
			}
		}
		else {
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
		// addition
		if (addition) {
			// extend core index and bicore index and set threshold
			g.left_index[u].push_back(0);
			g.right_index[v].push_back(0);
			int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
			int th_alpha, th_beta;
			if (r.first + 1 + r.second + 1 >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + 1 + r.second + 1 >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first + 1; th_beta = r.second + 1;
			}
			// compute alpha beta condition in advance
			vector<int> tau_alpha; vector<int> tau_beta;
			tau_alpha.push_back(-1); tau_beta.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v];
				}
				if (vbeta < 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = bound > vbeta ? vbeta : bound;
				tau_alpha.push_back(beta);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u];
				}
				if (ualpha < 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = bound > ualpha ? ualpha : bound;
				tau_beta.push_back(alpha);
			}
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				update_index_with_fixed_left_k_addition_with_limit_swap(g, alpha, tau_alpha[alpha], u, v);
				// restore graph
				int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
				fill_n(g.left_delete.begin(), g.left_delete.size(), false);
				for (vid_t u = 0; u < g.num_v1; u++) {
					g.degree_v1[u] = g.neighbor_v1[u].size();
					left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
				}
				g.v1_max_degree = left_degree_max;
				fill_n(g.right_delete.begin(), g.right_delete.size(), false);
				for (vid_t v = 0; v < g.num_v2; v++) {
					g.degree_v2[v] = g.neighbor_v2[v].size();
					right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
				}
				g.v2_max_degree = right_degree_max;
			}
			inv(g);
			for (int alpha = 1; alpha <= th_beta; alpha++) {
				update_index_with_fixed_left_k_addition_with_limit_swap(g, alpha, tau_beta[alpha], v, u);
				// restore graph
				int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
				fill_n(g.left_delete.begin(), g.left_delete.size(), false);
				for (vid_t u = 0; u < g.num_v1; u++) {
					g.degree_v1[u] = g.neighbor_v1[u].size();
					left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
				}
				g.v1_max_degree = left_degree_max;
				fill_n(g.right_delete.begin(), g.right_delete.size(), false);
				for (vid_t v = 0; v < g.num_v2; v++) {
					g.degree_v2[v] = g.neighbor_v2[v].size();
					right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
				}
				g.v2_max_degree = right_degree_max;
			}
			inv(g);
		}
		// deletion
		else {
			int oldusat = g.left_index[u].back();
			int oldvsat = g.right_index[v].back();
			int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
			int th_alpha, th_beta;
			if (r.first + r.second >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + r.second >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first; th_beta = r.second;
			}
			// compute alpha beta condition in advance
			vector<int> alpha_start; vector<int> beta_start; vector<int> alpha_start_bound; vector<int> beta_start_bound;
			alpha_start.push_back(-1); beta_start.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				// g.degree_v2[v] has already been decreased
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v] + 1;
				}
				if (vbeta == -1 || vbeta == 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = ubeta > vbeta ? vbeta : ubeta;
				alpha_start.push_back(beta);
				alpha_start_bound.push_back(bound);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u] + 1;
				}
				if (ualpha == -1 || ualpha == 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = valpha > ualpha ? ualpha : valpha;
				beta_start.push_back(alpha);
				beta_start_bound.push_back(bound);
			}
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				update_index_with_fixed_left_k_deletion_with_limit_swap(g, alpha, alpha_start[alpha], alpha_start_bound[alpha], u, v);
				// restore graph
				int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
				fill_n(g.left_delete.begin(), g.left_delete.size(), false);
				for (vid_t u = 0; u < g.num_v1; u++) {
					g.degree_v1[u] = g.neighbor_v1[u].size();
					left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
				}
				g.v1_max_degree = left_degree_max;
				fill_n(g.right_delete.begin(), g.right_delete.size(), false);
				for (vid_t v = 0; v < g.num_v2; v++) {
					g.degree_v2[v] = g.neighbor_v2[v].size();
					right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
				}
				g.v2_max_degree = right_degree_max;
			}
			inv(g);
			for (int alpha = 1; alpha <= th_beta; alpha++) {
				update_index_with_fixed_left_k_deletion_with_limit_swap(g, alpha, beta_start[alpha], beta_start_bound[alpha], v, u);
				// restore graph
				int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
				fill_n(g.left_delete.begin(), g.left_delete.size(), false);
				for (vid_t u = 0; u < g.num_v1; u++) {
					g.degree_v1[u] = g.neighbor_v1[u].size();
					left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
				}
				g.v1_max_degree = left_degree_max;
				fill_n(g.right_delete.begin(), g.right_delete.size(), false);
				for (vid_t v = 0; v < g.num_v2; v++) {
					g.degree_v2[v] = g.neighbor_v2[v].size();
					right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
				}
				g.v2_max_degree = right_degree_max;
			}
			inv(g);
			g.left_index[u].pop_back();
			g.right_index[v].pop_back();
		}
		auto end = chrono::system_clock::now();
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		chrono::duration<double> elapsed_seconds = end - start;
		return elapsed_seconds.count();
	}

	double update_bicore_index_with_limit_swap_base_opt(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition)
	{
		auto start = chrono::system_clock::now();
		// calculate swap threshold
		pair<int, int> r = calculate_Delta(g);
		// check whether operation is legal and insert/remove the edge
		if (addition) {
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						cout << "illegal insertion 1" << endl;
						return 0;
					}
				}
				g.addEdge(u, v);
			}
			else {
				cout << "illegal insertion 2" << endl;
				return 0;
			}
		}
		else {
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
		// addition
		if (addition) {
			// extend core index and bicore index and set threshold
			g.left_index[u].push_back(0);
			g.right_index[v].push_back(0);
			int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
			int th_alpha, th_beta;
			if (r.first + 1 + r.second + 1 >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + 1 + r.second + 1 >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first + 1; th_beta = r.second + 1;
			}
			// compute alpha beta condition in advance
			vector<int> tau_alpha; vector<int> tau_beta; vector<int> alpha_start_bound; vector<int> beta_start_bound;
			tau_alpha.push_back(-1); tau_beta.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v];
				}
				if (vbeta < 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = bound > vbeta ? vbeta : bound;
				tau_alpha.push_back(beta);
				alpha_start_bound.push_back(bound);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u];
				}
				if (ualpha < 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = bound > ualpha ? ualpha : bound;
				tau_beta.push_back(alpha);
				beta_start_bound.push_back(bound);
			}
			// compute a base core to reduce computation
			int base_alpha = 2;
			auto left_delete = g.left_delete; auto right_delete = g.right_delete;
			auto left_degree = g.degree_v1; auto right_degree = g.degree_v2;
			auto left_max_degree = g.v1_max_degree; auto right_max_degree = g.v2_max_degree;
			////
			compute_a_b_core(g, base_alpha, tau_alpha[th_alpha]);
			////
			swap(left_delete, g.left_delete); swap(right_delete, g.right_delete);
			swap(left_degree, g.degree_v1); swap(right_degree, g.degree_v2);
			swap(left_max_degree, g.v1_max_degree); swap(right_max_degree, g.v2_max_degree);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				update_index_with_fixed_left_k_addition_with_limit_swap(g, alpha, tau_alpha[alpha], u, v);
				// restore graph
				if (alpha < base_alpha - 1 || alpha == th_alpha) {
					int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
					fill_n(g.left_delete.begin(), g.left_delete.size(), false);
					for (vid_t u = 0; u < g.num_v1; u++) {
						g.degree_v1[u] = g.neighbor_v1[u].size();
						left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
					}
					g.v1_max_degree = left_degree_max;
					fill_n(g.right_delete.begin(), g.right_delete.size(), false);
					for (vid_t v = 0; v < g.num_v2; v++) {
						g.degree_v2[v] = g.neighbor_v2[v].size();
						right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
					}
					g.v2_max_degree = right_degree_max;
				}
				else {
					g.left_delete = left_delete; g.right_delete = right_delete;
					g.degree_v1 = left_degree; g.degree_v2 = right_degree;
					g.v1_max_degree = left_max_degree; g.v2_max_degree = right_max_degree;
				}
			}
			inv(g);
			// compute a base core to reduce computation
			left_delete = g.left_delete; right_delete = g.right_delete;
			left_degree = g.degree_v1; right_degree = g.degree_v2;
			left_max_degree = g.v1_max_degree; right_max_degree = g.v2_max_degree;
			////
			compute_a_b_core(g, base_alpha, tau_beta[th_beta]);
			////
			swap(left_delete, g.left_delete); swap(right_delete, g.right_delete);
			swap(left_degree, g.degree_v1); swap(right_degree, g.degree_v2);
			swap(left_max_degree, g.v1_max_degree); swap(right_max_degree, g.v2_max_degree);
			for (int alpha = 1; alpha <= th_beta; alpha++) {
				update_index_with_fixed_left_k_addition_with_limit_swap(g, alpha, tau_beta[alpha], v, u);
				// restore graph
				if (alpha < base_alpha - 1 || alpha == th_beta) {
					int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
					fill_n(g.left_delete.begin(), g.left_delete.size(), false);
					for (vid_t u = 0; u < g.num_v1; u++) {
						g.degree_v1[u] = g.neighbor_v1[u].size();
						left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
					}
					g.v1_max_degree = left_degree_max;
					fill_n(g.right_delete.begin(), g.right_delete.size(), false);
					for (vid_t v = 0; v < g.num_v2; v++) {
						g.degree_v2[v] = g.neighbor_v2[v].size();
						right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
					}
					g.v2_max_degree = right_degree_max;
				}
				else {
					g.left_delete = left_delete; g.right_delete = right_delete;
					g.degree_v1 = left_degree; g.degree_v2 = right_degree;
					g.v1_max_degree = left_max_degree; g.v2_max_degree = right_max_degree;
				}

			}
			inv(g);
		}
		// deletion
		else {
			int oldusat = g.left_index[u].back();
			int oldvsat = g.right_index[v].back();
			int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
			int th_alpha, th_beta;
			if (r.first + r.second >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + r.second >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first; th_beta = r.second;
			}
			// compute alpha beta condition in advance
			vector<int> tau_alpha; vector<int> tau_beta; vector<int> alpha_start_bound; vector<int> beta_start_bound;
			tau_alpha.push_back(-1); tau_beta.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				// g.degree_v2[v] has already been decreased
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v] + 1;
				}
				if (vbeta == -1 || vbeta == 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = ubeta > vbeta ? vbeta : ubeta;
				tau_alpha.push_back(beta);
				alpha_start_bound.push_back(bound);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u] + 1;
				}
				if (ualpha == -1 || ualpha == 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = valpha > ualpha ? ualpha : valpha;
				tau_beta.push_back(alpha);
				beta_start_bound.push_back(bound);
			}
			// compute a base core to reduce computation
			int base_alpha = 2;
			auto left_delete = g.left_delete; auto right_delete = g.right_delete;
			auto left_degree = g.degree_v1; auto right_degree = g.degree_v2;
			auto left_max_degree = g.v1_max_degree; auto right_max_degree = g.v2_max_degree;
			////
			compute_a_b_core(g, base_alpha, tau_alpha[th_alpha]);
			////
			swap(left_delete, g.left_delete); swap(right_delete, g.right_delete);
			swap(left_degree, g.degree_v1); swap(right_degree, g.degree_v2);
			swap(left_max_degree, g.v1_max_degree); swap(right_max_degree, g.v2_max_degree);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				update_index_with_fixed_left_k_deletion_with_limit_swap(g, alpha, tau_alpha[alpha], alpha_start_bound[alpha], u, v);
				// restore graph
				if (alpha < base_alpha - 1 || alpha == th_alpha) {
					int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
					fill_n(g.left_delete.begin(), g.left_delete.size(), false);
					for (vid_t u = 0; u < g.num_v1; u++) {
						g.degree_v1[u] = g.neighbor_v1[u].size();
						left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
					}
					g.v1_max_degree = left_degree_max;
					fill_n(g.right_delete.begin(), g.right_delete.size(), false);
					for (vid_t v = 0; v < g.num_v2; v++) {
						g.degree_v2[v] = g.neighbor_v2[v].size();
						right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
					}
					g.v2_max_degree = right_degree_max;
				}
				else {
					g.left_delete = left_delete; g.right_delete = right_delete;
					g.degree_v1 = left_degree; g.degree_v2 = right_degree;
					g.v1_max_degree = left_max_degree; g.v2_max_degree = right_max_degree;
				}
			}
			inv(g);
			// compute a base core to reduce computation
			left_delete = g.left_delete; right_delete = g.right_delete;
			left_degree = g.degree_v1; right_degree = g.degree_v2;
			left_max_degree = g.v1_max_degree; right_max_degree = g.v2_max_degree;
			////
			compute_a_b_core(g, base_alpha, tau_beta[th_beta]);
			////
			swap(left_delete, g.left_delete); swap(right_delete, g.right_delete);
			swap(left_degree, g.degree_v1); swap(right_degree, g.degree_v2);
			swap(left_max_degree, g.v1_max_degree); swap(right_max_degree, g.v2_max_degree);
			for (int alpha = 1; alpha <= th_beta; alpha++) {
				update_index_with_fixed_left_k_deletion_with_limit_swap(g, alpha, tau_beta[alpha], beta_start_bound[alpha], v, u);
				// restore graph
				if (alpha < base_alpha - 1 || alpha == th_beta) {
					int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
					fill_n(g.left_delete.begin(), g.left_delete.size(), false);
					for (vid_t u = 0; u < g.num_v1; u++) {
						g.degree_v1[u] = g.neighbor_v1[u].size();
						left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
					}
					g.v1_max_degree = left_degree_max;
					fill_n(g.right_delete.begin(), g.right_delete.size(), false);
					for (vid_t v = 0; v < g.num_v2; v++) {
						g.degree_v2[v] = g.neighbor_v2[v].size();
						right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
					}
					g.v2_max_degree = right_degree_max;
				}
				else {
					g.left_delete = left_delete; g.right_delete = right_delete;
					g.degree_v1 = left_degree; g.degree_v2 = right_degree;
					g.v1_max_degree = left_max_degree; g.v2_max_degree = right_max_degree;
				}
			}
			inv(g);
			g.left_index[u].pop_back();
			g.right_index[v].pop_back();
		}
		auto end = chrono::system_clock::now();
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		chrono::duration<double> elapsed_seconds = end - start;
		return elapsed_seconds.count();
	}

	double update_bicore_index_swap_with_dfs(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition) {
		auto start = chrono::system_clock::now();
		// calculate swap threshold
		pair<int, int> r = calculate_Delta(g);
		// check whether operation is legal and insert/remove the edge
		if (addition) {
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						cout << "illegal insertion 1" << endl;
						return 0;
					}
				}
				g.addEdge(u, v);
			}
			else {
				cout << "illegal insertion 2" << endl;
				return 0;
			}
		}
		else {
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
		if (addition) {
			// extend core index and bicore index and set threshold
			g.left_index[u].push_back(0);
			g.right_index[v].push_back(0);
			int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
			int th_alpha, th_beta;
			if (r.first + 1 + r.second + 1 >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + 1 + r.second + 1 >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first + 1; th_beta = r.second + 1;
			}
			// compute alpha beta condition in advance
			vector<int> tau_alpha; vector<int> tau_beta; vector<int> alpha_start_bound; vector<int> beta_start_bound;
			tau_alpha.push_back(-1); tau_beta.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v];
				}
				if (vbeta < 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = bound > vbeta ? vbeta : bound;
				tau_alpha.push_back(beta);
				alpha_start_bound.push_back(bound);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u];
				}
				if (ualpha < 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = bound > ualpha ? ualpha : bound;
				tau_beta.push_back(alpha);
				beta_start_bound.push_back(bound);
			}
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				update_index_with_fixed_left_k_addition_swap_dfs(g, alpha, tau_alpha[alpha], u, v);
			}
			inv(g);
			for (int alpha = 1; alpha <= th_beta; alpha++) {
				update_index_with_fixed_left_k_addition_swap_dfs(g, alpha, tau_beta[alpha], v, u);
			}
			inv(g);
		}
		else {
			int oldusat = g.left_index[u].back();
			int oldvsat = g.right_index[v].back();
			int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
			int th_alpha, th_beta;
			if (r.first + r.second >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + r.second >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first; th_beta = r.second;
			}
			// compute alpha beta condition in advance
			vector<int> alpha_start; vector<int> beta_start; vector<int> alpha_start_bound; vector<int> beta_start_bound;
			alpha_start.push_back(-1); beta_start.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				// g.degree_v2[v] has already been decreased
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v] + 1;
				}
				if (vbeta == -1 || vbeta == 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = ubeta > vbeta ? vbeta : ubeta;
				alpha_start.push_back(beta);
				alpha_start_bound.push_back(bound);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u] + 1;
				}
				if (ualpha == -1 || ualpha == 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = valpha > ualpha ? ualpha : valpha;
				beta_start.push_back(alpha);
				beta_start_bound.push_back(bound);
			}
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				update_index_with_fixed_left_k_deletion_swap_dfs(g, alpha, alpha_start[alpha], alpha_start_bound[alpha], u, v);
			}
			inv(g);
			for (int alpha = 1; alpha <= th_beta; alpha++) {
				update_index_with_fixed_left_k_deletion_swap_dfs(g, alpha, beta_start[alpha], beta_start_bound[alpha], v, u);
			}
			inv(g);
			g.left_index[u].pop_back();
			g.right_index[v].pop_back();
		}
		auto end = chrono::system_clock::now();
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		chrono::duration<double> elapsed_seconds = end - start;
		return elapsed_seconds.count();
	}


	int maxkcorenumber_test_optimal(BiGraph& g) {
		vector<vid_t> left_vertices_to_be_peeled;
		vector<vid_t> right_vertices_to_be_peeled;

		int left_remain_nodes_num = g.num_v1;
		vector<vid_t> left_remain_nodes; left_remain_nodes.resize(g.num_v1);
		for (int i = 0; i < left_remain_nodes.size(); i++) {
			left_remain_nodes[i] = i;
		}
		int left_remain_nodes_tmp_num = 0;
		vector<vid_t> left_remain_nodes_tmp; left_remain_nodes_tmp.resize(g.num_v1);

		int right_remain_nodes_num = g.num_v2;
		vector<vid_t> right_remain_nodes; right_remain_nodes.resize(g.num_v2);
		for (int i = 0; i < right_remain_nodes.size(); i++) {
			right_remain_nodes[i] = i;
		}
		int right_remain_nodes_tmp_num = 0;
		vector<vid_t> right_remain_nodes_tmp; right_remain_nodes_tmp.resize(g.num_v2);

		int kc = 1;
		for (kc = 1; kc <= g.v2_max_degree + 1; kc++) {

			bool stop = true;
			left_remain_nodes_tmp_num = 0;
			for (int i = 0; i < left_remain_nodes_num; i++) {
				vid_t u = left_remain_nodes[i];
				if (!g.left_delete[u]) {
					stop = false;
					left_remain_nodes_tmp[left_remain_nodes_tmp_num] = u;
					left_remain_nodes_tmp_num++;
					if (g.degree_v1[u] < kc) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
			}
			swap(left_remain_nodes, left_remain_nodes_tmp);
			left_remain_nodes_num = left_remain_nodes_tmp_num;
			if (stop) break;

			stop = true;
			right_remain_nodes_tmp_num = 0;
			for (int i = 0; i < right_remain_nodes_num; i++) {
				vid_t v = right_remain_nodes[i];
				if (!g.right_delete[v]) {
					stop = false;
					right_remain_nodes_tmp[right_remain_nodes_tmp_num] = v;
					right_remain_nodes_tmp_num++;
					if (g.degree_v2[v] < kc) {
						right_vertices_to_be_peeled.push_back(v);
					}
				}
			}
			swap(right_remain_nodes, right_remain_nodes_tmp);
			right_remain_nodes_num = right_remain_nodes_tmp_num;
			if (stop) break;

			while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
				// peel left
				for (auto j = left_vertices_to_be_peeled.begin(); j != left_vertices_to_be_peeled.end(); j++) {
					vid_t u = *j;
					if (g.left_delete[u]) continue;
					for (int k = 0; k < g.neighbor_v1[u].size(); k++) {
						vid_t v = g.neighbor_v1[u][k];
						if (g.right_delete[v]) continue;
						g.degree_v2[v]--;
						//g.degree_v1[u]--;
						if (g.degree_v2[v] == 0) {
							g.right_delete[v] = true;
						}
						if (g.degree_v2[v] < kc) {
							right_vertices_to_be_peeled.push_back(v);
						}
					}
					g.degree_v1[u] = 0;
					g.left_delete[u] = true;
				}
				left_vertices_to_be_peeled.clear();
				// peel right
				for (auto j = right_vertices_to_be_peeled.begin(); j != right_vertices_to_be_peeled.end(); j++) {
					vid_t v = *j;
					if (g.right_delete[v]) continue;
					for (int k = 0; k < g.neighbor_v2[v].size(); k++) {
						vid_t u = g.neighbor_v2[v][k];
						if (g.left_delete[u]) continue;
						//g.degree_v2[v]--;
						g.degree_v1[u]--;
						if (g.degree_v1[u] == 0) {
							g.left_delete[u] = true;
						}
						if (g.degree_v1[u] < kc) {
							left_vertices_to_be_peeled.push_back(u);
						}
					}
					g.degree_v2[v] = 0;
					g.right_delete[v] = true;
				}
				right_vertices_to_be_peeled.clear();
			}
		}
		// restore g
		fill_n(g.left_delete.begin(), g.left_delete.size(), false);
		for (vid_t u = 0; u < g.num_v1; u++) {
			g.degree_v1[u] = g.neighbor_v1[u].size();
		}
		fill_n(g.right_delete.begin(), g.right_delete.size(), false);
		for (vid_t v = 0; v < g.num_v2; v++) {
			g.degree_v2[v] = g.neighbor_v2[v].size();
		}
		return kc - 2;
	}

	void batch_update_index_with_fixed_alpha_addition(BiGraph& g, int alpha, int pi_alpha) {
		int dd_;
		int pre_left_k_ = alpha - 1;
		vector<vid_t> left_vertices_to_be_peeled;
		vector<vid_t> right_vertices_to_be_peeled;
		for (vid_t u = 0; u < g.getV1Num(); u++) {
			if (g.degree_v1[u] < alpha && !g.left_delete[u]) {
				left_vertices_to_be_peeled.push_back(u);
			}
		}
		int right_remain_nodes_num = g.num_v2;
		vector<vid_t> right_remain_nodes; right_remain_nodes.resize(g.num_v2);
		for (int i = 0; i < right_remain_nodes.size(); i++) {
			right_remain_nodes[i] = i;
		}
		int right_remain_nodes_tmp_num = 0;
		vector<vid_t> right_remain_nodes_tmp; right_remain_nodes_tmp.resize(g.num_v2);
		bool update_flag = false;
		for (int right_k = pi_alpha; right_k <= g.v2_max_degree + 1; right_k++) {
			if (right_k > pi_alpha) {
				update_flag = true;
			}
			int pre_ = right_k - 1;
			bool stop = true;
			right_remain_nodes_tmp_num = 0;
			for (int i = 0; i < right_remain_nodes_num; i++) {
				vid_t v = right_remain_nodes[i];
				if (!g.right_delete[v]) {
					stop = false;
					right_remain_nodes_tmp[right_remain_nodes_tmp_num] = v;
					right_remain_nodes_tmp_num++;
					if (g.degree_v2[v] < right_k) {
						right_vertices_to_be_peeled.push_back(v);
					}
				}
			}
			swap(right_remain_nodes, right_remain_nodes_tmp);
			right_remain_nodes_num = right_remain_nodes_tmp_num;
			if (stop) break;
			while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
				// peel left
				int oo_ = left_vertices_to_be_peeled.size();
				for (int j = 0; j < oo_; j++) {
					vid_t u = left_vertices_to_be_peeled[j];
					if (g.left_delete[u]) continue;
					vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						vid_t v = tmp_neigh_[k];
						if (g.right_delete[v]) continue;
						dd_ = --g.degree_v2[v];
						if (update_flag && dd_ == 0) {
							dyn_crossUpdate_addition(g, alpha, pre_, v);
							g.right_delete[v] = true;
						}
						if (dd_ == pre_) {
							right_vertices_to_be_peeled.push_back(v);
						}
					}
					g.degree_v1[u] = 0;
					g.left_delete[u] = true;
					if (update_flag) {
						g.left_index[u][alpha] = pre_;
					}
				}
				left_vertices_to_be_peeled.clear();
				// peel right
				oo_ = right_vertices_to_be_peeled.size();
				for (int j = 0; j < oo_; j++) {
					vid_t v = right_vertices_to_be_peeled[j];
					if (g.right_delete[v]) continue;
					vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						vid_t u = tmp_neigh_[k];
						if (g.left_delete[u]) continue;
						dd_ = --g.degree_v1[u];
						if (update_flag && dd_ == 0) {
							g.left_index[u][alpha] = pre_;
							g.left_delete[u] = true;
						}
						if (dd_ == pre_left_k_) {
							left_vertices_to_be_peeled.push_back(u);
						}
					}
					g.degree_v2[v] = 0;
					g.right_delete[v] = true;
					if (update_flag) {
						dyn_crossUpdate_addition(g, alpha, pre_, v);
					}
				}
				right_vertices_to_be_peeled.clear();
			}
		}
		// restore graph
		int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
		fill_n(g.left_delete.begin(), g.left_delete.size(), false);
		for (vid_t u = 0; u < g.num_v1; u++) {
			g.degree_v1[u] = g.neighbor_v1[u].size();
			left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
		}
		g.v1_max_degree = left_degree_max;
		fill_n(g.right_delete.begin(), g.right_delete.size(), false);
		for (vid_t v = 0; v < g.num_v2; v++) {
			g.degree_v2[v] = g.neighbor_v2[v].size();
			right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
		}
		g.v2_max_degree = right_degree_max;
	}

	void batch_update_index_with_fixed_alpha_deletion(BiGraph& g, int alpha, int pi_alpha) {
		int dd_;
		int pre_left_k_ = alpha - 1;
		vector<vid_t> left_vertices_to_be_peeled;
		vector<vid_t> right_vertices_to_be_peeled;
		for (vid_t u = 0; u < g.getV1Num(); u++) {
			if (g.degree_v1[u] < alpha && !g.left_delete[u]) {
				left_vertices_to_be_peeled.push_back(u);
			}
		}
		int right_remain_nodes_num = g.num_v2;
		vector<vid_t> right_remain_nodes; right_remain_nodes.resize(g.num_v2);
		for (int i = 0; i < right_remain_nodes.size(); i++) {
			right_remain_nodes[i] = i;
		}
		int right_remain_nodes_tmp_num = 0;
		vector<vid_t> right_remain_nodes_tmp; right_remain_nodes_tmp.resize(g.num_v2);
		bool update_flag = false;
		for (int right_k = 1; right_k <= pi_alpha + 1; right_k++) {
			if (right_k - 1 > 0) {
				update_flag = true;
			}
			int pre_ = right_k - 1;
			bool stop = true;
			right_remain_nodes_tmp_num = 0;
			for (int i = 0; i < right_remain_nodes_num; i++) {
				vid_t v = right_remain_nodes[i];
				if (!g.right_delete[v]) {
					stop = false;
					right_remain_nodes_tmp[right_remain_nodes_tmp_num] = v;
					right_remain_nodes_tmp_num++;
					if (g.degree_v2[v] < right_k) {
						right_vertices_to_be_peeled.push_back(v);
					}
				}
			}
			swap(right_remain_nodes, right_remain_nodes_tmp);
			right_remain_nodes_num = right_remain_nodes_tmp_num;
			if (stop) break;
			while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
				// peel left
				int oo_ = left_vertices_to_be_peeled.size();
				for (int j = 0; j < oo_; j++) {
					vid_t u = left_vertices_to_be_peeled[j];
					if (g.left_delete[u]) continue;
					vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						vid_t v = tmp_neigh_[k];
						if (g.right_delete[v]) continue;
						dd_ = --g.degree_v2[v];
						if (update_flag && dd_ == 0) {
							dyn_crossUpdate_deletion(g, alpha, right_k, v);
							g.right_delete[v] = true;
						}
						if (dd_ == pre_) {
							right_vertices_to_be_peeled.push_back(v);
						}
					}
					g.degree_v1[u] = 0;
					g.left_delete[u] = true;
					if (update_flag) {
						g.left_index[u][alpha] = pre_;
					}
				}
				left_vertices_to_be_peeled.clear();
				// peel right
				oo_ = right_vertices_to_be_peeled.size();
				for (int j = 0; j < oo_; j++) {
					vid_t v = right_vertices_to_be_peeled[j];
					if (g.right_delete[v]) continue;
					vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						vid_t u = tmp_neigh_[k];
						if (g.left_delete[u]) continue;
						dd_ = --g.degree_v1[u];
						if (update_flag && dd_ == 0) {
							g.left_index[u][alpha] = pre_;
							g.left_delete[u] = true;
						}
						if (dd_ == pre_left_k_) {
							left_vertices_to_be_peeled.push_back(u);
						}
					}
					g.degree_v2[v] = 0;
					g.right_delete[v] = true;
					if (update_flag) {
						dyn_crossUpdate_deletion(g, alpha, right_k, v);
					}
				}
				right_vertices_to_be_peeled.clear();
			}
		}
		// restore graph
		int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
		fill_n(g.left_delete.begin(), g.left_delete.size(), false);
		for (vid_t u = 0; u < g.num_v1; u++) {
			g.degree_v1[u] = g.neighbor_v1[u].size();
			left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
		}
		g.v1_max_degree = left_degree_max;
		fill_n(g.right_delete.begin(), g.right_delete.size(), false);
		for (vid_t v = 0; v < g.num_v2; v++) {
			g.degree_v2[v] = g.neighbor_v2[v].size();
			right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
		}
		g.v2_max_degree = right_degree_max;
	}

	double update_bicore_index_batch(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v, 
		batch_set I, batch_set R) {
		auto start = chrono::system_clock::now();
		for (auto e : I) {
			vid_t u = e.first; vid_t v = e.second;
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						cout << "illegal insertion 1" << endl;
						return 0;
					}
				}
				g.addEdge(u, v);
				g.left_index[u].push_back(0);
				g.right_index[v].push_back(0);
			}
			else {
				cout << "illegal insertion 2" << endl;
				return 0;
			}
		}
		int delta = maxkcorenumber_test_optimal(g);
		vector<int> pi; pi.resize(delta + 1); fill_n(pi.begin(), pi.size(), 1e9);
		for (auto e : I) {
			vid_t u = e.first; vid_t v = e.second;
			for (int alpha = 1; alpha <= delta; alpha++) {
				if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
				if (g.left_index[u].size() > alpha) {
					pi[alpha] = pi[alpha] <= g.left_index[u][alpha] ? pi[alpha] : g.left_index[u][alpha];
				}
				int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
				pi[alpha] = pi[alpha] <= vbeta ? pi[alpha] : vbeta;
			}
		}
		
		for (int alpha = 1; alpha <= delta; alpha++) {
			batch_update_index_with_fixed_alpha_addition(g, alpha, pi[alpha]);
		}

		inv(g);
		fill_n(pi.begin(), pi.size(), 1e9);
		for (auto e : I) {
			vid_t v = e.first; vid_t u = e.second;
			for (int alpha = 1; alpha <= delta; alpha++) {
				if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
				if (g.left_index[u].size() > alpha) {
					pi[alpha] = pi[alpha] <= g.left_index[u][alpha] ? pi[alpha] : g.left_index[u][alpha];
				}
				int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
				pi[alpha] = pi[alpha] <= vbeta ? pi[alpha] : vbeta;
			}
		}

		for (int alpha = 1; alpha <= delta; alpha++) {
			batch_update_index_with_fixed_alpha_addition(g, alpha, pi[alpha]);
		}
		inv(g);

		// remove edge set
		vector<int> pi_1; pi_1.resize(delta + 1); fill_n(pi_1.begin(), pi_1.size(), 0);
		for (auto e : R) {
			vid_t u = e.first; vid_t v = e.second;
			for (int alpha = 1; alpha <= delta; alpha++) {
				if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
				if (g.left_index[u].size() > alpha) {
					pi_1[alpha] = pi_1[alpha] >= g.left_index[u][alpha] ? pi_1[alpha] : g.left_index[u][alpha];
				}
				int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
				pi_1[alpha] = pi_1[alpha] >= vbeta ? pi_1[alpha] : vbeta;
			}
		}
		inv(g);
		vector<int> pi_2; pi_2.resize(delta + 1); fill_n(pi_2.begin(), pi_2.size(), 0);
		for (auto e : R) {
			vid_t v = e.first; vid_t u = e.second;
			for (int alpha = 1; alpha <= delta; alpha++) {
				if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
				if (g.left_index[u].size() > alpha) {
					pi_2[alpha] = pi_2[alpha] >= g.left_index[u][alpha] ? pi_2[alpha] : g.left_index[u][alpha];
				}
				int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
				pi_2[alpha] = pi_2[alpha] >= vbeta ? pi_2[alpha] : vbeta;
			}
		}
		inv(g);

		for (auto e : R) {
			vid_t u = e.first; vid_t v = e.second;
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				bool deleted = false;
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						deleted = true;
						g.deleteEdge(u, v);
						g.left_index[u].pop_back();
						g.right_index[v].pop_back();
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
		delta = maxkcorenumber_test_optimal(g) + 1;

		for (int alpha = 1; alpha <= delta; alpha++) {
			batch_update_index_with_fixed_alpha_deletion(g, alpha, pi_1[alpha]);
		}

		inv(g);
		for (int alpha = 1; alpha <= delta; alpha++) {
			batch_update_index_with_fixed_alpha_deletion(g, alpha, pi_2[alpha]);
		}
		inv(g);

		auto end = chrono::system_clock::now();
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		chrono::duration<double> elapsed_seconds = end - start;
		return elapsed_seconds.count();
	}

	//////////////////////////////
	// parallel index maintenaince
	//////////////////////////////
	void multithread_insertion_limit_swap_A(int start, int end, vid_t u, vid_t v, vector<int>& tau_alpha, BiGraph& g) {
		for (int alpha = start; alpha <= end; alpha++) {
			update_index_with_fixed_left_k_addition_with_limit_swap(g, alpha, tau_alpha[alpha], u, v);
			// restore graph
			int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
			fill_n(g.left_delete.begin(), g.left_delete.size(), false);
			for (vid_t u = 0; u < g.num_v1; u++) {
				g.degree_v1[u] = g.neighbor_v1[u].size();
				left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
			}
			g.v1_max_degree = left_degree_max;
			fill_n(g.right_delete.begin(), g.right_delete.size(), false);
			for (vid_t v = 0; v < g.num_v2; v++) {
				g.degree_v2[v] = g.neighbor_v2[v].size();
				right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
			}
			g.v2_max_degree = right_degree_max;
		}
	}

	void multithread_insertion_limit_swap_B(int start, int end, vid_t u, vid_t v, vector<int>& tau_beta, BiGraph& g) {
		inv(g);
		for (int alpha = start; alpha <= end; alpha++) {
			update_index_with_fixed_left_k_addition_with_limit_swap(g, alpha, tau_beta[alpha], v, u);
			// restore graph
			int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
			fill_n(g.left_delete.begin(), g.left_delete.size(), false);
			for (vid_t u = 0; u < g.num_v1; u++) {
				g.degree_v1[u] = g.neighbor_v1[u].size();
				left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
			}
			g.v1_max_degree = left_degree_max;
			fill_n(g.right_delete.begin(), g.right_delete.size(), false);
			for (vid_t v = 0; v < g.num_v2; v++) {
				g.degree_v2[v] = g.neighbor_v2[v].size();
				right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
			}
			g.v2_max_degree = right_degree_max;
		}
		inv(g);
	}

	//void merge_insertion_result(BiGraph& target_g, vector<BiGraph>& bigraphsA, vector<BiGraph>& bigraphsB) {
	//	for (int i = 0; i < bigraphsA.size(); i++) {
	//		BiGraph& g_tmp = bigraphsA[i];
	//		for (vid_t u = 0; u < target_g.num_v1; u++) {
	//			for (int alpha = 1; alpha < target_g.left_index[u].size(); alpha++) {
	//				if (target_g.left_index[u][alpha] < g_tmp.left_index[u][alpha]) {
	//					target_g.left_index[u][alpha] = g_tmp.left_index[u][alpha];
	//				}
	//			}
	//		}
	//		for (vid_t v = 0; v < target_g.num_v2; v++) {
	//			for (int beta = 1; beta < target_g.right_index[v].size(); beta++) {
	//				if (target_g.right_index[v][beta] < g_tmp.right_index[v][beta]) {
	//					target_g.right_index[v][beta] = g_tmp.right_index[v][beta];
	//				}
	//			}
	//		}
	//	}
	//	for (int i = 0; i < bigraphsB.size(); i++) {
	//		BiGraph& g_tmp = bigraphsB[i];
	//		for (vid_t u = 0; u < target_g.num_v1; u++) {
	//			for (int alpha = 1; alpha < target_g.left_index[u].size(); alpha++) {
	//				if (target_g.left_index[u][alpha] < g_tmp.left_index[u][alpha]) {
	//					target_g.left_index[u][alpha] = g_tmp.left_index[u][alpha];
	//				}
	//			}
	//		}
	//		for (vid_t v = 0; v < target_g.num_v2; v++) {
	//			for (int beta = 1; beta < target_g.right_index[v].size(); beta++) {
	//				if (target_g.right_index[v][beta] < g_tmp.right_index[v][beta]) {
	//					target_g.right_index[v][beta] = g_tmp.right_index[v][beta];
	//				}
	//			}
	//		}
	//	}
	//}

	void merge_insertion_result(BiGraph& target_g, vector<BiGraph>& bigraphs) {
		for (int i = 0; i < bigraphs.size(); i++) {
			BiGraph& g_tmp = bigraphs[i];
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

	void multithread_removal_limit_swap_A(int start, int end, vid_t u, vid_t v, vector<int>& tau_alpha, vector<int>& alpha_start_bound, BiGraph& g) {
		for (int alpha = start; alpha <= end; alpha++) {
			update_index_with_fixed_left_k_deletion_with_limit_swap(g, alpha, tau_alpha[alpha], alpha_start_bound[alpha], u, v);
			// restore graph
			int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
			fill_n(g.left_delete.begin(), g.left_delete.size(), false);
			for (vid_t u = 0; u < g.num_v1; u++) {
				g.degree_v1[u] = g.neighbor_v1[u].size();
				left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
			}
			g.v1_max_degree = left_degree_max;
			fill_n(g.right_delete.begin(), g.right_delete.size(), false);
			for (vid_t v = 0; v < g.num_v2; v++) {
				g.degree_v2[v] = g.neighbor_v2[v].size();
				right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
			}
			g.v2_max_degree = right_degree_max;
		}
	}

	void multithread_removal_limit_swap_B(int start, int end, vid_t u, vid_t v, vector<int>& tau_beta, vector<int> beta_start_bound, BiGraph& g) {
		inv(g);
		for (int alpha = start; alpha <= end; alpha++) {
			update_index_with_fixed_left_k_deletion_with_limit_swap(g, alpha, tau_beta[alpha], beta_start_bound[alpha], v, u);
			// restore graph
			int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
			fill_n(g.left_delete.begin(), g.left_delete.size(), false);
			for (vid_t u = 0; u < g.num_v1; u++) {
				g.degree_v1[u] = g.neighbor_v1[u].size();
				left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
			}
			g.v1_max_degree = left_degree_max;
			fill_n(g.right_delete.begin(), g.right_delete.size(), false);
			for (vid_t v = 0; v < g.num_v2; v++) {
				g.degree_v2[v] = g.neighbor_v2[v].size();
				right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
			}
			g.v2_max_degree = right_degree_max;
		}
		inv(g);
	}

	void merge_removal_result(BiGraph& target_g, vector<BiGraph>& bigraphs) {
		for (int i = 0; i < bigraphs.size(); i++) {
			BiGraph& g_tmp = bigraphs[i];
			for (vid_t u = 0; u < target_g.num_v1; u++) {
				for (int alpha = 1; alpha < target_g.left_index[u].size(); alpha++) {
					if (target_g.left_index[u][alpha] > g_tmp.left_index[u][alpha]) {
						target_g.left_index[u][alpha] = g_tmp.left_index[u][alpha];
					}
				}
			}
			for (vid_t v = 0; v < target_g.num_v2; v++) {
				for (int beta = 1; beta < target_g.right_index[v].size(); beta++) {
					if (target_g.right_index[v][beta] > g_tmp.right_index[v][beta]) {
						target_g.right_index[v][beta] = g_tmp.right_index[v][beta];
					}
				}
			}
		}
	}

	void multithread_insertion_limit_swap_dfs_A(int start, int end, vid_t u, vid_t v, vector<int>& tau_alpha, BiGraph& g) {
		for (int alpha = start; alpha <= end; alpha++) {
			update_index_with_fixed_left_k_addition_swap_dfs(g, alpha, tau_alpha[alpha], u, v);
		}
	}

	void multithread_insertion_limit_swap_dfs_B(int start, int end, vid_t u, vid_t v, vector<int>& tau_beta, BiGraph& g) {
		inv(g);
		for (int alpha = start; alpha <= end; alpha++) {
			update_index_with_fixed_left_k_addition_swap_dfs(g, alpha, tau_beta[alpha], v, u);
		}
		inv(g);
	}

	void multithread_removal_limit_swap_dfs_A(int start, int end, vid_t u, vid_t v, vector<int>& tau_alpha, vector<int>& alpha_start_bound, BiGraph& g) {
		for (int alpha = start; alpha <= end; alpha++) {
			update_index_with_fixed_left_k_deletion_swap_dfs(g, alpha, tau_alpha[alpha], alpha_start_bound[alpha], u, v);
		}
	}

	void multithread_removal_limit_swap_dfs_B(int start, int end, vid_t u, vid_t v, vector<int>& tau_beta, vector<int> beta_start_bound, BiGraph& g) {
		inv(g);
		for (int alpha = start; alpha <= end; alpha++) {
			update_index_with_fixed_left_k_deletion_swap_dfs(g, alpha, tau_beta[alpha], beta_start_bound[alpha], v, u);
		}
		inv(g);
	}


	// return running time
	//double update_bicore_index_with_limit_swap_parallel(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	//	vid_t u, vid_t v, bool addition, int threads) {
	//	vector<BiGraph> bigraphsA;
	//	vector<BiGraph> bigraphsB;
	//	vector<thread> threadgroupsA;
	//	vector<thread> threadgroupsB;
	//	auto start = chrono::system_clock::now();
	//	auto end = chrono::system_clock::now();
	//	// calculate swap threshold
	//	pair<int, int> r = calculate_Delta(g);
	//	// check whether operation is legal and insert/remove the edge
	//	if (addition) {
	//		if (u < g.num_v1 && v < g.num_v2) {
	//			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
	//			int ss = tmp_neigh_.size();
	//			for (int i = 0; i < ss; i++) {
	//				if (tmp_neigh_[i] == v) {
	//					cout << "illegal insertion 1" << endl;
	//					return 0;
	//				}
	//			}
	//			g.addEdge(u, v);
	//		}
	//		else {
	//			cout << "illegal insertion 2" << endl;
	//			return 0;
	//		}
	//	}
	//	else {
	//		if (u < g.num_v1 && v < g.num_v2) {
	//			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
	//			int ss = tmp_neigh_.size();
	//			bool deleted = false;
	//			for (int i = 0; i < ss; i++) {
	//				if (tmp_neigh_[i] == v) {
	//					deleted = true;
	//					g.deleteEdge(u, v);
	//					break;
	//				}
	//			}
	//			if (!deleted) {
	//				cout << "illegal deletion 1" << endl;
	//				return 0;
	//			}
	//		}
	//		else {
	//			cout << "illegal deletion 2" << endl;
	//			return 0;
	//		}
	//	}
	//	// addition
	//	if (addition) {
	//		// extend core index and bicore index and set threshold
	//		g.left_index[u].push_back(0);
	//		g.right_index[v].push_back(0);
	//		int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
	//		int th_alpha, th_beta;
	//		if (r.first + 1 + r.second + 1 >= max_alpha) {
	//			th_alpha = max_alpha; th_beta = 0;
	//		}
	//		else if (r.first + 1 + r.second + 1 >= max_beta) {
	//			th_alpha = 0; th_beta = max_beta;
	//		}
	//		else {
	//			th_alpha = r.first + 1; th_beta = r.second + 1;
	//		}
	//		// compute alpha beta condition in advance
	//		vector<int> tau_alpha; vector<int> tau_beta;
	//		tau_alpha.push_back(-1); tau_beta.push_back(-1);
	//		for (int alpha = 1; alpha <= th_alpha; alpha++) {
	//			int ubeta = g.left_index[u][alpha];
	//			int bound = calculate_bound_for_left_node(g, u, alpha);
	//			int vbeta = -1;
	//			bool changed = false;
	//			for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
	//				if (g.right_index[v][b] < alpha) {
	//					vbeta = b - 1;
	//					changed = true;
	//					break;
	//				}
	//			}
	//			if (!changed) {
	//				vbeta = g.degree_v2[v];
	//			}
	//			if (vbeta < 0) {
	//				cout << "error: vbeta" << endl;
	//			}
	//			int beta = bound > vbeta ? vbeta : bound;
	//			tau_alpha.push_back(beta);
	//		}
	//		for (int beta = 1; beta <= th_beta; beta++) {
	//			int valpha = g.right_index[v][beta];
	//			int bound = calculate_bound_for_right_node(g, v, beta);
	//			int ualpha = -1;
	//			bool changed = false;
	//			for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
	//				if (g.left_index[u][a] < beta) {
	//					ualpha = a - 1;
	//					changed = true;
	//					break;
	//				}
	//			}
	//			if (!changed) {
	//				ualpha = g.degree_v1[u];
	//			}
	//			if (ualpha < 0) {
	//				cout << "error: ualpha" << endl;
	//			}
	//			int alpha = bound > ualpha ? ualpha : bound;
	//			tau_beta.push_back(alpha);
	//		}
	//		// allocate workload
	//		if (th_beta == 0) {
	//			int gap = ceil(th_alpha * 1.0 / (threads));
	//			for (int i = 0; i < threads; i++) {
	//				bigraphsA.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads; i++) {
	//				threadgroupsA.push_back(thread(multithread_insertion_limit_swap_A, i * gap + 1, MIN(th_alpha, i * gap + gap), u, v, ref(tau_alpha), ref(bigraphsA[i])));
	//			}
	//			for (int i = 0; i < threadgroupsA.size(); i++) {
	//				threadgroupsA[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		else if (th_alpha == 0) {
	//			int gap = ceil(th_beta * 1.0 / (threads));
	//			for (int i = 0; i < threads; i++) {
	//				bigraphsB.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads; i++) {
	//				threadgroupsB.push_back(thread(multithread_insertion_limit_swap_B, i * gap + 1, MIN(th_beta, i * gap + gap), u, v, ref(tau_beta), ref(bigraphsB[i])));
	//			}
	//			for (int i = 0; i < threadgroupsB.size(); i++) {
	//				threadgroupsB[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		else {
	//			int gap = ceil(MAX(th_alpha, th_beta) * 1.0 / (threads / 2));
	//			for (int i = 0; i < threads / 2; i++) {
	//				bigraphsA.push_back(g);
	//				bigraphsB.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads / 2; i++) {
	//				threadgroupsA.push_back(thread(multithread_insertion_limit_swap_A, i * gap + 1, MIN(th_alpha, i * gap + gap), u, v, ref(tau_alpha), ref(bigraphsA[i])));
	//			}
	//			for (int i = 0; i < threads / 2; i++) {
	//				threadgroupsB.push_back(thread(multithread_insertion_limit_swap_B, i * gap + 1, MIN(th_beta, i * gap + gap), u, v, ref(tau_beta), ref(bigraphsB[i])));
	//			}
	//			for (int i = 0; i < threadgroupsA.size(); i++) {
	//				threadgroupsA[i].join();
	//			}
	//			for (int i = 0; i < threadgroupsB.size(); i++) {
	//				threadgroupsB[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		merge_insertion_result(g, bigraphsA);
	//	}
	//	// deletion
	//	else {
	//		int oldusat = g.left_index[u].back();
	//		int oldvsat = g.right_index[v].back();
	//		int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
	//		int th_alpha, th_beta;
	//		if (r.first + r.second >= max_alpha) {
	//			th_alpha = max_alpha; th_beta = 0;
	//		}
	//		else if (r.first + r.second >= max_beta) {
	//			th_alpha = 0; th_beta = max_beta;
	//		}
	//		else {
	//			th_alpha = r.first; th_beta = r.second;
	//		}
	//		// compute alpha beta condition in advance
	//		vector<int> alpha_start; vector<int> beta_start; vector<int> alpha_start_bound; vector<int> beta_start_bound;
	//		alpha_start.push_back(-1); beta_start.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
	//		for (int alpha = 1; alpha <= th_alpha; alpha++) {
	//			int ubeta = g.left_index[u][alpha];
	//			int bound = calculate_bound_for_left_node(g, u, alpha);
	//			int vbeta = -1;
	//			bool changed = false;
	//			// g.degree_v2[v] has already been decreased
	//			for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
	//				if (g.right_index[v][b] < alpha) {
	//					vbeta = b - 1;
	//					changed = true;
	//					break;
	//				}
	//			}
	//			if (!changed) {
	//				vbeta = g.degree_v2[v] + 1;
	//			}
	//			if (vbeta == -1 || vbeta == 0) {
	//				cout << "error: vbeta" << endl;
	//			}
	//			int beta = ubeta > vbeta ? vbeta : ubeta;
	//			alpha_start.push_back(beta);
	//			alpha_start_bound.push_back(bound);
	//		}
	//		for (int beta = 1; beta <= th_beta; beta++) {
	//			int valpha = g.right_index[v][beta];
	//			int bound = calculate_bound_for_right_node(g, v, beta);
	//			int ualpha = -1;
	//			bool changed = false;
	//			for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
	//				if (g.left_index[u][a] < beta) {
	//					ualpha = a - 1;
	//					changed = true;
	//					break;
	//				}
	//			}
	//			if (!changed) {
	//				ualpha = g.degree_v1[u] + 1;
	//			}
	//			if (ualpha == -1 || ualpha == 0) {
	//				cout << "error: ualpha" << endl;
	//			}
	//			int alpha = valpha > ualpha ? ualpha : valpha;
	//			beta_start.push_back(alpha);
	//			beta_start_bound.push_back(bound);
	//		}
	//		// allocate workload
	//		if (th_beta == 0) {
	//			int gap = ceil(th_alpha * 1.0 / (threads));
	//			for (int i = 0; i < threads; i++) {
	//				bigraphsA.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads; i++) {
	//				threadgroupsA.push_back(thread(multithread_removal_limit_swap_A, i * gap + 1, MIN(th_alpha, i * gap + gap), u, v, ref(alpha_start), ref(alpha_start_bound), ref(bigraphsA[i])));
	//			}
	//			for (int i = 0; i < threadgroupsA.size(); i++) {
	//				threadgroupsA[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		else if (th_alpha == 0) {
	//			int gap = ceil(th_beta * 1.0 / (threads));
	//			for (int i = 0; i < threads; i++) {
	//				bigraphsB.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads; i++) {
	//				threadgroupsB.push_back(thread(multithread_removal_limit_swap_B, i * gap + 1, MIN(th_beta, i * gap + gap), u, v, ref(beta_start), ref(beta_start_bound), ref(bigraphsB[i])));
	//			}
	//			for (int i = 0; i < threadgroupsB.size(); i++) {
	//				threadgroupsB[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		else {
	//			int gap = ceil(MAX(th_alpha, th_beta) * 1.0 / (threads / 2));
	//			for (int i = 0; i < threads / 2; i++) {
	//				bigraphsA.push_back(g);
	//				bigraphsB.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads / 2; i++) {
	//				threadgroupsA.push_back(thread(multithread_removal_limit_swap_A, i * gap + 1, MIN(th_alpha, i * gap + gap), u, v, ref(alpha_start), ref(alpha_start_bound), ref(bigraphsA[i])));
	//			}
	//			for (int i = 0; i < threads / 2; i++) {
	//				threadgroupsB.push_back(thread(multithread_removal_limit_swap_B, i * gap + 1, MIN(th_beta, i * gap + gap), u, v, ref(beta_start), ref(beta_start_bound), ref(bigraphsB[i])));
	//			}
	//			for (int i = 0; i < threadgroupsA.size(); i++) {
	//				threadgroupsA[i].join();
	//			}
	//			for (int i = 0; i < threadgroupsB.size(); i++) {
	//				threadgroupsB[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		merge_removal_result(g, bigraphsA);
	//		g.left_index[u].pop_back();
	//		g.right_index[v].pop_back();
	//	}
	//	//build_bicore_index(g, bicore_index_u, bicore_index_v);
	//	chrono::duration<double> elapsed_seconds = end - start;
	//	return elapsed_seconds.count();
	//}

	void multithread_insertion_limit_swap(int alpha_end, int beta_end, vid_t u, vid_t v, vector<int>& tau_alpha, vector<int>& tau_beta, BiGraph& g) {
		bool left_finish = false; bool right_finish = false;
		while (!left_finish) {
			g_mutex.lock();
			if (alpha_mutex >= alpha_end) {
				left_finish = true;
			}
			else {
				alpha_mutex++;
				//cout << "alpha_mutex : " << alpha_mutex << endl;
			}
			int alpha_locked = alpha_mutex;
			g_mutex.unlock();
			if (!left_finish) {
				//cout << "alpha_mutex : " << alpha_locked << endl;
				update_index_with_fixed_left_k_addition_with_limit_swap(g, alpha_locked, tau_alpha[alpha_locked], u, v);
				// restore graph
				int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
				fill_n(g.left_delete.begin(), g.left_delete.size(), false);
				for (vid_t u = 0; u < g.num_v1; u++) {
					g.degree_v1[u] = g.neighbor_v1[u].size();
					left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
				}
				g.v1_max_degree = left_degree_max;
				fill_n(g.right_delete.begin(), g.right_delete.size(), false);
				for (vid_t v = 0; v < g.num_v2; v++) {
					g.degree_v2[v] = g.neighbor_v2[v].size();
					right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
				}
				g.v2_max_degree = right_degree_max;
			}
			//g_mutex.unlock();
		}
		inv(g);
		while (!right_finish) {
			g_mutex.lock();
			if (beta_mutex >= beta_end) {
				right_finish = true;
			}
			else {
				beta_mutex++;
				//cout << "beta_mutex : " << beta_mutex << endl;
			}
			int beta_locked = beta_mutex;
			g_mutex.unlock();
			if (!right_finish) {
				update_index_with_fixed_left_k_addition_with_limit_swap(g, beta_locked, tau_beta[beta_locked], v, u);
				// restore graph
				int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
				fill_n(g.left_delete.begin(), g.left_delete.size(), false);
				for (vid_t u = 0; u < g.num_v1; u++) {
					g.degree_v1[u] = g.neighbor_v1[u].size();
					left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
				}
				g.v1_max_degree = left_degree_max;
				fill_n(g.right_delete.begin(), g.right_delete.size(), false);
				for (vid_t v = 0; v < g.num_v2; v++) {
					g.degree_v2[v] = g.neighbor_v2[v].size();
					right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
				}
				g.v2_max_degree = right_degree_max;
			}
			//g_mutex.unlock();
		}
		inv(g);
	}
	void multithread_removal_limit_swap(int alpha_end, int beta_end, vid_t u, vid_t v, vector<int>& tau_alpha, vector<int>& tau_beta, vector<int>& alpha_start_bound, vector<int>& beta_start_bound, BiGraph& g) {
		bool left_finish = false; bool right_finish = false;
		while (!left_finish) {
			g_mutex.lock();
			if (alpha_mutex >= alpha_end) {
				left_finish = true;
			}
			else {
				alpha_mutex++;
			}
			int alpha_locked = alpha_mutex;
			g_mutex.unlock();
			if (!left_finish) {
				update_index_with_fixed_left_k_deletion_with_limit_swap(g, alpha_locked, tau_alpha[alpha_locked], alpha_start_bound[alpha_locked], u, v);
				// restore graph
				int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
				fill_n(g.left_delete.begin(), g.left_delete.size(), false);
				for (vid_t u = 0; u < g.num_v1; u++) {
					g.degree_v1[u] = g.neighbor_v1[u].size();
					left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
				}
				g.v1_max_degree = left_degree_max;
				fill_n(g.right_delete.begin(), g.right_delete.size(), false);
				for (vid_t v = 0; v < g.num_v2; v++) {
					g.degree_v2[v] = g.neighbor_v2[v].size();
					right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
				}
				g.v2_max_degree = right_degree_max;
			}
		}
		inv(g);
		while (!right_finish) {
			g_mutex.lock();
			if (beta_mutex >= beta_end) {
				right_finish = true;
			}
			else {
				beta_mutex++;
			}
			int beta_locked = beta_mutex;
			g_mutex.unlock();
			if (!right_finish) {
				update_index_with_fixed_left_k_deletion_with_limit_swap(g, beta_locked, tau_beta[beta_locked], beta_start_bound[beta_locked], v, u);
				// restore graph
				int left_degree_max, right_degree_max; left_degree_max = -1; right_degree_max = -1;
				fill_n(g.left_delete.begin(), g.left_delete.size(), false);
				for (vid_t u = 0; u < g.num_v1; u++) {
					g.degree_v1[u] = g.neighbor_v1[u].size();
					left_degree_max = g.degree_v1[u] > left_degree_max ? g.degree_v1[u] : left_degree_max;
				}
				g.v1_max_degree = left_degree_max;
				fill_n(g.right_delete.begin(), g.right_delete.size(), false);
				for (vid_t v = 0; v < g.num_v2; v++) {
					g.degree_v2[v] = g.neighbor_v2[v].size();
					right_degree_max = g.degree_v2[v] > right_degree_max ? g.degree_v2[v] : right_degree_max;
				}
				g.v2_max_degree = right_degree_max;
			}
		}
		inv(g);
	}

	double update_bicore_index_with_limit_swap_parallel(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition, int threads) {
		vector<BiGraph> bigraphs;
		vector<thread> threadgroups;
		auto start = chrono::system_clock::now();
		auto end = chrono::system_clock::now();
		// calculate swap threshold
		pair<int, int> r = calculate_Delta(g);
		// check whether operation is legal and insert/remove the edge
		if (addition) {
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						cout << "illegal insertion 1" << endl;
						return 0;
					}
				}
				g.addEdge(u, v);
			}
			else {
				cout << "illegal insertion 2" << endl;
				return 0;
			}
		}
		else {
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
		// addition
		if (addition) {
			// extend core index and bicore index and set threshold
			g.left_index[u].push_back(0);
			g.right_index[v].push_back(0);
			int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
			int th_alpha, th_beta;
			if (r.first + 1 + r.second + 1 >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + 1 + r.second + 1 >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first + 1; th_beta = r.second + 1;
			}
			// compute alpha beta condition in advance
			vector<int> tau_alpha; vector<int> tau_beta;
			tau_alpha.push_back(-1); tau_beta.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v];
				}
				if (vbeta < 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = bound > vbeta ? vbeta : bound;
				tau_alpha.push_back(beta);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u];
				}
				if (ualpha < 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = bound > ualpha ? ualpha : bound;
				tau_beta.push_back(alpha);
			}
			// allocate workload
			//cout << "sss" << th_alpha << endl;
			alpha_mutex = 0;
			beta_mutex = 0;
			for (int i = 0; i < threads; i++) {
				bigraphs.push_back(g);
			}
			start = chrono::system_clock::now();
			for (int i = 0; i < threads; i++) {
				threadgroups.push_back(thread(multithread_insertion_limit_swap, th_alpha, th_beta, u, v, ref(tau_alpha), ref(tau_beta), ref(bigraphs[i])));
			}
			for (int i = 0; i < threadgroups.size(); i++) {
				threadgroups[i].join();
			}
			end = chrono::system_clock::now();
			merge_insertion_result(g, bigraphs);
		}
		// deletion
		else {
			int oldusat = g.left_index[u].back();
			int oldvsat = g.right_index[v].back();
			int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
			int th_alpha, th_beta;
			if (r.first + r.second >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + r.second >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first; th_beta = r.second;
			}
			// compute alpha beta condition in advance
			vector<int> alpha_start; vector<int> beta_start; vector<int> alpha_start_bound; vector<int> beta_start_bound;
			alpha_start.push_back(-1); beta_start.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				// g.degree_v2[v] has already been decreased
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v] + 1;
				}
				if (vbeta == -1 || vbeta == 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = ubeta > vbeta ? vbeta : ubeta;
				alpha_start.push_back(beta);
				alpha_start_bound.push_back(bound);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u] + 1;
				}
				if (ualpha == -1 || ualpha == 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = valpha > ualpha ? ualpha : valpha;
				beta_start.push_back(alpha);
				beta_start_bound.push_back(bound);
			}
			// allocate workload
			alpha_mutex = 0;
			beta_mutex = 0;
			for (int i = 0; i < threads; i++) {
				bigraphs.push_back(g);
			}
			start = chrono::system_clock::now();
			for (int i = 0; i < threads; i++) {
				threadgroups.push_back(thread(multithread_removal_limit_swap, th_alpha, th_beta, u, v, ref(alpha_start), ref(beta_start), ref(alpha_start_bound), ref(beta_start_bound), ref(bigraphs[i])));
			}
			for (int i = 0; i < threadgroups.size(); i++) {
				threadgroups[i].join();
			}
			end = chrono::system_clock::now();
			merge_removal_result(g, bigraphs);
			g.left_index[u].pop_back();
			g.right_index[v].pop_back();
		}
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
		chrono::duration<double> elapsed_seconds = end - start;
		return elapsed_seconds.count();
	}

	// return running time
	//double update_bicore_index_with_limit_swap_dfs_parallel(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	//	vid_t u, vid_t v, bool addition, int threads) {
	//	vector<BiGraph> bigraphsA;
	//	vector<BiGraph> bigraphsB;
	//	vector<thread> threadgroupsA;
	//	vector<thread> threadgroupsB;
	//	auto start = chrono::system_clock::now();
	//	auto end = chrono::system_clock::now();
	//	// calculate swap threshold
	//	pair<int, int> r = calculate_Delta(g);
	//	// check whether operation is legal and insert/remove the edge
	//	if (addition) {
	//		if (u < g.num_v1 && v < g.num_v2) {
	//			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
	//			int ss = tmp_neigh_.size();
	//			for (int i = 0; i < ss; i++) {
	//				if (tmp_neigh_[i] == v) {
	//					cout << "illegal insertion 1" << endl;
	//					return 0;
	//				}
	//			}
	//			g.addEdge(u, v);
	//		}
	//		else {
	//			cout << "illegal insertion 2" << endl;
	//			return 0;
	//		}
	//	}
	//	else {
	//		if (u < g.num_v1 && v < g.num_v2) {
	//			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
	//			int ss = tmp_neigh_.size();
	//			bool deleted = false;
	//			for (int i = 0; i < ss; i++) {
	//				if (tmp_neigh_[i] == v) {
	//					deleted = true;
	//					g.deleteEdge(u, v);
	//					break;
	//				}
	//			}
	//			if (!deleted) {
	//				cout << "illegal deletion 1" << endl;
	//				return 0;
	//			}
	//		}
	//		else {
	//			cout << "illegal deletion 2" << endl;
	//			return 0;
	//		}
	//	}
	//	// addition
	//	if (addition) {
	//		// extend core index and bicore index and set threshold
	//		g.left_index[u].push_back(0);
	//		g.right_index[v].push_back(0);
	//		int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
	//		int th_alpha, th_beta;
	//		if (r.first + 1 + r.second + 1 >= max_alpha) {
	//			th_alpha = max_alpha; th_beta = 0;
	//		}
	//		else if (r.first + 1 + r.second + 1 >= max_beta) {
	//			th_alpha = 0; th_beta = max_beta;
	//		}
	//		else {
	//			th_alpha = r.first + 1; th_beta = r.second + 1;
	//		}
	//		// compute alpha beta condition in advance
	//		vector<int> tau_alpha; vector<int> tau_beta;
	//		tau_alpha.push_back(-1); tau_beta.push_back(-1);
	//		for (int alpha = 1; alpha <= th_alpha; alpha++) {
	//			int ubeta = g.left_index[u][alpha];
	//			int bound = calculate_bound_for_left_node(g, u, alpha);
	//			int vbeta = -1;
	//			bool changed = false;
	//			for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
	//				if (g.right_index[v][b] < alpha) {
	//					vbeta = b - 1;
	//					changed = true;
	//					break;
	//				}
	//			}
	//			if (!changed) {
	//				vbeta = g.degree_v2[v];
	//			}
	//			if (vbeta < 0) {
	//				cout << "error: vbeta" << endl;
	//			}
	//			int beta = bound > vbeta ? vbeta : bound;
	//			tau_alpha.push_back(beta);
	//		}
	//		for (int beta = 1; beta <= th_beta; beta++) {
	//			int valpha = g.right_index[v][beta];
	//			int bound = calculate_bound_for_right_node(g, v, beta);
	//			int ualpha = -1;
	//			bool changed = false;
	//			for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
	//				if (g.left_index[u][a] < beta) {
	//					ualpha = a - 1;
	//					changed = true;
	//					break;
	//				}
	//			}
	//			if (!changed) {
	//				ualpha = g.degree_v1[u];
	//			}
	//			if (ualpha < 0) {
	//				cout << "error: ualpha" << endl;
	//			}
	//			int alpha = bound > ualpha ? ualpha : bound;
	//			tau_beta.push_back(alpha);
	//		}
	//		// allocate workload
	//		if (th_beta == 0) {
	//			int gap = ceil(th_alpha * 1.0 / (threads));
	//			for (int i = 0; i < threads; i++) {
	//				bigraphsA.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads; i++) {
	//				threadgroupsA.push_back(thread(multithread_insertion_limit_swap_dfs_A, i * gap + 1, MIN(th_alpha, i * gap + gap), u, v, ref(tau_alpha), ref(bigraphsA[i])));
	//			}
	//			for (int i = 0; i < threadgroupsA.size(); i++) {
	//				threadgroupsA[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		else if (th_alpha == 0) {
	//			int gap = ceil(th_beta * 1.0 / (threads));
	//			for (int i = 0; i < threads; i++) {
	//				bigraphsB.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads; i++) {
	//				threadgroupsB.push_back(thread(multithread_insertion_limit_swap_dfs_B, i * gap + 1, MIN(th_beta, i * gap + gap), u, v, ref(tau_beta), ref(bigraphsB[i])));
	//			}
	//			for (int i = 0; i < threadgroupsB.size(); i++) {
	//				threadgroupsB[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		else {
	//			int gap = ceil(MAX(th_alpha, th_beta) * 1.0 / (threads / 2));
	//			for (int i = 0; i < threads / 2; i++) {
	//				bigraphsA.push_back(g);
	//				bigraphsB.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads / 2; i++) {
	//				threadgroupsA.push_back(thread(multithread_insertion_limit_swap_dfs_A, i * gap + 1, MIN(th_alpha, i * gap + gap), u, v, ref(tau_alpha), ref(bigraphsA[i])));
	//			}
	//			for (int i = 0; i < threads / 2; i++) {
	//				threadgroupsB.push_back(thread(multithread_insertion_limit_swap_dfs_B, i * gap + 1, MIN(th_beta, i * gap + gap), u, v, ref(tau_beta), ref(bigraphsB[i])));
	//			}
	//			for (int i = 0; i < threadgroupsA.size(); i++) {
	//				threadgroupsA[i].join();
	//			}
	//			for (int i = 0; i < threadgroupsB.size(); i++) {
	//				threadgroupsB[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		merge_insertion_result(g, bigraphsA, bigraphsB);
	//	}
	//	// deletion
	//	else {
	//		int oldusat = g.left_index[u].back();
	//		int oldvsat = g.right_index[v].back();
	//		int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
	//		int th_alpha, th_beta;
	//		if (r.first + r.second >= max_alpha) {
	//			th_alpha = max_alpha; th_beta = 0;
	//		}
	//		else if (r.first + r.second >= max_beta) {
	//			th_alpha = 0; th_beta = max_beta;
	//		}
	//		else {
	//			th_alpha = r.first; th_beta = r.second;
	//		}
	//		// compute alpha beta condition in advance
	//		vector<int> alpha_start; vector<int> beta_start; vector<int> alpha_start_bound; vector<int> beta_start_bound;
	//		alpha_start.push_back(-1); beta_start.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
	//		for (int alpha = 1; alpha <= th_alpha; alpha++) {
	//			int ubeta = g.left_index[u][alpha];
	//			int bound = calculate_bound_for_left_node(g, u, alpha);
	//			int vbeta = -1;
	//			bool changed = false;
	//			// g.degree_v2[v] has already been decreased
	//			for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
	//				if (g.right_index[v][b] < alpha) {
	//					vbeta = b - 1;
	//					changed = true;
	//					break;
	//				}
	//			}
	//			if (!changed) {
	//				vbeta = g.degree_v2[v] + 1;
	//			}
	//			if (vbeta == -1 || vbeta == 0) {
	//				cout << "error: vbeta" << endl;
	//			}
	//			int beta = ubeta > vbeta ? vbeta : ubeta;
	//			alpha_start.push_back(beta);
	//			alpha_start_bound.push_back(bound);
	//		}
	//		for (int beta = 1; beta <= th_beta; beta++) {
	//			int valpha = g.right_index[v][beta];
	//			int bound = calculate_bound_for_right_node(g, v, beta);
	//			int ualpha = -1;
	//			bool changed = false;
	//			for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
	//				if (g.left_index[u][a] < beta) {
	//					ualpha = a - 1;
	//					changed = true;
	//					break;
	//				}
	//			}
	//			if (!changed) {
	//				ualpha = g.degree_v1[u] + 1;
	//			}
	//			if (ualpha == -1 || ualpha == 0) {
	//				cout << "error: ualpha" << endl;
	//			}
	//			int alpha = valpha > ualpha ? ualpha : valpha;
	//			beta_start.push_back(alpha);
	//			beta_start_bound.push_back(bound);
	//		}
	//		// allocate workload
	//		if (th_beta == 0) {
	//			int gap = ceil(th_alpha * 1.0 / (threads));
	//			for (int i = 0; i < threads; i++) {
	//				bigraphsA.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads; i++) {
	//				threadgroupsA.push_back(thread(multithread_removal_limit_swap_dfs_A, i * gap + 1, MIN(th_alpha, i * gap + gap), u, v, ref(alpha_start), ref(alpha_start_bound), ref(bigraphsA[i])));
	//			}
	//			for (int i = 0; i < threadgroupsA.size(); i++) {
	//				threadgroupsA[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		else if (th_alpha == 0) {
	//			int gap = ceil(th_beta * 1.0 / (threads));
	//			for (int i = 0; i < threads; i++) {
	//				bigraphsB.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads; i++) {
	//				threadgroupsB.push_back(thread(multithread_removal_limit_swap_dfs_B, i * gap + 1, MIN(th_beta, i * gap + gap), u, v, ref(beta_start), ref(beta_start_bound), ref(bigraphsB[i])));
	//			}
	//			for (int i = 0; i < threadgroupsB.size(); i++) {
	//				threadgroupsB[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		else {
	//			int gap = ceil(MAX(th_alpha, th_beta) * 1.0 / (threads / 2));
	//			for (int i = 0; i < threads / 2; i++) {
	//				bigraphsA.push_back(g);
	//				bigraphsB.push_back(g);
	//			}
	//			start = chrono::system_clock::now();
	//			for (int i = 0; i < threads / 2; i++) {
	//				threadgroupsA.push_back(thread(multithread_removal_limit_swap_dfs_A, i * gap + 1, MIN(th_alpha, i * gap + gap), u, v, ref(alpha_start), ref(alpha_start_bound), ref(bigraphsA[i])));
	//			}
	//			for (int i = 0; i < threads / 2; i++) {
	//				threadgroupsB.push_back(thread(multithread_removal_limit_swap_dfs_B, i * gap + 1, MIN(th_beta, i * gap + gap), u, v, ref(beta_start), ref(beta_start_bound), ref(bigraphsB[i])));
	//			}
	//			for (int i = 0; i < threadgroupsA.size(); i++) {
	//				threadgroupsA[i].join();
	//			}
	//			for (int i = 0; i < threadgroupsB.size(); i++) {
	//				threadgroupsB[i].join();
	//			}
	//			end = chrono::system_clock::now();
	//		}
	//		merge_removal_result(g, bigraphsA, bigraphsB);
	//		g.left_index[u].pop_back();
	//		g.right_index[v].pop_back();
	//	}
	//	//build_bicore_index(g, bicore_index_u, bicore_index_v);
	//	chrono::duration<double> elapsed_seconds = end - start;
	//	return elapsed_seconds.count();
	//}
	
	
	void multithread_insertion_limit_swap_dfs(int alpha_end, int beta_end, vid_t u, vid_t v, vector<int>& tau_alpha, vector<int>& tau_beta, BiGraph& g) {	
		bool left_finish = false; bool right_finish = false;
		while (!left_finish) {
			g_mutex.lock();
			if (alpha_mutex >= alpha_end) {
				left_finish = true;
			}
			else {
				alpha_mutex++;
				//cout << "alpha_mutex : " << alpha_mutex << endl;
			}
			int alpha_locked = alpha_mutex;
			g_mutex.unlock();
			if (!left_finish) {
				//cout << "alpha_mutex : " << alpha_locked << endl;
				update_index_with_fixed_left_k_addition_swap_dfs(g, alpha_locked, tau_alpha[alpha_locked], u, v);
			}
			//g_mutex.unlock();
		}
		inv(g);
		while (!right_finish) {
			g_mutex.lock();
			if (beta_mutex >= beta_end) {
				right_finish = true;
			}
			else {
				beta_mutex++;
				//cout << "beta_mutex : " << beta_mutex << endl;
			}
			int beta_locked = beta_mutex;
			g_mutex.unlock();
			if (!right_finish) {
				update_index_with_fixed_left_k_addition_swap_dfs(g, beta_locked, tau_beta[beta_locked], v, u);
			}
			//g_mutex.unlock();
		}
		inv(g);
	}
	void multithread_removal_limit_swap_dfs(int alpha_end, int beta_end, vid_t u, vid_t v, vector<int>& tau_alpha, vector<int>& tau_beta, vector<int>& alpha_start_bound, vector<int>& beta_start_bound, BiGraph& g) {
		bool left_finish = false; bool right_finish = false;
		while (!left_finish) {
			g_mutex.lock();
			if (alpha_mutex >= alpha_end) {
				left_finish = true;
			}
			else {
				alpha_mutex++;
			}
			int alpha_locked = alpha_mutex;
			g_mutex.unlock();
			if (!left_finish) {
				update_index_with_fixed_left_k_deletion_swap_dfs(g, alpha_locked, tau_alpha[alpha_locked], alpha_start_bound[alpha_locked], u, v);
			}
		}
		inv(g);
		while (!right_finish) {
			g_mutex.lock();
			if (beta_mutex >= beta_end) {
				right_finish = true;
			}
			else {
				beta_mutex++;
			}
			int beta_locked = beta_mutex;
			g_mutex.unlock();
			if (!right_finish) {
				update_index_with_fixed_left_k_deletion_swap_dfs(g, beta_locked, tau_beta[beta_locked], beta_start_bound[beta_locked], v, u);
			}
		}
		inv(g);
	}
	double update_bicore_index_with_limit_swap_dfs_parallel(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition, int threads) {
		vector<BiGraph> bigraphs;
		vector<thread> threadgroups;
		auto start = chrono::system_clock::now();
		auto end = chrono::system_clock::now();
		// calculate swap threshold
		pair<int, int> r = calculate_Delta(g);
		// check whether operation is legal and insert/remove the edge
		if (addition) {
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						cout << "illegal insertion 1" << endl;
						return 0;
					}
				}
				g.addEdge(u, v);
			}
			else {
				cout << "illegal insertion 2" << endl;
				return 0;
			}
		}
		else {
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
		// addition
		if (addition) {
			// extend core index and bicore index and set threshold
			g.left_index[u].push_back(0);
			g.right_index[v].push_back(0);
			int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
			int th_alpha, th_beta;
			if (r.first + 1 + r.second + 1 >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + 1 + r.second + 1 >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first + 1; th_beta = r.second + 1;
			}
			// compute alpha beta condition in advance
			vector<int> tau_alpha; vector<int> tau_beta;
			tau_alpha.push_back(-1); tau_beta.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v];
				}
				if (vbeta < 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = bound > vbeta ? vbeta : bound;
				tau_alpha.push_back(beta);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u];
				}
				if (ualpha < 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = bound > ualpha ? ualpha : bound;
				tau_beta.push_back(alpha);
			}
			// allocate workload
			//cout << "sss" << th_alpha << endl;
			alpha_mutex = 0;
			beta_mutex = 0;
			for (int i = 0; i < threads; i++) {
				bigraphs.push_back(g);
			}
			start = chrono::system_clock::now();
			for (int i = 0; i < threads; i++) {
				threadgroups.push_back(thread(multithread_insertion_limit_swap_dfs, th_alpha, th_beta, u, v, ref(tau_alpha), ref(tau_beta), ref(bigraphs[i])));
			}
			for (int i = 0; i < threadgroups.size(); i++) {
				threadgroups[i].join();
			}
			end = chrono::system_clock::now();
			merge_insertion_result(g, bigraphs);
		}
		// deletion
		else {
			int oldusat = g.left_index[u].back();
			int oldvsat = g.right_index[v].back();
			int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
			int th_alpha, th_beta;
			if (r.first + r.second >= max_alpha) {
				th_alpha = max_alpha; th_beta = 0;
			}
			else if (r.first + r.second >= max_beta) {
				th_alpha = 0; th_beta = max_beta;
			}
			else {
				th_alpha = r.first; th_beta = r.second;
			}
			// compute alpha beta condition in advance
			vector<int> alpha_start; vector<int> beta_start; vector<int> alpha_start_bound; vector<int> beta_start_bound;
			alpha_start.push_back(-1); beta_start.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
			for (int alpha = 1; alpha <= th_alpha; alpha++) {
				int ubeta = g.left_index[u][alpha];
				int bound = calculate_bound_for_left_node(g, u, alpha);
				int vbeta = -1;
				bool changed = false;
				// g.degree_v2[v] has already been decreased
				for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
					if (g.right_index[v][b] < alpha) {
						vbeta = b - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					vbeta = g.degree_v2[v] + 1;
				}
				if (vbeta == -1 || vbeta == 0) {
					cout << "error: vbeta" << endl;
				}
				int beta = ubeta > vbeta ? vbeta : ubeta;
				alpha_start.push_back(beta);
				alpha_start_bound.push_back(bound);
			}
			for (int beta = 1; beta <= th_beta; beta++) {
				int valpha = g.right_index[v][beta];
				int bound = calculate_bound_for_right_node(g, v, beta);
				int ualpha = -1;
				bool changed = false;
				for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
					if (g.left_index[u][a] < beta) {
						ualpha = a - 1;
						changed = true;
						break;
					}
				}
				if (!changed) {
					ualpha = g.degree_v1[u] + 1;
				}
				if (ualpha == -1 || ualpha == 0) {
					cout << "error: ualpha" << endl;
				}
				int alpha = valpha > ualpha ? ualpha : valpha;
				beta_start.push_back(alpha);
				beta_start_bound.push_back(bound);
			}
			// allocate workload
			alpha_mutex = 0;
			beta_mutex = 0;
			for (int i = 0; i < threads; i++) {
				bigraphs.push_back(g);
			}
			start = chrono::system_clock::now();
			for (int i = 0; i < threads; i++) {
				threadgroups.push_back(thread(multithread_removal_limit_swap_dfs, th_alpha, th_beta, u, v, ref(alpha_start), ref(beta_start), ref(alpha_start_bound), ref(beta_start_bound), ref(bigraphs[i])));
			}
			for (int i = 0; i < threadgroups.size(); i++) {
				threadgroups[i].join();
			}
			end = chrono::system_clock::now();
			merge_removal_result(g, bigraphs);
			g.left_index[u].pop_back();
			g.right_index[v].pop_back();
		}
		chrono::duration<double> elapsed_seconds = end - start;
		return elapsed_seconds.count();
	}

	void multithread_insertion_batch_A(int start, int end, vector<int>& pi_A, BiGraph& g) {
		for (int alpha = start; alpha <= end; alpha++) {
			batch_update_index_with_fixed_alpha_addition(g, alpha, pi_A[alpha]);
		}
	}

	void multithread_insertion_batch_B(int start, int end, vector<int>& pi_B, BiGraph& g) {
		inv(g);
		for (int alpha = start; alpha <= end; alpha++) {
			batch_update_index_with_fixed_alpha_addition(g, alpha, pi_B[alpha]);
		}
		inv(g);
	}

	void multithread_removal_batch_A(int start, int end, vector<int>& pi_A, BiGraph& g) {
		for (int alpha = start; alpha <= end; alpha++) {
			batch_update_index_with_fixed_alpha_deletion(g, alpha, pi_A[alpha]);
		}
	}

	void multithread_removal_batch_B(int start, int end, vector<int>& pi_B, BiGraph& g) {
		inv(g);
		for (int alpha = start; alpha <= end; alpha++) {
			batch_update_index_with_fixed_alpha_deletion(g, alpha, pi_B[alpha]);
		}
		inv(g);
	}

	void multithread_insertion_batch(int alpha_end, int beta_end, vector<int>& pi_A, vector<int>& pi_B, BiGraph& g) {
		bool left_finish = false; bool right_finish = false;
		while (!left_finish) {
			g_mutex.lock();
			if (alpha_mutex >= alpha_end) {
				left_finish = true;
			}
			else {
				alpha_mutex++;
				//cout << "alpha_mutex : " << alpha_mutex << endl;
			}
			int alpha_locked = alpha_mutex;
			g_mutex.unlock();
			if (!left_finish) {
				//cout << "alpha_mutex : " << alpha_locked << endl;
				batch_update_index_with_fixed_alpha_addition(g, alpha_locked, pi_A[alpha_locked]);
			}
			//g_mutex.unlock();
		}
		inv(g);
		while (!right_finish) {
			g_mutex.lock();
			if (beta_mutex >= beta_end) {
				right_finish = true;
			}
			else {
				beta_mutex++;
				//cout << "beta_mutex : " << beta_mutex << endl;
			}
			int beta_locked = beta_mutex;
			g_mutex.unlock();
			if (!right_finish) {
				batch_update_index_with_fixed_alpha_addition(g, beta_locked, pi_B[beta_locked]);
			}
			//g_mutex.unlock();
		}
		inv(g);
	}
	void multithread_removal_batch(int alpha_end, int beta_end, vector<int>& pi_A, vector<int>& pi_B, BiGraph& g) {
		bool left_finish = false; bool right_finish = false;
		while (!left_finish) {
			g_mutex.lock();
			if (alpha_mutex >= alpha_end) {
				left_finish = true;
			}
			else {
				alpha_mutex++;
			}
			int alpha_locked = alpha_mutex;
			g_mutex.unlock();
			if (!left_finish) {
				batch_update_index_with_fixed_alpha_deletion(g, alpha_locked, pi_A[alpha_locked]);
			}
		}
		inv(g);
		while (!right_finish) {
			g_mutex.lock();
			if (beta_mutex >= beta_end) {
				right_finish = true;
			}
			else {
				beta_mutex++;
			}
			int beta_locked = beta_mutex;
			g_mutex.unlock();
			if (!right_finish) {
				batch_update_index_with_fixed_alpha_deletion(g, beta_locked, pi_B[beta_locked]);
			}
		}
		inv(g);
	}

	double update_bicore_index_batch_parallel(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		batch_set I, batch_set R, int threads) {
		vector<BiGraph> bigraphs;
		vector<thread> threadgroups;
		auto start = chrono::system_clock::now();
		auto end = chrono::system_clock::now();
		for (auto e : I) {
			vid_t u = e.first; vid_t v = e.second;
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						cout << "illegal insertion 1" << endl;
						return 0;
					}
				}
				g.addEdge(u, v);
				g.left_index[u].push_back(0);
				g.right_index[v].push_back(0);
			}
			else {
				cout << "illegal insertion 2" << endl;
				return 0;
			}
		}
		int delta = maxkcorenumber_test_optimal(g);
		vector<int> pi_A; pi_A.resize(delta + 1); fill_n(pi_A.begin(), pi_A.size(), 1e9);
		for (auto e : I) {
			vid_t u = e.first; vid_t v = e.second;
			for (int alpha = 1; alpha <= delta; alpha++) {
				if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
				if (g.left_index[u].size() > alpha) {
					pi_A[alpha] = pi_A[alpha] <= g.left_index[u][alpha] ? pi_A[alpha] : g.left_index[u][alpha];
				}
				int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
				pi_A[alpha] = pi_A[alpha] <= vbeta ? pi_A[alpha] : vbeta;
			}
		}
		inv(g);
		vector<int> pi_B; pi_B.resize(delta + 1); fill_n(pi_B.begin(), pi_B.size(), 1e9);
		for (auto e : I) {
			vid_t v = e.first; vid_t u = e.second;
			for (int alpha = 1; alpha <= delta; alpha++) {
				if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
				if (g.left_index[u].size() > alpha) {
					pi_B[alpha] = pi_B[alpha] <= g.left_index[u][alpha] ? pi_B[alpha] : g.left_index[u][alpha];
				}
				int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
				pi_B[alpha] = pi_B[alpha] <= vbeta ? pi_B[alpha] : vbeta;
			}
		}
		inv(g);

		alpha_mutex = 0;
		beta_mutex = 0;
		for (int i = 0; i < threads; i++) {
			bigraphs.push_back(g);
		}
		start = chrono::system_clock::now();
		for (int i = 0; i < threads; i++) {
			threadgroups.push_back(thread(multithread_insertion_batch, delta, delta, ref(pi_A), ref(pi_B), ref(bigraphs[i])));
		}
		for (int i = 0; i < threadgroups.size(); i++) {
			threadgroups[i].join();
		}
		end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		double time = elapsed_seconds.count();
		merge_insertion_result(g, bigraphs);

		// remove part
		pi_A.resize(delta + 1); fill_n(pi_A.begin(), pi_A.size(), 0);
		for (auto e : R) {
			vid_t u = e.first; vid_t v = e.second;
			for (int alpha = 1; alpha <= delta; alpha++) {
				if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
				if (g.left_index[u].size() > alpha) {
					pi_A[alpha] = pi_A[alpha] >= g.left_index[u][alpha] ? pi_A[alpha] : g.left_index[u][alpha];
				}
				int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
				pi_A[alpha] = pi_A[alpha] >= vbeta ? pi_A[alpha] : vbeta;
			}
		}
		inv(g);
		pi_B.resize(delta + 1); fill_n(pi_B.begin(), pi_B.size(), 0);
		for (auto e : R) {
			vid_t v = e.first; vid_t u = e.second;
			for (int alpha = 1; alpha <= delta; alpha++) {
				if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
				if (g.left_index[u].size() > alpha) {
					pi_B[alpha] = pi_B[alpha] >= g.left_index[u][alpha] ? pi_B[alpha] : g.left_index[u][alpha];
				}
				int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
				pi_B[alpha] = pi_B[alpha] >= vbeta ? pi_B[alpha] : vbeta;
			}
		}
		inv(g);

		for (auto e : R) {
			vid_t u = e.first; vid_t v = e.second;
			if (u < g.num_v1 && v < g.num_v2) {
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				bool deleted = false;
				for (int i = 0; i < ss; i++) {
					if (tmp_neigh_[i] == v) {
						deleted = true;
						g.deleteEdge(u, v);
						g.left_index[u].pop_back();
						g.right_index[v].pop_back();
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
		delta = maxkcorenumber_test_optimal(g) + 1;
		bigraphs.clear();
		threadgroups.clear();
		alpha_mutex = 0;
		beta_mutex = 0;
		for (int i = 0; i < threads; i++) {
			bigraphs.push_back(g);
		}
		start = chrono::system_clock::now();
		for (int i = 0; i < threads; i++) {
			threadgroups.push_back(thread(multithread_removal_batch, delta, delta, ref(pi_A), ref(pi_B), ref(bigraphs[i])));
		}
		for (int i = 0; i < threadgroups.size(); i++) {
			threadgroups[i].join();
		}
		end = chrono::system_clock::now();
		elapsed_seconds = end - start;
		time += elapsed_seconds.count();
		merge_removal_result(g, bigraphs);
		return time;
	}

	//double update_bicore_index_batch_parallel(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	//	batch_set I, batch_set R, int threads) {
	//	vector<BiGraph> bigraphsA;
	//	vector<BiGraph> bigraphsB;
	//	vector<thread> threadgroupsA;
	//	vector<thread> threadgroupsB;
	//	auto start = chrono::system_clock::now();
	//	auto end = chrono::system_clock::now();
	//	for (auto e : I) {
	//		vid_t u = e.first; vid_t v = e.second;
	//		if (u < g.num_v1 && v < g.num_v2) {
	//			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
	//			int ss = tmp_neigh_.size();
	//			for (int i = 0; i < ss; i++) {
	//				if (tmp_neigh_[i] == v) {
	//					cout << "illegal insertion 1" << endl;
	//					return 0;
	//				}
	//			}
	//			g.addEdge(u, v);
	//			g.left_index[u].push_back(0);
	//			g.right_index[v].push_back(0);
	//		}
	//		else {
	//			cout << "illegal insertion 2" << endl;
	//			return 0;
	//		}
	//	}
	//	int delta = maxkcorenumber_test_optimal(g);
	//	vector<int> pi_A; pi_A.resize(delta + 1); fill_n(pi_A.begin(), pi_A.size(), 1e9);
	//	for (auto e : I) {
	//		vid_t u = e.first; vid_t v = e.second;
	//		for (int alpha = 1; alpha <= delta; alpha++) {
	//			if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
	//			if (g.left_index[u].size() > alpha) {
	//				pi_A[alpha] = pi_A[alpha] <= g.left_index[u][alpha] ? pi_A[alpha] : g.left_index[u][alpha];
	//			}
	//			int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
	//			pi_A[alpha] = pi_A[alpha] <= vbeta ? pi_A[alpha] : vbeta;
	//		}
	//	}
	//	inv(g);
	//	vector<int> pi_B; pi_B.resize(delta + 1); fill_n(pi_B.begin(), pi_B.size(), 1e9);
	//	for (auto e : I) {
	//		vid_t v = e.first; vid_t u = e.second;
	//		for (int alpha = 1; alpha <= delta; alpha++) {
	//			if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
	//			if (g.left_index[u].size() > alpha) {
	//				pi_B[alpha] = pi_B[alpha] <= g.left_index[u][alpha] ? pi_B[alpha] : g.left_index[u][alpha];
	//			}
	//			int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
	//			pi_B[alpha] = pi_B[alpha] <= vbeta ? pi_B[alpha] : vbeta;
	//		}
	//	}
	//	inv(g);

	//	int gap = ceil(delta * 1.0 / (threads / 2));
	//	for (int i = 0; i < threads / 2; i++) {
	//		bigraphsA.push_back(g);
	//		bigraphsB.push_back(g);
	//	}
	//	start = chrono::system_clock::now();
	//	for (int i = 0; i < threads / 2; i++) {
	//		threadgroupsA.push_back(thread(multithread_insertion_batch_A, i * gap + 1, MIN(delta, i * gap + gap), ref(pi_A), ref(bigraphsA[i])));
	//	}
	//	for (int i = 0; i < threads / 2; i++) {
	//		threadgroupsB.push_back(thread(multithread_insertion_batch_B, i * gap + 1, MIN(delta, i * gap + gap), ref(pi_B), ref(bigraphsB[i])));
	//	}
	//	for (int i = 0; i < threadgroupsA.size(); i++) {
	//		threadgroupsA[i].join();
	//	}
	//	for (int i = 0; i < threadgroupsB.size(); i++) {
	//		threadgroupsB[i].join();
	//	}
	//	end = chrono::system_clock::now();
	//	chrono::duration<double> elapsed_seconds = end - start;
	//	double time = elapsed_seconds.count();
	//	merge_insertion_result(g, bigraphsA);

	//	// remove part
	//	pi_A.resize(delta + 1); fill_n(pi_A.begin(), pi_A.size(), 0);
	//	for (auto e : R) {
	//		vid_t u = e.first; vid_t v = e.second;
	//		for (int alpha = 1; alpha <= delta; alpha++) {
	//			if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
	//			if (g.left_index[u].size() > alpha) {
	//				pi_A[alpha] = pi_A[alpha] >= g.left_index[u][alpha] ? pi_A[alpha] : g.left_index[u][alpha];
	//			}
	//			int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
	//			pi_A[alpha] = pi_A[alpha] >= vbeta ? pi_A[alpha] : vbeta;
	//		}
	//	}
	//	inv(g);
	//	pi_B.resize(delta + 1); fill_n(pi_B.begin(), pi_B.size(), 0);
	//	for (auto e : R) {
	//		vid_t v = e.first; vid_t u = e.second;
	//		for (int alpha = 1; alpha <= delta; alpha++) {
	//			if (g.degree_v1[u] < alpha) break; // has no effect on abcore whose a > alpha
	//			if (g.left_index[u].size() > alpha) {
	//				pi_B[alpha] = pi_B[alpha] >= g.left_index[u][alpha] ? pi_B[alpha] : g.left_index[u][alpha];
	//			}
	//			int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
	//			pi_B[alpha] = pi_B[alpha] >= vbeta ? pi_B[alpha] : vbeta;
	//		}
	//	}
	//	inv(g);

	//	for (auto e : R) {
	//		vid_t u = e.first; vid_t v = e.second;
	//		if (u < g.num_v1 && v < g.num_v2) {
	//			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
	//			int ss = tmp_neigh_.size();
	//			bool deleted = false;
	//			for (int i = 0; i < ss; i++) {
	//				if (tmp_neigh_[i] == v) {
	//					deleted = true;
	//					g.deleteEdge(u, v);
	//					g.left_index[u].pop_back();
	//					g.right_index[v].pop_back();
	//					break;
	//				}
	//			}
	//			if (!deleted) {
	//				cout << "illegal deletion 1" << endl;
	//				return 0;
	//			}
	//		}
	//		else {
	//			cout << "illegal deletion 2" << endl;
	//			return 0;
	//		}
	//	}
	//	delta = maxkcorenumber_test_optimal(g) + 1;
	//	bigraphsA.clear();
	//	bigraphsB.clear();
	//	threadgroupsA.clear();
	//	threadgroupsB.clear();
	//	gap = ceil(delta * 1.0 / (threads / 2));
	//	for (int i = 0; i < threads / 2; i++) {
	//		bigraphsA.push_back(g);
	//		bigraphsB.push_back(g);
	//	}
	//	start = chrono::system_clock::now();
	//	for (int i = 0; i < threads / 2; i++) {
	//		threadgroupsA.push_back(thread(multithread_removal_batch_A, i * gap + 1, MIN(delta, i * gap + gap), ref(pi_A), ref(bigraphsA[i])));
	//	}
	//	for (int i = 0; i < threads / 2; i++) {
	//		threadgroupsB.push_back(thread(multithread_removal_batch_B, i * gap + 1, MIN(delta, i * gap + gap), ref(pi_B), ref(bigraphsB[i])));
	//	}
	//	for (int i = 0; i < threadgroupsA.size(); i++) {
	//		threadgroupsA[i].join();
	//	}
	//	for (int i = 0; i < threadgroupsB.size(); i++) {
	//		threadgroupsB[i].join();
	//	}
	//	end = chrono::system_clock::now();
	//	elapsed_seconds = end - start;
	//	time += elapsed_seconds.count();
	//	merge_removal_result(g, bigraphsA);
	//	return time;
	//}
}