#include "dynamic.h"
#include <cmath>
#include <algorithm>
using namespace std;

Subgraph::Subgraph() {
	left_side_subgraph_neighbors.clear();
	right_side_subgraph_neighbors.clear();
	left2score.clear();
	right2score.clear();
	left_remain.clear();
	right_remain.clear();
}

bool Subgraph::selfcheck() {
	for (auto it = left_side_subgraph_neighbors.begin(); it != left_side_subgraph_neighbors.end(); it++) {
		if (left2score.find(it->first) == left2score.end()) return false;
		if (left_remain.find(it->first) == left_remain.end()) return false;
	}
	for (auto it = right_side_subgraph_neighbors.begin(); it != right_side_subgraph_neighbors.end(); it++) {
		if (right2score.find(it->first) == right2score.end()) return false;
		if (right_remain.find(it->first) == right_remain.end()) return false;
	}
	return true;
}

void _update_x_k_with_k_x(BiGraph& g, int left_k, int k_x, vid_t v) {
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

// consider if v is newly added!!!!!!!
int calculate_right_side_k_x_with_fixed_left_k(BiGraph& g, int left_k, vid_t v) {
	int right_k = 1;
	for (right_k = 1; right_k < g.right_index[v].size(); right_k++) {
		int l_k = g.right_index[v][right_k];
		if (l_k < left_k) break;
	}
	return right_k - 1;
}

int calculate_u_lower_bound(BiGraph& g, int left_k, vid_t u) {
	vector<int> bins;
	bins.resize(g.v2_max_degree + 1);
	fill_n(bins.begin(), bins.size(), 0);
	vector<vid_t> neigh = g.getV1Neighbors(u);
	int degree = g.getV1Degree(u);
	for (int i = 0; i < degree; i++) {
		// right_k is zero if v(*,1) < left_k
		int right_k = calculate_right_side_k_x_with_fixed_left_k(g, left_k, neigh[i]);
		bins[right_k]++;
	}
	int right_k = 0;
	for (right_k = bins.size() - 1; right_k > 0; right_k--) {
		if (bins[right_k] >= left_k) {
			break;
		}
		bins[right_k - 1] += bins[right_k];
	}
	return right_k;
}

int compute_score_for_left_side_with_fixed_left_k(BiGraph& g, unordered_map<vid_t, int>& right_side2k_x, int left_k, vid_t u) {
	//int right_k = g.left_index[u][left_k];
	int right_k = g.get_left_index_with_fixed_left_k(u, left_k);
	int score = 0;
	vector<vid_t> neigh = g.getV1Neighbors(u);
	int degree = g.getV1Degree(u);
	for (int i = 0; i < degree; i++) {
		vid_t v = neigh[i];
		if (right_side2k_x.find(v) == right_side2k_x.end()) {
			right_side2k_x[v] = calculate_right_side_k_x_with_fixed_left_k(g, left_k, v);
		}
		if (right_k <= right_side2k_x[v]) {
			score++;
		}
	}
	return score;
}

int compute_score_for_right_side_with_fixed_left_k(BiGraph& g, unordered_map<vid_t, int>& right_side2k_x, int left_k, vid_t v) {
	if (right_side2k_x.find(v) == right_side2k_x.end()) {
		right_side2k_x[v] = calculate_right_side_k_x_with_fixed_left_k(g, left_k, v);
	}
	int right_k = right_side2k_x[v];
	int score = 0;
	vector<vid_t> neigh = g.getV2Neighbors(v);
	int degree = g.getV2Degree(v);
	for (int i = 0; i < degree; i++) {
		vid_t u = neigh[i];
		if (right_k <= g.get_left_index_with_fixed_left_k(u, left_k)) {
			score++;
		}
	}
	return score;
}

/*void add_left_side_neighbors_to_subgraph(BiGraph& g, vector<int>& v_tmp_k_star, int left_k, vid_t u, queue<pair<vid_t, bool>>& q, ) {
	left_visited[u] = true;
	vector<vid_t> subneigh;
	vector<vid_t> neigh = g.getV1Neighbors(u);
	int degree = g.getV1Degree(u);
	for (int i = 0; i < degree; i++) {
		if (v_tmp_k_star)
	}
}*/

void find_vertices_with_right_k_increased(Subgraph& subgraph, int left_k, int k_x_selected, vid_t start, bool start_from_left_side) {
	bool unchanged = false;
	while (!unchanged) {
		unchanged = true;
		for (auto it = subgraph.left2score.begin(); it != subgraph.left2score.end(); it++) {
			vid_t u = it->first;
			if (!subgraph.left_remain[u]) continue;
			if (it->second < left_k) { // new condition!
				unchanged = false;
				subgraph.left_remain[u] = false;
				vector<vid_t> neigh = subgraph.left_side_subgraph_neighbors[u];
				for (int i = 0; i < neigh.size(); i++) {
					vid_t v = neigh[i];
					if (!subgraph.right_remain[v]) continue;
					subgraph.right2score[v]--;
				}
			}
		}
		for (auto it = subgraph.right2score.begin(); it != subgraph.right2score.end(); it++) {
			vid_t v = it->first;
			if (!subgraph.right_remain[v]) continue;
			if (it->second <= k_x_selected) {
				unchanged = false;
				subgraph.right_remain[v] = false;
				vector<vid_t> neigh = subgraph.right_side_subgraph_neighbors[v];
				for (int i = 0; i < neigh.size(); i++) {
					vid_t u = neigh[i];
					if (!subgraph.left_remain[u]) continue;
					subgraph.left2score[u]--;
				}
			}
		}
	}
	// compute connected graph from start
	unordered_set<int> left_increment;
	unordered_set<int> right_increment;
	queue<pair<vid_t, bool>> q;
	if (start_from_left_side) {		
		if (subgraph.left_remain[start]) {
			q.push(make_pair(start, LEFT));
			left_increment.insert(start);
		}
		
	}
	else {
		if (subgraph.right_remain[start]) {
			q.push(make_pair(start, RIGHT));
			right_increment.insert(start);
		}
	}
	while (!q.empty()) {
		pair<vid_t, bool> p = q.front();
		q.pop();
		if (p.second == LEFT) {
			vid_t u = p.first;
			vector<vid_t> neigh = subgraph.left_side_subgraph_neighbors[u];
			for (int i = 0; i < neigh.size(); i++) {
				vid_t v = neigh[i];
				if (subgraph.right_remain[v] && right_increment.find(v) == right_increment.end()) {
					q.push(make_pair(v, RIGHT));
					right_increment.insert(v);
				}
			}
		}
		else {
			vid_t v = p.first;
			vector<vid_t> neigh = subgraph.right_side_subgraph_neighbors[v];
			for (int i = 0; i < neigh.size(); i++) {
				vid_t u = neigh[i];
				if (subgraph.left_remain[u] && left_increment.find(u) == left_increment.end()) {
					q.push(make_pair(u, LEFT));
					left_increment.insert(u);
				}
			}
		}
	}
	for (auto it = subgraph.left_remain.begin(); it != subgraph.left_remain.end(); it++) {
		if (left_increment.find(it->first) != left_increment.end()) {
			subgraph.left_remain[it->first] = true;
		}
		else {
			subgraph.left_remain[it->first] = false;
		}
	}
	for (auto it = subgraph.right_remain.begin(); it != subgraph.right_remain.end(); it++) {
		if (right_increment.find(it->first) != right_increment.end()) {
			subgraph.right_remain[it->first] = true;
		}
		else {
			subgraph.right_remain[it->first] = false;
		}
	}
}

void update_index_with_fixed_left_k(BiGraph& g, vector<int>& v_tmp_k_star, int left_k, vid_t u, vid_t v) {
	unordered_map<vid_t, int> right_side2k_x;
	Subgraph subgraph;
	vector<bool> left_visited;
	vector<bool> right_visited;
	left_visited.resize(g.getV1Num());
	fill_n(left_visited.begin(), left_visited.size(), false);
	right_visited.resize(g.getV2Num());
	fill_n(right_visited.begin(), right_visited.size(), false);
	int k_x_selected = -1;
	// subgraph starts from u
	if (g.left_index[u][left_k] <= v_tmp_k_star[left_k]) {
		k_x_selected = g.left_index[u][left_k];		
		/*if (k_x_selected == 0) { // left is 0 and right maybe 0
			subgraph.left_remain[u] = true;
			return;
		}*/
		queue<pair<vid_t, bool>> q; // bool = false is left side otherwise right side
		q.push(make_pair(u, LEFT));
		left_visited[u] = true;
		while (!q.empty()) {
			pair<vid_t, bool> p = q.front();
			q.pop();
			// left side node
			if (p.second == LEFT) {
				vid_t u = p.first;
				vector<vid_t> subneigh;
				vector<vid_t> neigh = g.getV1Neighbors(u);
				int degree = g.getV1Degree(u);
				for (int i = 0; i < degree; i++) {
					vid_t v = neigh[i];
					if (right_side2k_x.find(v) == right_side2k_x.end()) {
						int right_side_k_x = calculate_right_side_k_x_with_fixed_left_k(g, left_k, v);
						right_side2k_x[v] = right_side_k_x;
					}
					int k_x = right_side2k_x[v];
					if (k_x == k_x_selected) {
						subneigh.push_back(v);
						if (!right_visited[v]) {
							q.push(make_pair(v, RIGHT));
							right_visited[v] = true;
						}
					}					
				}
				int score = compute_score_for_left_side_with_fixed_left_k(g, right_side2k_x, left_k, u);
				subgraph.left2score[u] = score;
				subgraph.left_side_subgraph_neighbors[u] = subneigh;
				subgraph.left_remain[u] = true;
			}
			// right side node
			else {
				vid_t v = p.first;
				vector<vid_t> subneigh;
				vector<vid_t> neigh = g.getV2Neighbors(v);
				int degree = g.getV2Degree(v);
				for (int i = 0; i < degree; i++) {
					vid_t u = neigh[i];
					int k_x = g.get_left_index_with_fixed_left_k(u, left_k);
					if (k_x == k_x_selected) {
						subneigh.push_back(u);
						if (!left_visited[u]) {
							q.push(make_pair(u, LEFT));
							left_visited[u] = true;
						}
					}
				}
				int score = compute_score_for_right_side_with_fixed_left_k(g, right_side2k_x, left_k, v);
				subgraph.right2score[v] = score;
				subgraph.right_side_subgraph_neighbors[v] = subneigh;
				subgraph.right_remain[v] = true;
			}
		}
		// sub graph created!
		find_vertices_with_right_k_increased(subgraph, left_k, k_x_selected, u, true);
	}
	// subgraph starts from v
	else {
		k_x_selected = v_tmp_k_star[left_k];
		if (k_x_selected == 0) { // left is not zero right is zero
			subgraph.right_remain[v] = true;
			for (auto it = subgraph.left_remain.begin(); it != subgraph.left_remain.end(); it++) {
				if (it->second) {
					g.left_index[it->first][left_k] = k_x_selected + 1;
				}
			}
			for (auto it = subgraph.right_remain.begin(); it != subgraph.right_remain.end(); it++) {
				if (it->second) {
					_update_x_k_with_k_x(g, left_k, k_x_selected + 1, it->first);
				}
			}
			return;
		}
		queue<pair<vid_t, bool>> q; // bool = false is left side otherwise right side
		q.push(make_pair(v, RIGHT));
		right_visited[v] = true;
		while (!q.empty()) {
			pair<vid_t, bool> p = q.front();
			q.pop();
			// left side node
			if (p.second == LEFT) {
				vid_t u = p.first;
				vector<vid_t> subneigh;
				vector<vid_t> neigh = g.getV1Neighbors(u);
				int degree = g.getV1Degree(u);
				for (int i = 0; i < degree; i++) {
					vid_t v = neigh[i];
					if (right_side2k_x.find(v) == right_side2k_x.end()) {
						int right_side_k_x = calculate_right_side_k_x_with_fixed_left_k(g, left_k, v);
						right_side2k_x[v] = right_side_k_x;
					}
					int k_x = right_side2k_x[v];
					if (k_x == k_x_selected) {
						subneigh.push_back(v);
						if (!right_visited[v]) {
							q.push(make_pair(v, RIGHT));
							right_visited[v] = true;
						}
					}
				}
				int score = compute_score_for_left_side_with_fixed_left_k(g, right_side2k_x, left_k, u);
				subgraph.left2score[u] = score;
				subgraph.left_side_subgraph_neighbors[u] = subneigh;
				subgraph.left_remain[u] = true;
			}
			// right side node
			else {
				vid_t v = p.first;
				vector<vid_t> subneigh;
				vector<vid_t> neigh = g.getV2Neighbors(v);
				int degree = g.getV2Degree(v);
				for (int i = 0; i < degree; i++) {
					vid_t u = neigh[i];
					int k_x = g.get_left_index_with_fixed_left_k(u, left_k);
					if (k_x == k_x_selected) {
						subneigh.push_back(u);
						if (!left_visited[u]) {
							q.push(make_pair(u, LEFT));
							left_visited[u] = true;
						}
					}
				}
				int score = compute_score_for_right_side_with_fixed_left_k(g, right_side2k_x, left_k, v);
				subgraph.right2score[v] = score;
				subgraph.right_side_subgraph_neighbors[v] = subneigh;
				subgraph.right_remain[v] = true;
			}
		}
		// sub graph created!
		find_vertices_with_right_k_increased(subgraph, left_k, k_x_selected, v, false);
	}
	for (auto it = subgraph.left_remain.begin(); it != subgraph.left_remain.end(); it++) {
		if (it->second) {
			g.left_index[it->first][left_k] = k_x_selected + 1;
		}
	}
	for (auto it = subgraph.right_remain.begin(); it != subgraph.right_remain.end(); it++) {
		if (it->second) {
			_update_x_k_with_k_x(g, left_k, k_x_selected + 1, it->first);
		}
	}
}

// find the entire subgraph, adjust score using converagence algorithm, increase core value finally
void dynamic_KKCore_incremental_basic(BiGraph& g, vid_t u, vid_t v) {
	if (g.isEdge(u, v)) {
		return;
	}
	g.addEdge(u, v);
	g.left_index[u].push_back(0);
	for (int left_k = 1; left_k <= g.getV1Degree(u); left_k++) {
		int lb = calculate_u_lower_bound(g, left_k, u);
		g.left_index[u][left_k] = lb;
	}
	

	g.right_index[v].push_back(0);

	vector<int> v_tmp_k_star;
	v_tmp_k_star.resize(g.getV1Degree(u) + 1);
	fill_n(v_tmp_k_star.begin(), v_tmp_k_star.size(), 0);
	for (int left_k = 1; left_k < v_tmp_k_star.size(); left_k++) {
		v_tmp_k_star[left_k] = calculate_right_side_k_x_with_fixed_left_k(g, left_k, v);
	}
	
	for (int left_k = 1; left_k <= g.getV1Degree(u); left_k++) {
		update_index_with_fixed_left_k(g, v_tmp_k_star, left_k, u, v);
	}
	// compute right side
}

void dynamic_KKCore_incremental_eviction(BiGraph& g, vid_t u, vid_t v) {
	
}

//////////////////////////////////////////////////////////////////////
/////////////////////// paper dynamic part ///////////////////////////
//////////////////////////////////////////////////////////////////////

// error if vbeta < 0
int compute_right_beta_with_fixed_alpha(BiGraph& g, int alpha, vid_t v) {
	int vbeta = -1;
	for (int b = 1; b < g.degree_v2[v]; b++) {
		if (g.right_index[v][b] < alpha) {
			vbeta = b - 1;
		}
	}
	return vbeta;
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

int calculate_bound_for_left_node(BiGraph& g, vid_t u, int alpha) {
	// deletion operation when removing level
	if (alpha > g.degree_v1[u]) {
		return g.left_index[u][alpha];
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

// correct based on the fact that sat increases at most 1 even if we set COREdegree,beta to 1
void rearrange_bicore_index_addition(vector<vector<bicore_index_block*>>& bicore_index, int fat, int oldsat, int newsat, vid_t target_node) {
	// extend fat
	if (bicore_index.size() == 0) {
		bicore_index.resize(1);
	}
	if (bicore_index.size() == fat) {
		bicore_index_block* newblock_ = new bicore_index_block;
		newblock_->nodeset.push_back(target_node);
		newblock_->next = NULL;
		vector<bicore_index_block*> newsat_arr; newsat_arr.resize(1); newsat_arr.push_back(newblock_);
		bicore_index.push_back(newsat_arr);
	}
	// the first time to this level (blind addition no need to consider over limit)
	else if (oldsat < 1) {
		if (bicore_index[fat].size() == newsat) {
			bicore_index_block* newblock_ = new bicore_index_block;
			newblock_->next = NULL;
			bicore_index[fat].push_back(newblock_);

		}
		bicore_index_block* newblock = bicore_index[fat][newsat];
		newblock->nodeset.push_back(target_node);
		if (newblock->nodeset.size() == 1) {
			for (int oo = newsat - 1; oo > 0; oo--) {
				if (bicore_index[fat][oo]->nodeset.size() == 0) {
					bicore_index[fat][oo]->next = newblock;
				}
				else {
					bicore_index[fat][oo]->next = newblock;
					break;
				}
			}
		}
	}
	else {
		if (bicore_index[fat].size() == newsat) {
			bicore_index_block* newblock_ = new bicore_index_block;
			newblock_->next = NULL;
			bicore_index[fat].push_back(newblock_);

		}
		bicore_index_block* newblock = bicore_index[fat][newsat];
		bicore_index_block* oldblock = bicore_index[fat][oldsat];
		newblock->nodeset.push_back(target_node);
		auto it = find(oldblock->nodeset.begin(), oldblock->nodeset.end(), target_node);
		oldblock->nodeset.erase(it);
		if (newblock->nodeset.size() == 1) {
			for (int oo = newsat - 1; oo > 0; oo--) {
				if (bicore_index[fat][oo]->nodeset.size() == 0) {
					bicore_index[fat][oo]->next = newblock;
				}
				else {
					bicore_index[fat][oo]->next = newblock;
					break;
				}
			}
		}
		if (oldblock->nodeset.size() == 0) {
			for (int oo = oldsat - 1; oo > 0; oo--) {
				if (bicore_index[fat][oo]->nodeset.size() == 0) {
					bicore_index[fat][oo]->next = oldblock->next;
				}
				else {
					bicore_index[fat][oo]->next = oldblock->next;
					break;
				}
			}
		}
	}
}

void rearrange_bicore_index_deletion(vector<vector<bicore_index_block*>>& bicore_index, int fat, int oldsat, int newsat, vid_t target_node) {
	bicore_index_block* oldblock = bicore_index[fat][oldsat];
	auto it = find(oldblock->nodeset.begin(), oldblock->nodeset.end(), target_node);
	oldblock->nodeset.erase(it);
	// update newsat
	if (newsat >= 1) {
		bicore_index_block* newblock = bicore_index[fat][newsat];
		newblock->nodeset.push_back(target_node);
		if (newblock->nodeset.size() == 1) {
			for (int oo = newsat - 1; oo > 0; oo--) {
				if (bicore_index[fat][oo]->nodeset.size() == 0) {
					bicore_index[fat][oo]->next = newblock;
				}
				else {
					bicore_index[fat][oo]->next = newblock;
					break;
				}
			}
		}
	}
	// delete bicore index
	if (oldblock->nodeset.size() == 0) {
		// delete the last empty block
		if (oldsat == bicore_index[fat].size() - 1) {
			int j = -1;
			for (j = oldsat; j > 0; j--) {
				if (bicore_index[fat][j]->nodeset.empty()) {
					delete bicore_index[fat][j];
					bicore_index[fat].pop_back();
				}
				else {
					break;
				}
			}
			// delete empty level
			if (bicore_index[fat].size() == 1) {
				//delete bicore_index[fat][0];
				bicore_index.pop_back();
			}
			else {
				bicore_index[fat][j]->next = NULL;
			}
		}
		else {
			for (int oo = oldsat - 1; oo > 0; oo--) {
				if (bicore_index[fat][oo]->nodeset.size() == 0) {
					bicore_index[fat][oo]->next = oldblock->next;
				}
				else {
					bicore_index[fat][oo]->next = oldblock->next;
					break;
				}
			}
		}
	}
}

void dyn_crossUpdate_addition(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, int k_x, vid_t v) {
	for (int beta = k_x; beta > 0; beta--) {
		int oldalpha = g.right_index[v][beta];
		if (oldalpha < alpha) {
			g.right_index[v][beta] = alpha;
			rearrange_bicore_index_addition(bicore_index_v, beta, oldalpha, alpha, v);
		}
		else {
			break;
		}
	}
}

void dyn_crossUpdate_deletion(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, int k_x, vid_t v) {
	int newalpha = alpha - 1;
	int truedegree = g.neighbor_v2[v].size();
	for (int i = k_x + 1; i <= truedegree; i++) {
		int oldalpha = g.right_index[v][i];
		if (oldalpha > newalpha) {
			g.right_index[v][i] = newalpha;
			rearrange_bicore_index_deletion(bicore_index_v, i, oldalpha, newalpha, v);
		}
		else {
			break;
		}
	}
}


// without +1 limit
void alphaPeel_for_addition_without_limit(vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int left_k, int beta, BiGraph& g) {
	int dd_;
	int pre_left_k_ = left_k - 1;
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;
	int right_remain_nodes_num = g.num_v2;
	vector<vid_t> right_remain_nodes; right_remain_nodes.resize(g.num_v2);
	for (int i = 0; i < right_remain_nodes.size(); i++) {
		right_remain_nodes[i] = i;
	}
	int right_remain_nodes_tmp_num = 0;
	vector<vid_t> right_remain_nodes_tmp; right_remain_nodes_tmp.resize(g.num_v2);
	for (int right_k = beta + 1; right_k <= g.v2_max_degree + 1; right_k++) {
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
					if (pre_ > beta && dd_ == 0) {
						// core part
						dyn_crossUpdate_addition(g, bicore_index_v, left_k, pre_, v);
						g.right_delete[v] = true;
					}
					if (dd_ == pre_) {
						right_vertices_to_be_peeled.push_back(v);
					}
				}
				g.degree_v1[u] = 0;
				g.left_delete[u] = true;
				if (pre_ > beta) {
					// core part
					int oldbeta = g.left_index[u][left_k];
					if (pre_ > oldbeta) {
						g.left_index[u][left_k] = pre_;
						rearrange_bicore_index_addition(bicore_index_u, left_k, oldbeta, pre_, u);
					}
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
					if (pre_ > beta && dd_ == 0) {
						g.left_delete[u] = true;
						int oldbeta = g.left_index[u][left_k];
						if (pre_ > oldbeta) {
							// core part
							g.left_index[u][left_k] = pre_;
							rearrange_bicore_index_addition(bicore_index_u, left_k, oldbeta, pre_, u);
						}
					}
					if (dd_ == pre_left_k_) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.degree_v2[v] = 0;
				g.right_delete[v] = true;
				if (pre_ > beta) {
					// core part
					dyn_crossUpdate_addition(g, bicore_index_v, left_k, pre_, v);
				}
			}
			right_vertices_to_be_peeled.clear();
		}
	}
}

void update_index_with_fixed_left_k_addition_without_limit(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, vid_t u, vid_t v) {
	int ubeta = g.left_index[u][alpha];
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
	// compute subgraph
	// intialize g.left_delete and g.right_deletion
	fill_n(g.left_delete.begin(), g.left_delete.size(), true);
	fill_n(g.right_delete.begin(), g.right_delete.size(), true);
	int beta = ubeta > vbeta ? vbeta : ubeta;
	if (beta > 0) {
		retrieve_via_bicore_index_inverse(g, bicore_index_u, bicore_index_v, g.left_delete, g.right_delete, alpha, beta);
	}
	else {
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (g.degree_v1[u] >= alpha) {
				g.left_delete[u] = false;
				vector<vid_t>& neigh = g.neighbor_v1[u];
				int ww = neigh.size();
				for (int i = 0; i < ww; i++) {
					g.right_delete[neigh[i]] = false;
				}
			}
		}
	}
	for (vid_t u_ = 0; u_ < g.num_v1; u_++) {
		if (!g.left_delete[u_]) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u_];
			int ss = tmp_neigh_.size();
			for (int k = 0; k < ss; k++) {
				if (g.right_delete[tmp_neigh_[k]]) {
					g.degree_v1[u_]--;
				}
			}
		}
	}
	for (vid_t v_ = 0; v_ < g.num_v2; v_++) {
		if (!g.right_delete[v_]) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v_];
			int ss = tmp_neigh_.size();
			for (int k = 0; k < ss; k++) {
				if (g.left_delete[tmp_neigh_[k]]) {
					g.degree_v2[v_]--;
				}
			}
		}
	}
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
	// peel subgraph
	alphaPeel_for_addition_without_limit(bicore_index_u, bicore_index_v, alpha, beta, g);
	// restore g
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

bool deletion_alphaPeel_with_update_index(vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int left_k, int beta, BiGraph& g, vid_t du, vid_t dv) {
	bool legal = true;
	int dd_;
	int pre_left_k_ = left_k - 1;
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;
	int right_k = beta;
	int pre_ = right_k - 1;
	bool stop = true;
	// incorrect in swap
	/*if (!g.left_delete[du]) {
		if (g.degree_v1[du] < left_k) {
			left_vertices_to_be_peeled.push_back(du);
		}
	}
	if (!g.right_delete[dv]) {
		if (g.degree_v2[dv] < right_k) {
			right_vertices_to_be_peeled.push_back(dv);
		}
	}*/
	for (vid_t u = 0; u < g.num_v1; u++) {
		if (!g.left_delete[u]) {
			if (g.degree_v1[u] < left_k) {
				left_vertices_to_be_peeled.push_back(u);
			}
		}
	}
	for (vid_t v = 0; v < g.num_v2; v++) {
		if (!g.right_delete[v]) {
			if (g.degree_v2[v] < right_k) {
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
				dd_ = --g.degree_v2[v];
				if (pre_ < beta && dd_ == 0) {
					// core part
					legal = false;
					dyn_crossUpdate_deletion(g, bicore_index_v, left_k, pre_, v);
					g.right_delete[v] = true;
				}
				if (dd_ == pre_) {
					right_vertices_to_be_peeled.push_back(v);
				}
			}
			g.degree_v1[u] = 0;
			g.left_delete[u] = true;
			// core part
			int oldbeta = g.left_index[u][left_k];
			if (pre_ < oldbeta) {
				if (left_k <= g.neighbor_v1[u].size()) {
					legal = false;
					g.left_index[u][left_k] = pre_;
					rearrange_bicore_index_deletion(bicore_index_u, left_k, oldbeta, pre_, u);
				}
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
				if (pre_ < beta && dd_ == 0) {
					g.left_delete[u] = true;
					int oldbeta = g.left_index[u][left_k];
					if (pre_ < oldbeta) {
						// core part
						if (left_k <= g.neighbor_v1[u].size()) {
							legal = false;
							g.left_index[u][left_k] = pre_;
							rearrange_bicore_index_deletion(bicore_index_u, left_k, oldbeta, pre_, u);
						}
					}
				}
				if (dd_ == pre_left_k_) {
					left_vertices_to_be_peeled.push_back(u);
				}
			}
			g.degree_v2[v] = 0;
			g.right_delete[v] = true;
			legal = false;
			dyn_crossUpdate_deletion(g, bicore_index_v, left_k, pre_, v);
		}
		right_vertices_to_be_peeled.clear();
	}
	return legal;
}

void update_index_with_fixed_left_k_deletion_without_limit(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, vid_t u, vid_t v) {
	int ubeta = g.left_index[u][alpha];
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
	bool flag_ = true;
	// backup g
	vector<int> degree_v1_backup = g.degree_v1;
	vector<int> degree_v2_backup = g.degree_v2;
	int v1_degree_max_backup = g.v1_max_degree;
	int v2_degree_max_backup = g.v2_max_degree;
	while (flag_ && beta != 0) {
		// compute subgraph
		// intialize g.left_delete and g.right_deletion
		fill_n(g.left_delete.begin(), g.left_delete.size(), true);
		fill_n(g.right_delete.begin(), g.right_delete.size(), true);
		retrieve_via_bicore_index_inverse(g, bicore_index_u, bicore_index_v, g.left_delete, g.right_delete, alpha, beta);
		// optimize
		if (g.num_v1 < g.num_v2) {
			fill_n(g.degree_v2.begin(), g.degree_v2.size(), 0);
			for (vid_t u_ = 0; u_ < g.num_v1; u_++) {
				if (!g.left_delete[u_]) {
					vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u_];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						g.degree_v2[tmp_neigh_[k]]++;
						if (g.right_delete[tmp_neigh_[k]]) {
							g.degree_v1[u_]--;
						}
					}
				}
			}
		}
		else {
			fill_n(g.degree_v1.begin(), g.degree_v1.size(), 0);
			for (vid_t v_ = 0; v_ < g.num_v2; v_++) {
				if (!g.right_delete[v_]) {
					vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v_];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						g.degree_v1[tmp_neigh_[k]]++;
						if (g.left_delete[tmp_neigh_[k]]) {
							g.degree_v2[v_]--;
						}
					}
				}
			}
		}
		g.v1_max_degree = 0;
		g.v2_max_degree = 0;
		for (vid_t u = 0; u < g.degree_v1.size(); u++) {
			if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
		}
		for (vid_t v = 0; v < g.degree_v2.size(); v++) {
			if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
		}
		// peel subgraph
		flag_ = !deletion_alphaPeel_with_update_index(bicore_index_u, bicore_index_v, alpha, beta, g, u, v);
		beta--;
		// restore g
		g.v1_max_degree = v1_degree_max_backup; g.v2_max_degree = v2_degree_max_backup;
		g.degree_v1 = degree_v1_backup; g.degree_v2 = degree_v2_backup;
	}
}

void update_bicore_index_without_limit(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, vid_t u, vid_t v, bool addition) {
	// check whether operation is legal
	if (addition) {
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
				return;
			}
		}
		else {
			cout << "illegal deletion 2" << endl;
			return;
		}
	}
	// addition
	if (addition) {
		// extend core index and bicore index
		g.left_index[u].push_back(0);
		g.right_index[v].push_back(0);
		int max_alpha = g.degree_v1[u];
		for (int alpha = 1; alpha <= max_alpha; alpha++) {
			update_index_with_fixed_left_k_addition_without_limit(g, bicore_index_u, bicore_index_v, alpha, u, v);
		}
	}
	// deletion
	else {
		int oldusat = g.left_index[u].back();
		int oldvsat = g.right_index[v].back();
		int max_alpha = g.degree_v1[u];
		for (int alpha = max_alpha + 1; alpha > 0; alpha--) {
			update_index_with_fixed_left_k_deletion_without_limit(g, bicore_index_u, bicore_index_v, alpha, u, v);
		}
		g.left_index[u].pop_back();
		g.right_index[v].pop_back();
		rearrange_bicore_index_deletion(bicore_index_u, g.degree_v1[u] + 1, oldusat, -1, u);
		rearrange_bicore_index_deletion(bicore_index_v, g.degree_v2[v] + 1, oldvsat, -1, v);
	}
}

// with +1 limit
void alphaPeel_for_addition_with_limit(vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int left_k, int beta, BiGraph& g) {
	int dd_;
	int pre_left_k_ = left_k - 1;
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;
	int right_k = beta + 1;
	int pre_ = right_k - 1;
	bool stop = true;
	for (vid_t u = 0; u < g.num_v1; u++) {
		if (!g.left_delete[u]) {
			if (g.degree_v1[u] < left_k) {
				left_vertices_to_be_peeled.push_back(u);
			}
		}
	}
	for (vid_t v = 0; v < g.num_v2; v++) {
		if (!g.right_delete[v]) {
			if (g.degree_v2[v] < right_k) {
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
				dd_ = --g.degree_v2[v];
				if (dd_ == 0) {
					// core part
					g.right_delete[v] = true;
				}
				if (dd_ == pre_) {
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
				dd_ = --g.degree_v1[u];
				if (dd_ == 0) {
					g.left_delete[u] = true;
				}
				if (dd_ == pre_left_k_) {
					left_vertices_to_be_peeled.push_back(u);
				}
			}
			g.degree_v2[v] = 0;
			g.right_delete[v] = true;
		}
		right_vertices_to_be_peeled.clear();
	}
	for (vid_t u = 0; u < g.num_v1; u++) {
		if (!g.left_delete[u]) {
			int oldbeta = g.left_index[u][left_k];
			if (oldbeta < right_k) {
				g.left_index[u][left_k] = right_k;
				rearrange_bicore_index_addition(bicore_index_u, left_k, oldbeta, right_k, u);
			}
		}
	}
	for (vid_t v = 0; v < g.num_v2; v++) {
		if (!g.right_delete[v]) {
			dyn_crossUpdate_addition(g, bicore_index_v, left_k, right_k, v);
		}
	}
}

void update_index_with_fixed_left_k_addition_with_limit(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, vid_t u, vid_t v) {
	int ubeta = g.left_index[u][alpha];
	int bound = calculate_bound_for_left_node(g, u, alpha);
	if (bound > ubeta) {
		g.left_index[u][alpha] = bound;
		rearrange_bicore_index_addition(bicore_index_u, alpha, ubeta, bound, u);
		ubeta = bound;
	}
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
	if (vbeta == -1) {
		cout << "error: vbeta" << endl;
	}
	// compute subgraph
	// intialize g.left_delete and g.right_deletion
	fill_n(g.left_delete.begin(), g.left_delete.size(), true);
	fill_n(g.right_delete.begin(), g.right_delete.size(), true);
	int beta = ubeta > vbeta ? vbeta : ubeta;
	if (beta > 0) {
		retrieve_via_bicore_index_inverse(g, bicore_index_u, bicore_index_v, g.left_delete, g.right_delete, alpha, beta);
	}
	else {
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (g.degree_v1[u] >= alpha) {
				g.left_delete[u] = false;
				vector<vid_t>& neigh = g.neighbor_v1[u];
				int ww = neigh.size();
				for (int i = 0; i < ww; i++) {
					g.right_delete[neigh[i]] = false;
				}
			}
		}
	}
	for (vid_t u_ = 0; u_ < g.num_v1; u_++) {
		if (!g.left_delete[u_]) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u_];
			int ss = tmp_neigh_.size();
			for (int k = 0; k < ss; k++) {
				if (g.right_delete[tmp_neigh_[k]]) {
					g.degree_v1[u_]--;
				}
			}
		}
	}
	for (vid_t v_ = 0; v_ < g.num_v2; v_++) {
		if (!g.right_delete[v_]) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v_];
			int ss = tmp_neigh_.size();
			for (int k = 0; k < ss; k++) {
				if (g.left_delete[tmp_neigh_[k]]) {
					g.degree_v2[v_]--;
				}
			}
		}
	}
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
	// peel subgraph
	alphaPeel_for_addition_with_limit(bicore_index_u, bicore_index_v, alpha, beta, g);
	// restore g
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

void update_index_with_fixed_left_k_deletion_with_limit(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, vid_t u, vid_t v) {
	int ubeta = g.left_index[u][alpha];
	int bound = calculate_bound_for_left_node(g, u, alpha);
	g.left_index[u][alpha] = bound;
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
	bool flag_ = true;
	int oldbeta = beta;
	// backup g
	vector<int> degree_v1_backup = g.degree_v1;
	vector<int> degree_v2_backup = g.degree_v2;
	int v1_degree_max_backup = g.v1_max_degree;
	int v2_degree_max_backup = g.v2_max_degree;
	while (flag_ && beta != 0) {
		// compute subgraph
		// intialize g.left_delete and g.right_deletion
		fill_n(g.left_delete.begin(), g.left_delete.size(), true);
		fill_n(g.right_delete.begin(), g.right_delete.size(), true);
		retrieve_via_bicore_index_inverse(g, bicore_index_u, bicore_index_v, g.left_delete, g.right_delete, alpha, beta);
		// optimize
		if (g.num_v1 < g.num_v2) {
			fill_n(g.degree_v2.begin(), g.degree_v2.size(), 0);
			for (vid_t u_ = 0; u_ < g.num_v1; u_++) {
				if (!g.left_delete[u_]) {
					vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u_];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						g.degree_v2[tmp_neigh_[k]]++;
						if (g.right_delete[tmp_neigh_[k]]) {
							g.degree_v1[u_]--;
						}
					}
				}
			}
		}
		else {
			fill_n(g.degree_v1.begin(), g.degree_v1.size(), 0);
			for (vid_t v_ = 0; v_ < g.num_v2; v_++) {
				if (!g.right_delete[v_]) {
					vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v_];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						g.degree_v1[tmp_neigh_[k]]++;
						if (g.left_delete[tmp_neigh_[k]]) {
							g.degree_v2[v_]--;
						}
					}
				}
			}
		}
		g.v1_max_degree = 0;
		g.v2_max_degree = 0;
		for (vid_t u = 0; u < g.degree_v1.size(); u++) {
			if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
		}
		for (vid_t v = 0; v < g.degree_v2.size(); v++) {
			if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
		}
		// peel subgraph
		flag_ = !deletion_alphaPeel_with_update_index(bicore_index_u, bicore_index_v, alpha, beta, g, u, v);
		beta--;
		// restore g
		g.v1_max_degree = v1_degree_max_backup; g.v2_max_degree = v2_degree_max_backup;
		g.degree_v1 = degree_v1_backup; g.degree_v2 = degree_v2_backup;
		// only one more layer
		if (beta + 1 > oldbeta) {
			break;
		}
	}
	if (bound < ubeta) {
		rearrange_bicore_index_deletion(bicore_index_u, alpha, ubeta, bound, u);
	}
}

void update_bicore_index_with_limit(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, vid_t u, vid_t v, bool addition) {
	// check whether operation is legal
	if (addition) {
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
				return;
			}
		}
		else {
			cout << "illegal deletion 2" << endl;
			return;
		}
	}
	// addition
	if (addition) {
		// extend core index
		g.left_index[u].push_back(0);
		g.right_index[v].push_back(0);
		int max_alpha = g.degree_v1[u];
		for (int alpha = 1; alpha <= max_alpha; alpha++) {
			update_index_with_fixed_left_k_addition_with_limit(g, bicore_index_u, bicore_index_v, alpha, u, v);
		}
	}
	// deletion
	else {
		int oldusat = g.left_index[u].back();
		int oldvsat = g.right_index[v].back();
		int max_alpha = g.degree_v1[u];
		for (int alpha = max_alpha + 1; alpha > 0; alpha--) {
			update_index_with_fixed_left_k_deletion_with_limit(g, bicore_index_u, bicore_index_v, alpha, u, v);
		}
		g.left_index[u].pop_back();
		g.right_index[v].pop_back();
		rearrange_bicore_index_deletion(bicore_index_u, g.degree_v1[u] + 1, oldusat, -1, u);
		rearrange_bicore_index_deletion(bicore_index_v, g.degree_v2[v] + 1, oldvsat, -1, v);
	}
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

void update_index_with_fixed_left_k_addition_without_limit_swap(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, int alpha_start, vid_t u, vid_t v) {
	int beta = alpha_start;
	// compute subgraph
	// intialize g.left_delete and g.right_deletion
	fill_n(g.left_delete.begin(), g.left_delete.size(), true);
	fill_n(g.right_delete.begin(), g.right_delete.size(), true);
	if (beta > 0) {
		retrieve_via_bicore_index_inverse(g, bicore_index_u, bicore_index_v, g.left_delete, g.right_delete, alpha, beta);
	}
	else {
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (g.degree_v1[u] >= alpha) {
				g.left_delete[u] = false;
				vector<vid_t>& neigh = g.neighbor_v1[u];
				int ww = neigh.size();
				for (int i = 0; i < ww; i++) {
					g.right_delete[neigh[i]] = false;
				}
			}
		}
	}
	for (vid_t u_ = 0; u_ < g.num_v1; u_++) {
		if (!g.left_delete[u_]) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u_];
			int ss = tmp_neigh_.size();
			for (int k = 0; k < ss; k++) {
				if (g.right_delete[tmp_neigh_[k]]) {
					g.degree_v1[u_]--;
				}
			}
		}
	}
	for (vid_t v_ = 0; v_ < g.num_v2; v_++) {
		if (!g.right_delete[v_]) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v_];
			int ss = tmp_neigh_.size();
			for (int k = 0; k < ss; k++) {
				if (g.left_delete[tmp_neigh_[k]]) {
					g.degree_v2[v_]--;
				}
			}
		}
	}
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
	// peel subgraph
	alphaPeel_for_addition_without_limit(bicore_index_u, bicore_index_v, alpha, beta, g);
	// restore g
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

void update_index_with_fixed_left_k_deletion_without_limit_swap(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, int alpha_start, vid_t u, vid_t v) {
	int beta = alpha_start;
	bool flag_ = true;
	// backup g
	vector<int> degree_v1_backup = g.degree_v1;
	vector<int> degree_v2_backup = g.degree_v2;
	int v1_degree_max_backup = g.v1_max_degree;
	int v2_degree_max_backup = g.v2_max_degree;
	while (flag_ && beta != 0) {
		// compute subgraph
		// intialize g.left_delete and g.right_deletion
		fill_n(g.left_delete.begin(), g.left_delete.size(), true);
		fill_n(g.right_delete.begin(), g.right_delete.size(), true);
		retrieve_via_bicore_index_inverse(g, bicore_index_u, bicore_index_v, g.left_delete, g.right_delete, alpha, beta);
		// optimize
		if (g.num_v1 < g.num_v2) {
			fill_n(g.degree_v2.begin(), g.degree_v2.size(), 0);
			for (vid_t u_ = 0; u_ < g.num_v1; u_++) {
				if (!g.left_delete[u_]) {
					vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u_];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						g.degree_v2[tmp_neigh_[k]]++;
						if (g.right_delete[tmp_neigh_[k]]) {
							g.degree_v1[u_]--;
						}
					}
				}
			}
		}
		else {
			fill_n(g.degree_v1.begin(), g.degree_v1.size(), 0);
			for (vid_t v_ = 0; v_ < g.num_v2; v_++) {
				if (!g.right_delete[v_]) {
					vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v_];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						g.degree_v1[tmp_neigh_[k]]++;
						if (g.left_delete[tmp_neigh_[k]]) {
							g.degree_v2[v_]--;
						}
					}
				}
			}
		}
		g.v1_max_degree = 0;
		g.v2_max_degree = 0;
		for (vid_t u = 0; u < g.degree_v1.size(); u++) {
			if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
		}
		for (vid_t v = 0; v < g.degree_v2.size(); v++) {
			if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
		}
		// peel subgraph
		flag_ = !deletion_alphaPeel_with_update_index(bicore_index_u, bicore_index_v, alpha, beta, g, u, v);
		beta--;
		// restore g
		g.v1_max_degree = v1_degree_max_backup; g.v2_max_degree = v2_degree_max_backup;
		g.degree_v1 = degree_v1_backup; g.degree_v2 = degree_v2_backup;
	}
}


void update_bicore_index_without_limit_swap(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, vid_t u, vid_t v, bool addition) {
	// calculate swap threshold
	pair<int, int> r = calculate_Delta(bicore_index_u);
	// check whether operation is legal
	if (addition) {
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
				return;
			}
		}
		else {
			cout << "illegal deletion 2" << endl;
			return;
		}
	}
	// addition
	if (addition) {
		// extend core index and bicore index
		g.left_index[u].push_back(0);
		g.right_index[v].push_back(0);
		int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
		int th_alpha, th_beta;
		if (r.first + 1 >= max_alpha) {
			th_alpha = max_alpha; th_beta = 0;
		}
		else if (r.second + 1 >= max_beta) {
			th_alpha = 0; th_beta = max_beta;
		}
		else {
			th_alpha = r.first + 1; th_beta = r.second + 1;
		}
		// compute alpha beta condition in advance
		vector<int> alpha_start; vector<int> beta_start; alpha_start.push_back(-1); beta_start.push_back(-1);
		for (int alpha = 1; alpha <= th_alpha; alpha++) {
			int ubeta = g.left_index[u][alpha];
			int vbeta = -1;
			bool changed = false;
			for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
				if (g.right_index[v][b] < alpha) {
					vbeta = b - 1;
					changed = true;
					break;
				}
			}
			/*if (!changed) {
				vbeta = g.degree_v2[v];
			}
			if (vbeta < 0) {
				cout << "error: vbeta" << endl;
			}*/
			int beta = ubeta > vbeta ? vbeta : ubeta;
			alpha_start.push_back(beta);
		}
		for (int beta = 1; beta <= th_beta; beta++) {
			int valpha = g.right_index[v][beta];
			int ualpha = -1;
			bool changed = false;
			for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
				if (g.left_index[u][a] < beta) {
					ualpha = a - 1;
					changed = true;
					break;
				}
			}
			/*if (!changed) {
				ualpha = g.degree_v1[u];
			}
			if (ualpha < 0) {
				cout << "error: ualpha" << endl;
			}*/
			int alpha = valpha > ualpha ? ualpha : valpha;
			beta_start.push_back(alpha);
		}
		for (int alpha = 1; alpha <= th_alpha; alpha++) {
			update_index_with_fixed_left_k_addition_without_limit(g, bicore_index_u, bicore_index_v, alpha, u, v);
		}
		swap(bicore_index_u, bicore_index_v);
		inv(g);
		for (int alpha = 1; alpha <= th_beta; alpha++) {
			update_index_with_fixed_left_k_addition_without_limit(g, bicore_index_u, bicore_index_v, alpha, v, u);
		}
		swap(bicore_index_u, bicore_index_v);
		inv(g);
	}
	// deletion
	else {
		int oldusat = g.left_index[u].back();
		int oldvsat = g.right_index[v].back();
		int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
		int th_alpha, th_beta;
		if (r.first >= max_alpha) {
			th_alpha = max_alpha; th_beta = 0;
		}
		else if (r.second >= max_beta) {
			th_alpha = 0; th_beta = max_beta;
		}
		else {
			th_alpha = r.first; th_beta = r.second;
		}
		// compute alpha beta condition in advance
		vector<int> alpha_start; vector<int> beta_start; alpha_start.push_back(-1); beta_start.push_back(-1);
		for (int alpha = 1; alpha <= th_alpha; alpha++) {
			int ubeta = g.left_index[u][alpha];
			int vbeta = g.degree_v2[v] + 1;
			bool changed = false;
			// g.degree_v2[v] has already been decreased
			for (int b = 1; b <= g.degree_v2[v] + 1; b++) {
				if (g.right_index[v][b] < alpha) {
					vbeta = b - 1;
					changed = true;
					break;
				}
			}
			/*if (!changed) {
				vbeta = g.degree_v2[v] + 1;
			}
			if (vbeta == -1 || vbeta == 0) {
				cout << "error: vbeta" << endl;
			}*/
			int beta = ubeta > vbeta ? vbeta : ubeta;
			alpha_start.push_back(beta);
		}
		for (int beta = 1; beta <= th_beta; beta++) {
			int valpha = g.right_index[v][beta];
			int ualpha = g.degree_v1[u] + 1;
			bool changed = false;
			for (int a = 1; a <= g.degree_v1[u] + 1; a++) {
				if (g.left_index[u][a] < beta) {
					ualpha = a - 1;
					changed = true;
					break;
				}
			}
			/*if (!changed) {
				ualpha = g.degree_v1[u] + 1;
			}
			if (ualpha == -1 || ualpha == 0) {
				cout << "error: ualpha" << endl;
			}*/
			int alpha = valpha > ualpha ? ualpha : valpha;
			beta_start.push_back(alpha);
		}
		for (int alpha = 1; alpha <= th_alpha; alpha++) {
			update_index_with_fixed_left_k_deletion_without_limit_swap(g, bicore_index_u, bicore_index_v, alpha, alpha_start[alpha], u, v);
		}
		swap(bicore_index_u, bicore_index_v);
		inv(g);
		for (int alpha = 1; alpha <= th_beta; alpha++) {
			update_index_with_fixed_left_k_deletion_without_limit_swap(g, bicore_index_u, bicore_index_v, alpha, beta_start[alpha], v, u);
		}
		swap(bicore_index_u, bicore_index_v);
		inv(g);
		g.left_index[u].pop_back();
		g.right_index[v].pop_back();
		rearrange_bicore_index_deletion(bicore_index_u, g.degree_v1[u] + 1, oldusat, -1, u);
		rearrange_bicore_index_deletion(bicore_index_v, g.degree_v2[v] + 1, oldvsat, -1, v);
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

int calculate_bound_for_right_node(BiGraph& g, vid_t v, int beta) {
	// deletion operation when removing level
	if (beta > g.degree_v2[v]) {
		return g.right_index[v][beta];
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

void update_index_with_fixed_left_k_addition_with_limit_swap(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, int alpha_start, int start_bound, vid_t u, vid_t v) {
	int ubeta = g.left_index[u][alpha];
	int bound = start_bound;
	if (bound > ubeta) {
		g.left_index[u][alpha] = bound;
		rearrange_bicore_index_addition(bicore_index_u, alpha, ubeta, bound, u);
	}
	int beta = alpha_start;
	// compute subgraph
	// intialize g.left_delete and g.right_deletion
	fill_n(g.left_delete.begin(), g.left_delete.size(), true);
	fill_n(g.right_delete.begin(), g.right_delete.size(), true);
	if (beta > 0) {
		retrieve_via_bicore_index_inverse(g, bicore_index_u, bicore_index_v, g.left_delete, g.right_delete, alpha, beta);
	}
	else {
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (g.degree_v1[u] >= alpha) {
				g.left_delete[u] = false;
				vector<vid_t>& neigh = g.neighbor_v1[u];
				int ww = neigh.size();
				for (int i = 0; i < ww; i++) {
					g.right_delete[neigh[i]] = false;
				}
			}
		}
	}
	for (vid_t u_ = 0; u_ < g.num_v1; u_++) {
		if (!g.left_delete[u_]) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u_];
			int ss = tmp_neigh_.size();
			for (int k = 0; k < ss; k++) {
				if (g.right_delete[tmp_neigh_[k]]) {
					g.degree_v1[u_]--;
				}
			}
		}
	}
	for (vid_t v_ = 0; v_ < g.num_v2; v_++) {
		if (!g.right_delete[v_]) {
			vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v_];
			int ss = tmp_neigh_.size();
			for (int k = 0; k < ss; k++) {
				if (g.left_delete[tmp_neigh_[k]]) {
					g.degree_v2[v_]--;
				}
			}
		}
	}
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
	// peel subgraph
	alphaPeel_for_addition_with_limit(bicore_index_u, bicore_index_v, alpha, beta, g);
	// restore g
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

void update_index_with_fixed_left_k_deletion_with_limit_swap(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, int alpha, int alpha_start, int start_bound, vid_t u, vid_t v) {
	int ubeta = g.left_index[u][alpha];
	int beta = alpha_start;
	int bound = start_bound;
	if (g.left_index[u][alpha] > bound) {
		g.left_index[u][alpha] = bound;
	}	
	int vbeta = -1;
	bool changed = false;
	bool flag_ = true;
	int oldbeta = beta;
	// backup g
	vector<int> degree_v1_backup = g.degree_v1;
	vector<int> degree_v2_backup = g.degree_v2;
	int v1_degree_max_backup = g.v1_max_degree;
	int v2_degree_max_backup = g.v2_max_degree;
	while (flag_ && beta != 0) {
		// compute subgraph
		// intialize g.left_delete and g.right_deletion
		fill_n(g.left_delete.begin(), g.left_delete.size(), true);
		fill_n(g.right_delete.begin(), g.right_delete.size(), true);
		//retrieve_via_bicore_index_inverse(g, bicore_index_u, bicore_index_v, g.left_delete, g.right_delete, alpha, beta);
		if (bicore_index_u.size() <= alpha) return;
		if (bicore_index_u[alpha].size() <= beta) return;
		bicore_index_block* block = bicore_index_u[alpha][beta];
		while (block != NULL) {
			for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
				g.left_delete[*i] = false;
			}
			block = block->next;
		}
		block = bicore_index_v[beta][alpha];
		while (block != NULL) {
			for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
				g.right_delete[*i] = false;
			}
			block = block->next;
		}
		// optimize
		if (g.num_v1 < g.num_v2) {
			fill_n(g.degree_v2.begin(), g.degree_v2.size(), 0);
			for (vid_t u_ = 0; u_ < g.num_v1; u_++) {
				if (!g.left_delete[u_]) {
					vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u_];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						g.degree_v2[tmp_neigh_[k]]++;
						if (g.right_delete[tmp_neigh_[k]]) {
							g.degree_v1[u_]--;
						}
					}
				}
			}
		}
		else {
			fill_n(g.degree_v1.begin(), g.degree_v1.size(), 0);
			for (vid_t v_ = 0; v_ < g.num_v2; v_++) {
				if (!g.right_delete[v_]) {
					vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v_];
					int ss = tmp_neigh_.size();
					for (int k = 0; k < ss; k++) {
						g.degree_v1[tmp_neigh_[k]]++;
						if (g.left_delete[tmp_neigh_[k]]) {
							g.degree_v2[v_]--;
						}
					}
				}
			}
		}
		g.v1_max_degree = 0;
		g.v2_max_degree = 0;
		for (vid_t u = 0; u < g.degree_v1.size(); u++) {
			if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
		}
		for (vid_t v = 0; v < g.degree_v2.size(); v++) {
			if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
		}
		// peel subgraph
		flag_ = !deletion_alphaPeel_with_update_index(bicore_index_u, bicore_index_v, alpha, beta, g, u, v);
		beta--;
		// restore g
		g.v1_max_degree = v1_degree_max_backup; g.v2_max_degree = v2_degree_max_backup;
		g.degree_v1 = degree_v1_backup; g.degree_v2 = degree_v2_backup;
		// only one more layer
		if (beta + 1 < oldbeta) {
			break;
		}
	}
	if (bound < ubeta) {
		rearrange_bicore_index_deletion(bicore_index_u, alpha, ubeta, bound, u);
	}
}

void update_bicore_index_with_limit_swap(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, vid_t u, vid_t v, bool addition) {
	// calculate swap threshold
	pair<int, int> r = calculate_Delta(bicore_index_u);
	// check whether operation is legal
	if (addition) {
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
				return;
			}
		}
		else {
			cout << "illegal deletion 2" << endl;
			return;
		}
	}
	// addition
	if (addition) {
		// extend core index and bicore index
		g.left_index[u].push_back(0);
		g.right_index[v].push_back(0);
		int max_alpha = g.degree_v1[u]; int max_beta = g.degree_v2[v];
		int th_alpha, th_beta;
		if (r.first + 1 >= max_alpha) {
			th_alpha = max_alpha; th_beta = 0;
		}
		else if (r.second + 1 >= max_beta) {
			th_alpha = 0; th_beta = max_beta;
		}
		else {
			th_alpha = r.first + 1; th_beta = r.second + 1;
		}
		// compute alpha beta condition in advance
		vector<int> alpha_start; vector<int> beta_start; vector<int> alpha_start_bound; vector<int> beta_start_bound;
		alpha_start.push_back(-1); beta_start.push_back(-1); alpha_start_bound.push_back(-1); beta_start_bound.push_back(-1);
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
				ualpha = g.degree_v1[u];
			}
			if (ualpha < 0) {
				cout << "error: ualpha" << endl;
			}
			int alpha = bound > ualpha ? ualpha : bound;
			beta_start.push_back(alpha);
			beta_start_bound.push_back(bound);
		}
		for (int alpha = 1; alpha <= th_alpha; alpha++) {
			update_index_with_fixed_left_k_addition_with_limit_swap(g, bicore_index_u, bicore_index_v, alpha, alpha_start[alpha], alpha_start_bound[alpha], u, v);
		}
		swap(bicore_index_u, bicore_index_v);
		inv(g);
		for (int alpha = 1; alpha <= th_beta; alpha++) {
			update_index_with_fixed_left_k_addition_with_limit_swap(g, bicore_index_u, bicore_index_v, alpha, beta_start[alpha], beta_start_bound[alpha], v, u);
		}
		swap(bicore_index_u, bicore_index_v);
		inv(g);
	}
	// deletion
	else {
		int oldusat = g.left_index[u].back();
		int oldvsat = g.right_index[v].back();
		int max_alpha = g.left_index[u].size() - 1; int max_beta = g.right_index[v].size() - 1;
		int th_alpha, th_beta;
		if (r.first >= max_alpha) {
			th_alpha = max_alpha; th_beta = 0;
		}
		else if (r.second >= max_beta) {
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
			update_index_with_fixed_left_k_deletion_with_limit_swap(g, bicore_index_u, bicore_index_v, alpha, alpha_start[alpha], alpha_start_bound[alpha], u, v);
		}
		swap(bicore_index_u, bicore_index_v);
		inv(g);
		for (int alpha = 1; alpha <= th_beta; alpha++) {
			update_index_with_fixed_left_k_deletion_with_limit_swap(g, bicore_index_u, bicore_index_v, alpha, beta_start[alpha], beta_start_bound[alpha], v, u);
		}
		swap(bicore_index_u, bicore_index_v);
		inv(g);
		g.left_index[u].pop_back();
		g.right_index[v].pop_back();
		rearrange_bicore_index_deletion(bicore_index_u, g.degree_v1[u] + 1, oldusat, -1, u);
		rearrange_bicore_index_deletion(bicore_index_v, g.degree_v2[v] + 1, oldvsat, -1, v);
	}
}