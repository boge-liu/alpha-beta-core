#include "random.h"
#include <random>
#include <sstream>
#include "gephi.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

void partial_random_graph(int left_size, int right_size, int num_edge, string output_location) {
	random_device rd;
	mt19937 gen(rd());
	BiGraph g;
	g.init(left_size, right_size);
	g.addEdge(3, 0); g.addEdge(3, 1); g.addEdge(3, 2); g.addEdge(3, 3); g.addEdge(3, 4); g.addEdge(3, 5); g.addEdge(3, 6); //g.addEdge(3, 7);
	g.addEdge(0, 4); g.addEdge(1, 4); g.addEdge(2, 4); g.addEdge(4, 4); g.addEdge(5, 4); g.addEdge(6, 4); //g.addEdge(7, 4);
	num_edge -= 13;
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
	g.dir = output_location;
	g.print();
}

/* real edge num is: num_edge + left_size + right_size complexity: O(m^2)
   output_location should end with \\
*/
void total_random_graph(int left_size, int right_size, int num_edge, string output_location) {
	random_device rd;
	mt19937 gen(rd());
	BiGraph g;
	g.init(left_size, right_size);
	uniform_int_distribution<int> dis_right(0, right_size - 1);
	uniform_int_distribution<int> dis_left(0, left_size - 1);
	for (int i = 0; i < num_edge;){
		int u = dis_left(gen);
		int v = dis_right(gen);
		if (!g.isEdge(u, v)) {
			g.addEdge(u, v);
			i++;
		}
	}
	g.dir = output_location;
	g.print();
}

void generate_Erdos_Renyi_random_graph(int left_size, int right_size, int avg_deg) {
	int num_edge = (left_size + right_size) * avg_deg / 2;
	random_device rd;
	mt19937 gen(rd());
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
	stringstream ss;
	ss << "../randomGraph/ER_" << "l_" << left_size << "_r_" << right_size << "_avg_" << avg_deg;
	g.dir = ss.str();
	g.print();
}

int binary_search(vector<int>& arr, int val) {
	int start, middle, end;

	start = 0;
	end = arr.size() - 1;
	middle = (start + end) / 2;
	if (arr[start] > val) {
		return start;
	}
	if (val > arr[end]) {
		return end;
	}
	while (start < end)
	{
		if (arr[middle] < val)
		{
			start = middle + 1;
		}
		else if (val < arr[middle])
		{
			end = middle;
		}
		else
		{
			return middle;
		}

		if (middle == (start + end) / 2)
		{
			break;
		}

		middle = (start + end) / 2;
	}

	return middle;
}

void generate_Barabasi_Albert_random_graph(int left_size, int right_size, double avg_deg) {
	int bonus = 2;
	BiGraph g;
	g.init(left_size, right_size);
	int num_edge = (left_size + right_size) * avg_deg / 2;
	vector<int> left_linking_pro; left_linking_pro.resize(left_size);
	left_linking_pro[0] = g.degree_v1[0] * bonus + 1;
	for (int i = 1; i < left_linking_pro.size(); i++) {
		left_linking_pro[i] = left_linking_pro[i-1] + g.degree_v1[i] * bonus + 1;
	}
	vector<int> right_linking_pro; right_linking_pro.resize(right_size);
	right_linking_pro[0] = g.degree_v2[0] * bonus + 1;
	for (int i = 1; i < right_linking_pro.size(); i++) {
		right_linking_pro[i] = right_linking_pro[i - 1] + g.degree_v2[i] * bonus + 1;
	}
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis_left(1, left_linking_pro[left_linking_pro.size() - 1]);
	uniform_int_distribution<int> dis_right(1, right_linking_pro[right_linking_pro.size() - 1]);
	for (int edges = 0; edges < num_edge;) {		
		int up = dis_left(gen);
		int vp = dis_right(gen);
		vid_t u = binary_search(left_linking_pro, up);
		vid_t v = binary_search(right_linking_pro, vp);
		if (!g.isEdge(u, v)) {
			g.addEdge(u, v);
			edges++;
		}
		// update left_linking_pro and right_linking_pro
		if (edges % 1000 == 0) {
			//cout << edges << endl;
			left_linking_pro[0] = g.degree_v1[0] * bonus + 1;
			for (int i = 1; i < left_linking_pro.size(); i++) {
				left_linking_pro[i] = left_linking_pro[i - 1] + g.degree_v1[i] * bonus + 1;
			}
			right_linking_pro[0] = g.degree_v2[0] * bonus + 1;
			for (int i = 1; i < right_linking_pro.size(); i++) {
				right_linking_pro[i] = right_linking_pro[i - 1] + g.degree_v2[i] * bonus + 1;
			}
			uniform_int_distribution<int> dis_left__(1, left_linking_pro[left_linking_pro.size() - 1]);
			dis_left = dis_left__;
			uniform_int_distribution<int> dis_right__(1, right_linking_pro[right_linking_pro.size() - 1]);
			dis_right = dis_right__;
		}
	}
	//printGraph(g);
	stringstream ss;
	ss << "../randomGraph/BA_" << "l_" << left_size << "_r_" << right_size << "_avg_" << avg_deg;
	g.dir = ss.str();
	g.print();
}

void generate_random_graph_and_stream(string graph_location, string stream_location) {
	int left_node_num = 100; int right_node_num = 100; int init_edge_num = 8000;
	int stream_addition_edge_num = 3000; int stream_deletion_edge_num = 3000;
	random_device rd;
	mt19937 gen(rd());
	BiGraph g;
	g.init(left_node_num, right_node_num);
	uniform_int_distribution<int> dis_right(0, right_node_num - 1);
	uniform_int_distribution<int> dis_left(0, left_node_num - 1);
	for (int i = 0; i < init_edge_num;) {
		int u = dis_left(gen);
		int v = dis_right(gen);
		if (!g.isEdge(u, v)) {
			g.addEdge(u, v);
			i++;
		}
	}
	g.dir = graph_location;
	g.print();
	ofstream edge_stream_file;
	edge_stream_file.open(stream_location);
	uniform_int_distribution<int> add_or_delete(0, 1); // 0 delete, 1 add
	int addition = 0; int deletion = 0;
	while (addition < stream_addition_edge_num || deletion < stream_deletion_edge_num) {
		if (add_or_delete(gen) == 1 && addition < stream_addition_edge_num) {
			int u = dis_left(gen);
			int v = dis_right(gen);
			if (!g.isEdge(u, v)) {
				g.addEdge(u, v);
				edge_stream_file << u << ' ' << v << ' ' << 1 << endl;
				addition++;
			}
		}
		else if (deletion < stream_deletion_edge_num) {
			int u = dis_left(gen);
			int v = dis_right(gen);
			if (g.isEdge(u, v)) {
				g.deleteEdge(u, v);
				edge_stream_file << u << ' ' << v << ' ' << 0 << endl;
				deletion++;
			}
		}
	}
	edge_stream_file.close();
}

void generate_random_addition_stream_for_real_graph(string graph_location, string stream_location, int edge_stream_num) {
	BiGraph g(graph_location);
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis_right(0, g.num_v2 - 1);
	uniform_int_distribution<int> dis_left(0, g.num_v1 - 1);
	ofstream edge_stream_file;
	edge_stream_file.open(stream_location);
	for (int i = 0; i < edge_stream_num; ) {
		vid_t u = dis_left(gen);
		vid_t v = dis_right(gen);
		if (!g.isEdge(u, v)) {
			g.addEdge(u, v);
			edge_stream_file << u << ' ' << v << ' ' << 1 << endl;
			i++;
		}
	}
	edge_stream_file.close();
}

void generate_random_deletion_stream_for_real_graph(string graph_location, string stream_location, int edge_stream_num) {
	BiGraph g(graph_location);
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dis_right(0, g.num_v2 - 1);
	uniform_int_distribution<int> dis_left(0, g.num_v1 - 1);
	ofstream edge_stream_file;
	edge_stream_file.open(stream_location);
	for (int i = 0; i < edge_stream_num; ) {
		vid_t u = dis_left(gen);
		if (g.neighbor_v1[u].size() == 0) {
			continue;
		}
		uniform_int_distribution<int> dis_neigh(0, g.neighbor_v1[u].size() - 1);
		vid_t tt = dis_neigh(gen);
		vid_t v = g.neighbor_v1[u][tt];
		g.deleteEdge(u, v);
		edge_stream_file << u << ' ' << v << ' ' << 0 << endl;
		i++;
	}
	edge_stream_file.close();
}