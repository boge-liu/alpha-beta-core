#pragma once
#ifndef RANDOM_H
#define RANDOM_H
#include "bigraph.h"

extern void partial_random_graph(int left_size, int right_size, int num_edge, std::string output_location);

extern void total_random_graph(int left_size, int right_size, int num_edge, std::string output_location);

extern void generate_Erdos_Renyi_random_graph(int left_size, int right_size, int avg_deg);

extern void generate_Barabasi_Albert_random_graph(int left_size, int right_size, double avg_deg);

extern void generate_random_graph_and_stream(std::string graph_location, std::string stream_location);

extern void generate_random_addition_stream_for_real_graph(std::string graph_location, std::string stream_location, int edge_stream_num);

extern void generate_random_deletion_stream_for_real_graph(std::string graph_location, std::string stream_location, int edge_stream_num);

#endif // !RANDOM_H
