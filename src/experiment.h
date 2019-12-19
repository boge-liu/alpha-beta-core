#pragma once
#ifndef EXPERIMENT_H
#define EXPERIMENT_H
#include "bigraph.h"
#include "Multithread.h"

extern std::vector<std::string> DATASET;

extern std::vector<std::string> ALIAS;

extern int DYNAMIC_TEST_NUM;

extern bool DELETEC;
extern bool APC;
extern long long DELETECOUNT;
extern long long APCOUNT;

extern void core2size3d(int i);

extern void timetest_basic();

extern void timetest_basic_single(std::string graph_name, std::string graph_alias);

extern void timetest_swap();

extern void timetest_swap_single(std::string graph_name, std::string graph_alias);

extern void timetest_kcore();

extern void timetest_kcore_single(std::string graph_name, std::string graph_alias);

extern void scalability_sample();

extern void scalabilityTest_basic(int ds);

extern void scalabilityTest_swap(int ds);

extern void scalabilityTest_kcore(int ds);

extern void coresize();

extern void coresize_cutoff(int ds, int cut_begin, int cut_end);

extern void peeltimes();

extern void bicore_index_vs_peeling();

extern void core_index_vs_peeling();

extern void generate_ER();

extern void generate_BA();

extern void vary_average_deg_basic();

extern void vary_average_deg_swap();

extern void vary_average_deg_kcore();

extern void graph_stas();

extern void memory_size_and_retrieve_test();

extern void memory_size_and_all_dataset();

extern void index_detail_print(std::string GRAPH_NAME);

extern void generate_random_graph_and_stream(std::string graph_location, std::string stream_location);

extern void generate_edge_streams();

extern void dynamic_test_all_addition();

extern void dynamic_test_all_deletion();

extern void dynamic_test_scala_vary_node_addition();

extern void dynamic_test_scala_vary_node_deletion();

extern void dynamic_test_scala_vary_edge_addition();

extern void dynamic_test_scala_vary_edge_deletion();

extern void dynamic_degree_out_put(std::string graph_location, std::string edge_stream_location);

extern void dynamic_test_scala_vary_node_deletion_single_test();

extern double dynamic_test_removal_single_test(Datasets ds, Index_update iu, int threads = 4);

extern void dynamic_test_removal_single_test_for_optimize_purpose(Datasets ds, Index_update iu, int threads = 4);

extern double dynamic_test_insertion_single_test(Datasets ds, Index_update iu, int threads = 4);

extern void dynamic_test_insertion_single_test_for_optimize_purpose(Datasets ds, Index_update iu, int threads = 4);

extern void dynamic_test_all_insertion();

extern void dynamic_test_all_removal();

extern void dynamic_test_all_batch();

extern void dynamic_test_scala_vary_node_insertion_single_test(Datasets ds, Index_update iu, int threads = 4);

extern void dynamic_test_scala_vary_node_removal_single_test(Datasets ds, Index_update iu, int threads = 4);

extern void dynamic_test_scala_vary_edge_insertion_single_test(Datasets ds, Index_update iu, int threads = 4);

extern void dynamic_test_scala_vary_edge_removal_single_test(Datasets ds, Index_update iu, int threads = 4);

extern void dynamic_test_scala_vary_node_insertion();

extern void dynamic_test_scala_vary_node_removal();

extern void dynamic_test_scala_vary_edge_insertion();

extern void dynamic_test_scala_vary_edge_removal();

extern void index_construction_scalability();

extern void parallel_vary_core_insertion();

extern void parallel_vary_core_removal();

// ext use

extern void ext_exp(std::string graph_location, int threads, std::string output_location);

#endif // !GEPHI_H