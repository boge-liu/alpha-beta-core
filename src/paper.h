#pragma once
#ifndef PAPER_H
#define PAPER_H
#include "bigraph.h"
#include "incremental_test.h"
#include "dynamic.h"
#include "Multithread.h"
#include "dyn_rebuild.h"
#include <fstream>

extern bool DELETEC;
extern bool APC;
extern long long DELETECOUNT;
extern long long APCOUNT;

extern void alphaCopyPeel_for_kcore(int left_k, BiGraph& g);

extern void __coreIndexBasic(BiGraph& g);

extern void coreIndexBasic(BiGraph& g);

extern void coreIndexBasic_del_count(BiGraph& g);

extern void coreIndexBasicWithCopyPeel_del_count(BiGraph& g);

extern void __coreIndexSwap(BiGraph& g);

extern void coreIndexSwap(BiGraph& g);

extern void coreIndexSwap_del_count(BiGraph& g);

extern int coreIndexKCore_REAL(BiGraph& g);

extern int coreIndexKCore(BiGraph& g);

extern int coreIndexKCore_with_pruning_rules(BiGraph& g);

extern int coreIndexKCore_del_count(BiGraph& g);

extern int coreIndexKCore_with_alpha_peel_and_cross_updating_del_count(BiGraph& g);

extern void alpha_copy_peel_vs_alpha_peel();

extern void findabcore_core_index(int alpha, int beta, BiGraph& g, std::vector<bool>& left_sub, std::vector<bool>& right_sub);

extern void findabcore_peeling(int alpha, int beta, BiGraph& g, std::vector<bool>& left_sub, std::vector<bool>& right_sub);

extern bool test_abcore_equal(std::vector<bool>& left1, std::vector<bool>& right1, std::vector<bool>& left2, std::vector<bool>& right2);

extern void coreIndexSwap_test_optimal(BiGraph& g);

extern int maxkcorenumber_test_optimal(BiGraph& g);

extern int maxkcorenumber(BiGraph& g);

extern void index_construction_test(std::string graph_location);

extern void test_correctness_of_bicore_index_space_saver(std::string graph_location);

// naive increamental without +1

extern void dynamic_test(std::string graph_location, std::string stream_location, Index_update iu);

#endif // !GEPHI_H