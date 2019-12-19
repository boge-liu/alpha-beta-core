#pragma once
#ifndef STATIC_H
#define STATIC_H
#include "bigraph.h"

extern void static_construct_index_with_peeling_algorithm_basic(BiGraph& g);

extern void static_construct_index_with_convergence_algorithm_basic(BiGraph& g);

extern void static_construct_index_with_peeling_algorithm_acc(BiGraph& g);

extern void static_construct_index_with_convergence_algorithm_acc(BiGraph& g);

extern void static_construct_index_with_peeling_algorithm_long_yuan(BiGraph& g);

extern void static_construct_index_with_peeling_algorithm_acc_no_edge_deletion(BiGraph& g);

extern void static_construct_index_with_peeling_algorithm_long_yuan_no_edge_deletion(BiGraph& g);

extern void static_construct_index_with_convergence_algorithm_long_yuan(BiGraph& g);
#endif