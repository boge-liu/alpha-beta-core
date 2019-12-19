#pragma once
#ifndef INCREMENTAL_TEST_H
#define INCREMENTAL_TEST_H
#include "dynamic.h"
#include "static.h"
#include <random>

extern void random_incremental_correctness_test(BiGraph& g, int num_insert_edges);
extern bool isIndexEqual(BiGraph& g1, BiGraph& g2);
#endif // !INCREMENTAL_TEST_H
