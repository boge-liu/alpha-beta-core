#pragma once
#ifndef MISCHIEVOUS_H
#define MISCHIEVOUS_H
#include "bigraph.h"
extern std::vector<std::string> DATASET;

extern std::vector<std::string> ALIAS;

extern void mischievous();

extern void mischievousReal();

extern void test_bicore_index_correctness(std::string graph_location, int alpha, int beta);

#endif // !MISCHIEVOUS_H