#pragma once
#ifndef MULTITHREAD
#define MULTITHREAD
#include "bigraph.h"
#include "paper.h"

extern void multithread_index_construction(BiGraph& g, int threads);

extern void multithread_test(BiGraph& g, int threads);

extern void pure_multithread_test(BiGraph& g, int threads, std::ofstream& output_file);

#endif // !MULTITHREAD
