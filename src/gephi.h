#pragma once
#ifndef GEPHI_H
#define GEPHI_H
#include "bigraph.h"

extern void outputKKCore(std::string output_location, BiGraph& g, int left_k, int right_k);

extern void printGraph(BiGraph& g);
#endif // !GEPHI_H
