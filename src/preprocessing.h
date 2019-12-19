#pragma once
#ifndef PREPROCESSING_H
#define PREPROCESSING_H
#include "bigraph.h"
#include <iostream>
#include <fstream>
#include <map>
extern void preprocessing_raw_graph(std::string input_file_location, std::string output_file_location);

extern void preprocessing_raw_graph_amazon_ratings(std::string input_file_location, std::string output_file_location);

extern void preprocessing_raw_graph_live_journal(std::string input_file_location, std::string output_file_location);

extern void preprocessing_raw_graph_with_duplicate_edges(std::string input_file_location, std::string output_file_location);

extern void preprocessing_raw_graph_with_duplicate_edges_real(std::string input_file_location, std::string output_file_location);
#endif // !PREPROCESSING_H
