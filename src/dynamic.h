#pragma once
#ifndef DYNAMIC_H
#define DYNAMIC_H
#include "bigraph.h"
#include <iostream>
#include <queue>
#define LEFT false
#define RIGHT true
class Subgraph {
public:
	std::unordered_map<vid_t, std::vector<vid_t>> left_side_subgraph_neighbors;
	std::unordered_map<vid_t, std::vector<vid_t>> right_side_subgraph_neighbors;
	std::unordered_map<vid_t, int> left2score;
	std::unordered_map<vid_t, int> right2score;
	std::unordered_map<vid_t, bool> left_remain;
	std::unordered_map<vid_t, bool> right_remain;
	
public:
	bool selfcheck();
	Subgraph();
	~Subgraph() {}
};

extern void dynamic_KKCore_incremental_basic(BiGraph& g, vid_t u, vid_t v);

// naive increamental without +1
extern void rearrange_bicore_index_addition(std::vector<std::vector<bicore_index_block*>>& bicore_index, int fat, int oldsat, int newsat, vid_t target_node);

extern void update_index_with_fixed_left_k_addition_without_limit(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	int alpha, vid_t u, vid_t v);

extern void update_index_with_fixed_left_k_deletion_without_limit(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	int alpha, vid_t u, vid_t v);

extern void update_index_with_fixed_left_k_addition_with_limit(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	int alpha, vid_t u, vid_t v);

extern void update_index_with_fixed_left_k_deletion_with_limit(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	int alpha, vid_t u, vid_t v);

extern void update_bicore_index_without_limit(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v, 
	vid_t u, vid_t v, bool addition);

extern void update_bicore_index_with_limit(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v, 
	vid_t u, vid_t v, bool addition);

extern void update_bicore_index_without_limit_swap(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	vid_t u, vid_t v, bool addition);

extern void update_bicore_index_with_limit_swap(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	vid_t u, vid_t v, bool addition);

#endif // !DYNAMIC_H
