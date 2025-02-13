#pragma once
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include "class.h"

//  영역(region)별 그룹(group)별 단면적 저장 (unordered_map 사용)
std::unordered_map<int, std::vector<std::vector<double>>> crossSections;

//  노드 저장소 (1D, 2D, 3D)
std::vector<MultiGroupNode*> nodeGrid1D;
std::vector<std::vector<MultiGroupNode*>> nodeGrid2D;
std::vector<std::vector<std::vector<MultiGroupNode*>>> nodeGrid3D;

double* product_matrix(double** A, double* C, int ng);
void initializer(const std::string& filename);