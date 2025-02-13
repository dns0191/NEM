#pragma once
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include "class.h"



double* product_matrix(double** A, double* C, int ng);
void initializeNodesFromInput(const std::string& filename);
void debugPrintNodes();