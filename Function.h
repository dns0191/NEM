#pragma once
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include "class.h"



double* product_matrix(double** A, double* C, int ng);
double initializeNodesFromInput(const std::string& filename);
void debugPrintNodes();
bool totalConvergence(double ERROR);