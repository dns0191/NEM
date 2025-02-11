#pragma once
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <vector>
double* product_matrix(double** A, double* C, int ng);
void initializer(const std::string& filename);