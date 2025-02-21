#pragma once
#include "class.h"

double initializeNodesFromInput(const std::string& filename);
void debugPrintNodes();
bool totalConvergence(double ERROR);