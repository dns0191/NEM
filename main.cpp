#include <iostream>
#include <fstream>
#include <iomanip>
#include "class.h"
#include "Function.h"


int main() {
    // debug.out 파일 초기화
    std::ofstream debugFile("debug.out", std::ios_base::trunc);
    if (!debugFile.is_open()) {
        return 1;
    }
    debugFile.close();
    const double error = initializeNodesFromInput("input.inp");
    std::cout << "nodeGrid1D size: " << nodeGrid1D.size() << std::endl;
    std::cout << "nodeGrid2D size: " << nodeGrid2D.size() << ", " << nodeGrid2D[0].size() << std::endl;

    int step = 0;
    debugPrintNodes();
    std::ofstream outputFile("output.out");
    if (outputFile.is_open()) {
        outputFile << "     ";
        for (const auto& row : nodeGrid2D) {
            for (const auto& node : row) {
                if (node != nullptr) {
                    outputFile << std::setw(12) << " " << std::setw(3) << node->getId();
                }
            }
        }
        outputFile << "\n";

        while (!totalConvergence(error)) {
            outputFile << std::setw(5) << step;
            for (const auto& row : nodeGrid2D) {
                for (const auto& node : row) {
                    if (node != nullptr) {
                        outputFile << std::setw(15) << std::scientific << node->getCurrent(1);
                        node->runNEM();
                        node->getNodeInformation();
                    }
                }
            }
            outputFile << "\n";
            step += 1;
        }
        outputFile.close();
    }
    return 0;
}
