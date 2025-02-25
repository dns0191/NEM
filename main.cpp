#include <iostream>
#include <fstream>
#include <iomanip>
#include "class.h"
#include "Function.h"


int main() {
    // debug.out 파일 초기화
    std::ofstream debugFile("debug.out", std::ios_base::trunc);
    if (!debugFile.is_open()) {
        std::cerr << "Unable to open debug.out file.\n";
        return 1;
    }
    debugFile.close();
    const double error = initializeNodesFromInput("input.inp");
    int step = 0;
    debugPrintNodes();
    std::ofstream outputFile("output.out");
    if (outputFile.is_open()) {
        outputFile << "     ";
        for (const auto& row : nodeGrid2D) {
            for (const auto& node : row) {
                if (node != nullptr) {
                    outputFile << std::setw(15) << " " << std::setw(3) << node->getId() << std::setw(12) << " ";
                }
            }
        }
        outputFile << "\n";

        while (!totalConvergence(error)) {
            outputFile << std::setw(5) << step;
            for (const auto& row : nodeGrid2D) {
                for (const auto& node : row) {
                    if (node != nullptr) {
                        node->runNEM();
                        node->getNodeInformation();
                    	outputFile << std::setw(15) << std::scientific << node->getFlux(1);
                        
                    }
                }
            }

            // 그룹마다 전체 flux_avg의 최대값 찾기
            int numberOfGroups = nodeGrid2D[0][0]->getNumberOfGroups();
            double* max_flux_avg = new double[numberOfGroups];
            std::fill(max_flux_avg, max_flux_avg + numberOfGroups, 0.0);

            for (const auto& row : nodeGrid2D) {
                for (const auto& node : row) {
                    if (node != nullptr) {
                        for (int i = 0; i < numberOfGroups; i++) {
							if (abs(node->getFlux(i)) > max_flux_avg[i])
							{
                                max_flux_avg[i] = node->getFlux(i);
							}
                        }
                    }
                }
            }

            // 각 노드의 flux_avg 값을 최대값으로 나누어 정규화
            for (auto& row : nodeGrid2D) {
                for (MultiGroupNode* node : row) {
                    if (node != nullptr) {
                        for (int i = 0; i < numberOfGroups; i++) {
                            node->normalizeFluxAvg(max_flux_avg[i], i);
                        }
                    }
                }
            }

            delete[] max_flux_avg;
            outputFile << "\n";
            step += 1;
        }
        outputFile.close();
    }
    return 0;
}
