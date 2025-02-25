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
    const double error = initializeNodesFromInput("input.txt");
    std::cout << "nodeGrid1D size: " << nodeGrid1D.size() << "\n";
    std::cout << "nodeGrid2D size: " << nodeGrid2D.size() << ", " << nodeGrid2D[0].size() << "\n";

    debugPrintNodes();
    std::ofstream outputFile("output.out");
    if (!outputFile.is_open()) {
        std::cerr << "Failed to open output.out file.\n";
        return 1;
    }
    if (outputFile.is_open()) {
        int step = 0;
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
            double* max_flux_avg = new double[2] {0.0, 0.0};

            // 모든 노드의 flux_avg 중 최대값 찾기
            for (const auto& row : nodeGrid2D) {
                for (const auto& node : row) {
					for (int i = 0; i < node->getNumberOfGroups(); i++) {
						if (node->getFlux(i) > max_flux_avg[i])
							max_flux_avg[i] =  node->getFlux(i);
					}
                }
            }

            // 모든 노드의 flux_avg 값을 최대값으로 나누어 정규화
            for (auto& row : nodeGrid2D) {
                for (MultiGroupNode* node : row) {
					for (int i = 0; i < node->getNumberOfGroups(); i++) {
						node->normalizeFluxAvg(max_flux_avg[i]);
					}
                }
            }

            for (const auto& row : nodeGrid2D) {
                for (const auto& node : row) {
                    if (node != nullptr) {
                        outputFile << std::setw(15) << std::scientific << node->getFlux(1);
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
