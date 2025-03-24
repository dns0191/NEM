#include <iostream>
#include <fstream>
#include <iomanip>
#include "class.h"
#include "Function.h"

int main() {
    // debug.out 파일 초기화
    std::ofstream debugFile("debug.txt", std::ios_base::trunc);
    if (!debugFile.is_open()) {
        return 1;
    }
    debugFile.close();
    const double error = initializeNodesFromInput("input.txt");
    std::cout << "nodeGrid1D size: " << nodeGrid1D.size() << "\n";
    std::cout << "nodeGrid2D size: " << nodeGrid2D.size() << ", " << nodeGrid2D[0].size() << "\n";

    debugPrintNodes();
    std::ofstream outputFile("output.txt");
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
                else {
                    continue;
                }
            }
        }
        outputFile << "\n";

        while (!totalConvergence(error)) {
            double* max_flux_avg = new double[2] {0.0, 0.0};
            double* max_out_current = new double[2] {0.0, 0.0};

            std::vector<std::vector<Eigen::VectorXd>> new_flux_avg(nodeGrid2D.size(), std::vector<Eigen::VectorXd>(nodeGrid2D[0].size()));
            std::vector<std::vector<std::vector<Eigen::MatrixXd>>> new_out_current(nodeGrid2D.size(), std::vector<std::vector<Eigen::MatrixXd>>(nodeGrid2D[0].size()));

            for (auto& row : nodeGrid2D) {
                for (auto& node : row) {
                    if (node != nullptr && node->getId() % 2 == 0)
                        node->runNEM();
                }
                for (auto& node : row) {
                    if (node != nullptr && node->getId() % 2 != 0)
                        node->runNEM();
                }
            }

            outputFile << std::setw(5) << step;
            for (auto& row : nodeGrid2D) {
                for (auto& node : row) {
                    if (node != nullptr) {
                        outputFile << std::setw(15) << std::scientific << node->getFlux(0);
						node->getNodeInformationRaw();
                    }
                }
            }
            outputFile << "\n";
            step ++;
        }
        outputFile.close();
    }
    return 0;
}
