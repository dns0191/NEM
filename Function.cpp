#include "Function.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <Eigen/Dense>

double initializeNodesFromInput(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open input file.");
    }

    int DIM = 1;
    int GROUP_NUM = 1;
    double error = 0.0;
    std::vector<int> BENCH;
    std::vector<double> nodeWidths;
    std::vector<double> avgFluxValues;
    std::string line;
    std::string currentKey;

    size_t x_size = 0;
    size_t y_size = 0;

    while (std::getline(file, line)) {
        if (line.empty()) {
            currentKey.clear();
            continue;
        }

        std::istringstream iss(line);
        if (currentKey.empty()) {
            iss >> currentKey;
        }

        if (currentKey == "DIM") {
            iss >> DIM;
        }
        else if (currentKey == "GROUP_NUM") {
            iss >> GROUP_NUM;
        }
        else if (currentKey == "XS") {
            int region_id;
            iss >> region_id;
            std::vector<std::vector<double>> xs_values(GROUP_NUM, std::vector<double>(4, 0.0));
            for (int i = 0; i < 4; ++i) {
                for (int g = 0; g < GROUP_NUM; ++g) {
                    iss >> xs_values[g][i];
                }
            }
            crossSections[region_id] = xs_values;
        }
        else if (currentKey == "NODE_WIDTH") {
            double width;
            while (iss >> width) {
                nodeWidths.push_back(width);
            }
        }
        else if (currentKey == "BENCH") {
            int id;
            size_t current_row_size = 0;
            while (iss >> id) {
                BENCH.push_back(id);
                current_row_size++;
            }
            if (current_row_size > 0) {
                if (y_size == 0) {
                    y_size = current_row_size;
                }
                else if (y_size != current_row_size) {
                    throw std::runtime_error("Error: Inconsistent row sizes in BENCH data.");
                }
                x_size++;
            }
        }
        else if (currentKey == "AVG_FLUX") {
            double flux;
            while (iss >> flux) {
                avgFluxValues.push_back(flux);
            }
        }
        else if (currentKey == "CONVERGENCE")
        {
            iss >> error;
        }
    }
    file.close();

    if (DIM == 2) {
        if (x_size == 0 || y_size == 0) {
            throw std::runtime_error("Error: Invalid grid size (x_size or y_size is zero).");
        }
        if (BENCH.size() != x_size * y_size) {
            throw std::runtime_error("Error: BENCH size mismatch. Expected " +
                std::to_string(x_size * y_size) + ", got " + std::to_string(BENCH.size()));
        }
        if (nodeWidths.size() < x_size) {
            throw std::runtime_error("Error: Insufficient nodeWidths size. Expected at least " +
                std::to_string(x_size) + ", got " + std::to_string(nodeWidths.size()));
        }

        nodeGrid2D.resize(x_size, std::vector<MultiGroupNode*>(y_size, nullptr));

        for (size_t x = 0; x < x_size; ++x) {
            for (size_t y = 0; y < y_size; ++y) {
                size_t index = x * y_size + y;
                int region = BENCH[index];

                if (region == 0) {
                    continue;
                }
                else {
                    try {
                        Eigen::VectorXd node_widths = Eigen::VectorXd::Map(nodeWidths.data(), DIM);
                        nodeGrid2D[x][y] = new MultiGroupNode(static_cast<int>(index), region, GROUP_NUM, 2, node_widths);
                        nodeGrid2D[x][y]->setFluxAvg(avgFluxValues);
                    }
                    catch (const std::exception& e) {
                        std::cerr << "Error initializing node (" << x << ", " << y << "): " << e.what() << "\n";
                        throw;
                    }
                }
            }
        }

        // 경계 조건 설정
        for (size_t x = 0; x < x_size; ++x) {
            for (size_t y = 0; y < y_size; ++y) {
                if (nodeGrid2D[x][y] != nullptr) {
                    // 기본적으로 모든 면을 반사 조건으로 설정
                    nodeGrid2D[x][y]->setBoundaryCondition(0, LEFT_SIDE, REFLECTIVE);
                    nodeGrid2D[x][y]->setBoundaryCondition(0, RIGHT_SIDE, REFLECTIVE);
                    nodeGrid2D[x][y]->setBoundaryCondition(1, LEFT_SIDE, REFLECTIVE);
                    nodeGrid2D[x][y]->setBoundaryCondition(1, RIGHT_SIDE, REFLECTIVE);

                    // BENCH에 0이 있는 경우 이웃한 면을 진공 조건으로 설정
                    if (x > 0 && BENCH[(x - 1) * y_size + y] == 0) {
                        nodeGrid2D[x][y]->setBoundaryCondition(0, LEFT_SIDE, VACUUM);
                    }
                    if (x < x_size - 1 && BENCH[(x + 1) * y_size + y] == 0) {
                        nodeGrid2D[x][y]->setBoundaryCondition(0, RIGHT_SIDE, VACUUM);
                    }
                    if (y > 0 && BENCH[x * y_size + (y - 1)] == 0) {
                        nodeGrid2D[x][y]->setBoundaryCondition(1, LEFT_SIDE, VACUUM);
                    }
                    if (y < y_size - 1 && BENCH[x * y_size + (y + 1)] == 0) {
                        nodeGrid2D[x][y]->setBoundaryCondition(1, RIGHT_SIDE, VACUUM);
                    }
                }
            }
        }
    }

    return error;
}




void debugPrintNodes() {
    std::cout << "==== Node Data Debugging ====" << "\n";
    if (!nodeGrid1D.empty()) {
        std::cout << " 1D Node Grid Detected (" << nodeGrid1D.size() << " nodes)\n";
        for (const auto& node : nodeGrid1D) {
            if (node == nullptr) continue;
            node->getNodeInformation();
        }
    }
    if (!nodeGrid2D.empty()) {
        std::cout << " 2D Node Grid Detected (" << nodeGrid2D.size() << " x " << nodeGrid2D[0].size() << ")\n";
        for (const auto& row : nodeGrid2D) {
            for (const auto& node : row) {
                if (node == nullptr) continue;
                node->getNodeInformation();
            }
        }
    }
    if (!nodeGrid3D.empty()) {
        std::cout << "Unimplemented" << "\n";
    }
    std::cout << "\n";
    std::cout << "==== Cross Sections Data ====" << "\n";
    for (const auto& entry : crossSections) {
        const int region_id = entry.first;
        if (region_id <= 0) continue;
        const auto& xs_values = entry.second;
        std::cout << "Region ID: " << region_id << "\n";
        for (size_t g = 0; g < xs_values.size(); ++g) {
            std::cout << "  Group " << g + 1 << ": ";
            for (const auto& value : xs_values[g]) {
                std::cout << value << " ";
            }
            std::cout << "\n";
        }
    }
    std::cout << "==== End of Debugging ====" << "\n";
}


bool totalConvergence(double ERROR)
{
    for (const auto& row : nodeGrid2D) {
        for (const auto& node : row) {
            if (node != nullptr && !node->checkConvergence(ERROR)) {
                return false;
            }
        }
    }
	std::cout << "Converged!" << "\n";
    return true;
}
