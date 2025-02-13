#include "Function.h"

double* product_matrix(double** A, double* C, int ng)
{
    double* result = new double[ng];

    for (int i = 0; i < ng; i++)
    {
        result[i] = 0.0;
        for (int j = 0; j < ng; j++)
        {
            result[i] += A[i][j] * C[j];
        }
    }
    return result;
}


void initializeNodesFromInput(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open input file.");
    }

    int DIM = 1;
    int GROUP_NUM = 1;
    std::vector<int> BENCH;
    std::vector<double> nodeWidths;
    bool readingXS = false;
    bool readingNodeWidth = false;
    bool readingBench = false;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "DIM") {
            iss >> DIM;
        }
        else if (key == "GROUP_NUM") {
            iss >> GROUP_NUM;
        }
        else if (key == "NODE_WIDTH") {
            readingNodeWidth = true;
        }
        else if (key == "BENCH") {
            readingBench = true;
        }
        else if (key == "XS") {
            readingXS = true;
        }
        else if (readingNodeWidth) {
            double width;
            while (iss >> width) {
                nodeWidths.push_back(width);
            }
            readingNodeWidth = false; // 읽기 완료 후 false로 설정
        }
        else if (readingBench) {
            int id;
            while (iss >> id) {
                BENCH.push_back(id);
            }
            readingBench = false; // 읽기 완료 후 false로 설정
        }
        else if (readingXS && !key.empty()) {
            int region_id = std::stoi(key);
            std::vector<std::vector<double>> xs_values(GROUP_NUM, std::vector<double>(4, 0.0));
            for (int i = 0; i < 4; ++i) {
                for (int g = 0; g < GROUP_NUM; ++g) {
                    iss >> xs_values[g][i];
                }
            }
            crossSections[region_id] = xs_values;
        }
    }
    file.close();

    // ✅ 차원별 노드 생성 (오버플로 방지 추가)
    if (DIM == 1) {
        nodeGrid1D.resize(BENCH.size(), nullptr);
        for (size_t i = 0; i < BENCH.size(); ++i) {
            int region = BENCH[i];
            nodeGrid1D[i] = new MultiGroupNode(static_cast<int>(i), region, GROUP_NUM, 1, nodeWidths.data());
        }
    }
    else if (DIM == 2) {
        if (nodeWidths.empty()) {
            throw std::runtime_error("Error: NODE_WIDTH is missing in input file.");
        }
        size_t x_size = nodeWidths.size();
        if (x_size == 0) {
            throw std::runtime_error("Error: x_size is 0. Possible NODE_WIDTH reading error.");
        }
        if (BENCH.size() % x_size != 0) {
            throw std::runtime_error("Error: BENCH size is not a multiple of x_size. Check input file.");
        }
        size_t y_size = BENCH.size() / x_size;
        nodeGrid2D.resize(x_size, std::vector<MultiGroupNode*>(y_size, nullptr));

        for (size_t x = 0; x < x_size; ++x) {
            for (size_t y = 0; y < y_size; ++y) {
                size_t index = x * y_size + y;
                int region = BENCH[index];
                nodeGrid2D[x][y] = new MultiGroupNode(static_cast<int>(index), region, GROUP_NUM, 2, nodeWidths.data());
            }
        }
    }
    else if (DIM == 3) {
        if (nodeWidths.empty()) {
            throw std::runtime_error("Error: NODE_WIDTH is missing in input file.");
        }
        size_t x_size = nodeWidths.size();
        size_t y_size = nodeWidths.size();
        if (x_size == 0 || y_size == 0) {
            throw std::runtime_error("Error: x_size or y_size is 0. Check NODE_WIDTH in input file.");
        }
        if (BENCH.size() % (x_size * y_size) != 0) {
            throw std::runtime_error("Error: BENCH size is not a multiple of x_size * y_size.");
        }
        size_t z_size = BENCH.size() / (x_size * y_size);

        nodeGrid3D.resize(x_size, std::vector<std::vector<MultiGroupNode*>>(
            y_size, std::vector<MultiGroupNode*>(z_size, nullptr)));

        for (size_t x = 0; x < x_size; ++x) {
            for (size_t y = 0; y < y_size; ++y) {
                for (size_t z = 0; z < z_size; ++z) {
                    size_t index = x * y_size * z_size + y * z_size + z;
                    int region = BENCH[index];
                    nodeGrid3D[x][y][z] = new MultiGroupNode(static_cast<int>(index), region, GROUP_NUM, 3, nodeWidths.data());
                }
            }
        }
    }
}





// ✅ 초기화된 노드 데이터 확인
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
        for (size_t x = 0; x < nodeGrid2D.size(); ++x) {
            for (size_t y = 0; y < nodeGrid2D[x].size(); ++y) {
                const MultiGroupNode* node = nodeGrid2D[x][y];
                if (node == nullptr) continue;
                node->getNodeInformation();
            }
        }
    }

    if (!nodeGrid3D.empty()) {
        std::cout << " 3D Node Grid Detected (" << nodeGrid3D.size() << " x " << nodeGrid3D[0].size() << " x " << nodeGrid3D[0][0].size() << ")\n";
        for (size_t x = 0; x < nodeGrid3D.size(); ++x) {
            for (size_t y = 0; y < nodeGrid3D[x].size(); ++y) {
                for (size_t z = 0; z < nodeGrid3D[x][y].size(); ++z) {
                    const MultiGroupNode* node = nodeGrid3D[x][y][z];
                    if (node == nullptr) continue;
                    node->getNodeInformation();
                }
            }
        }
    }

    std::cout << "==== Cross Sections Data ====" << "\n";
    for (const auto& entry : crossSections) {
        int region_id = entry.first;
        const auto& xs_values = entry.second;
        std::cout << "Region ID: " << region_id << "\n";
        for (size_t g = 0; g < xs_values.size(); ++g) {
            std::cout << "  Group " << g+1 << ": ";
            for (size_t i = 0; i < xs_values[g].size(); ++i) {
                std::cout << xs_values[g][i] << " ";
            }
            std::cout << "\n";
        }
    }

    std::cout << "==== End of Debugging ====" << "\n";
}

