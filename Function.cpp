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

    int DIM = 1;  // 기본 차원
    int GROUP_NUM = 1;
    std::vector<int> BENCH;
    std::vector<double> nodeWidths;
    bool readingXS = false;

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
            double width;
            while (iss >> width) {
                nodeWidths.push_back(width);
            }
        }
        else if (key == "BENCH") {
            int id;
            while (iss >> id) {
                BENCH.push_back(id);
            }
        }
        else if (key == "XS") {
            readingXS = true;
        }
        else if (readingXS && !key.empty()) {
            int region_id = std::stoi(key);
            std::vector<std::vector<double>> xs_values(GROUP_NUM, std::vector<double>(8, 0.0));  // 기본 0으로 초기화
            for (int g = 0; g < GROUP_NUM; ++g) {
                for (int i = 0; i < 8; ++i) {
                    iss >> xs_values[g][i];
                }
            }
            crossSections[region_id] = xs_values;
        }
    }
    file.close();

    // ✅ 차원에 따라 노드 저장소 초기화
    if (DIM == 1) {
        nodeGrid1D.resize(BENCH.size(), nullptr);
        for (size_t i = 0; i < BENCH.size(); ++i) {
            int region = BENCH[i];  // 노드의 영역 ID
            nodeGrid1D[i] = new MultiGroupNode(i, region, GROUP_NUM, 1);
        }
    }
    else if (DIM == 2) {
        int x_size = nodeWidths.size();
        int y_size = BENCH.size() / x_size;

        nodeGrid2D.resize(x_size, std::vector<MultiGroupNode*>(y_size, nullptr));
        for (int x = 0; x < x_size; ++x) {
            for (int y = 0; y < y_size; ++y) {
                int region = BENCH[x * y_size + y];  // 노드의 영역 ID
                nodeGrid2D[x][y] = new MultiGroupNode(x * y_size + y, region, GROUP_NUM, 2);
            }
        }
    }
    else if (DIM == 3) {
        int x_size = nodeWidths.size();
        int y_size = nodeWidths.size();
        int z_size = BENCH.size() / (x_size * y_size);

        nodeGrid3D.resize(x_size, std::vector<std::vector<MultiGroupNode*>>(
            y_size, std::vector<MultiGroupNode*>(z_size, nullptr)));

        for (int x = 0; x < x_size; ++x) {
            for (int y = 0; y < y_size; ++y) {
                for (int z = 0; z < z_size; ++z) {
                    int region = BENCH[x * y_size * z_size + y * z_size + z];  // 노드의 영역 ID
                    nodeGrid3D[x][y][z] = new MultiGroupNode(
                        x * y_size * z_size + y * z_size + z, region, GROUP_NUM, 3);
                }
            }
        }
    }
}
