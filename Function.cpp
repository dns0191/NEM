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

void initializer(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open input file.");
    }

    std::string line;
    std::vector<double> node_widths;
    std::vector<int> bench_ids;

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string section;
        iss >> section;
        if (section == "XS")
        {
            double xs;

        }
        else if (section == "NODE_WIDTH")
        {
            double width;
            while (iss >> width)
            {
                if (width != 0.0)
                {
                    node_widths.push_back(width);
                }
            }
        }
        else if (section == "BENCH")
        {
            int id;
            while (iss >> id)
            {
                bench_ids.push_back(id);
            }
        }
    }

    file.close();
}
