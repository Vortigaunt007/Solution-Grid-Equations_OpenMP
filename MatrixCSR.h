#ifndef MATRIXCSR_H
#define MATRIXCSR_H

#include "Grid.h"

class MatrixCSR
{
public:
    std::vector<int> AI;
    std::vector<int> AJ;
    std::vector<double> A;

    MatrixCSR();
    MatrixCSR(std::vector<int> &AI, std::vector<int> &AJ, std::vector<double> &A);

    void Matrix_Portrait(Grid);
    void Fill_Matrix(double f(int, int));

    void saveToFile(std::string filename) const;
    void printMatrixCSR_Size() const;
    void printMatrixCSR_A() const;
};

void inverse_Matrix(MatrixCSR &A, MatrixCSR &A_1);
double InitDiagElement(double);
double InitOffDiagElement(int i, int j);

#endif // MATRIXCSR_H
