#ifndef SOLVER_H
#define SOLVER_H

#include "Vector.h"
#include "MatrixCSR.h"
#include "Grid.h"
#include "BasicOperation.h"

class Solver
{
    MatrixCSR &A;
    Vector &b;
    MatrixCSR JacobiPre;

    const int maxiter = (A.AI.size() - 1)*log((A.AI.size() - 1)); // N * ln(N) * const
    const double eps = 0.0000000001;

    double timeSolveCG;
public:
    Solver(MatrixCSR &, Vector &);

    void InitJacobiPreconditioner();

    Vector ConjugateGradient(int nt);

    void saveToFileData(std::string filename) const;
};

#endif // SOLVER_H
