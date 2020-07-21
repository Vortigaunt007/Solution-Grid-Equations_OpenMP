#ifndef BASICOPERATION_H
#define BASICOPERATION_H

#include "CommandLineArg.h"
#include "Vector.h"
#include "MatrixCSR.h"

class BasicOperation
{
    int nt = 1;
public:
    long long  countSummAlpha = 0, countMult = 0, countDot = 0, countSolver = 0;
    long long  countAISummAlpha = 0, countAIMult = 0, countAIDot = 0;
    long long  countAISolver = 0;
    double timeSummAlpha = 0, timeMult = 0, timeDot = 0;

    BasicOperation();

    bool checkSizes(Vector &v1, Vector &v2);

    void resetCounters();

    double dotProduct(Vector &v1, Vector &v2);
    void summAlpha(double alpha, Vector &v1_in, double beta, Vector &v2_in, Vector &out);
    void mult_Matrix_Vector(MatrixCSR &, Vector &, Vector &);
};

#endif // BASICOPERATION_H
