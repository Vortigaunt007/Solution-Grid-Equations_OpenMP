#include "BasicOperation.h"

BasicOperation::BasicOperation()
{

}

bool BasicOperation::checkSizes(Vector &v1, Vector &v2)
{
    if(v1.getSize() != v2.getSize()) {
        std::cout << "Sizes aren't equale!" << std::endl;
        return false;
    } else
        return true;
}

void BasicOperation::resetCounters()
{
    countDot = 0;
    countMult = 0;
    countSolver = 0;
    countSummAlpha = 0;

    timeDot = 0;
    timeMult = 0;
    timeSummAlpha = 0;
}

void BasicOperation::mult_Matrix_Vector(MatrixCSR &M, Vector &in, Vector &out)
{

    const int n = (int)M.AI.size() - 1;
    if(n != in.getSize()) {
        std::cout << "Function mult_Matrix_Vector: sizes aren't equale" << std::endl;
        return;
    }

    double t = omp_get_wtime();
    long long count = 0;

    #pragma omp parallel for reduction(+:count)
    for(int i = 0; i < in.getSize(); i++) { // string number
        double res = 0.0;
        const int j_begin = M.AI[i];
        const int j_end = M.AI[i+1];
        for(int j_ = j_begin; j_ < j_end; j_++) {
            int j = M.AJ[j_]; // column number

           // std::cout << " i = " << i << " j = " << j  << std::endl;
            double a_ij = M.A[j_];
            count += 2;
            res += a_ij * in[j];
        }
        out[i] = res;
    }

    timeMult += omp_get_wtime() - t;
    countAIMult += 5 * (count / 2);
    countAISolver += 5 * (count / 2);
    countMult += count;
    countSolver += count;
}

void BasicOperation::summAlpha(double alpha, Vector &v1_in, double beta, Vector &v2_in, Vector &out)
{
    if(!checkSizes(v1_in, v2_in))
        return;

    double t = omp_get_wtime();

    #pragma omp parallel for
    for(int i = 0; i < v1_in.getSize(); i++) {
        out[i] = alpha * v1_in[i] + beta * v2_in[i];
    }

    timeSummAlpha += omp_get_wtime() - t;
    countAISummAlpha += 5 * v1_in.getSize();
    countAISolver += 5 * v1_in.getSize();
    countSolver += 3 * v1_in.getSize();
    countSummAlpha += 3 * v1_in.getSize();
}

double BasicOperation::dotProduct(Vector &v1, Vector &v2)
{
    if(!checkSizes(v1, v2))
        return 0.0;

    double result = 0.0;
    double t = omp_get_wtime();

    #pragma omp parallel for reduction(+:result)
    for(int i = 0; i < v1.getSize(); i++) {
        result += v1[i] * v2[i];
    }

    timeDot += omp_get_wtime() - t;
    countAIDot += 4 * v1.getSize();
    countAISolver += 4 * v1.getSize();
    countSolver += 2 * v1.getSize();
    countDot += 2 * v1.getSize();

    return result;
}
