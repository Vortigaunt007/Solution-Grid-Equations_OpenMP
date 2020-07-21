#include "CommandLineArg.h"
#include "Grid.h"
#include "MatrixCSR.h"
#include "Vector.h"
#include "Solver.h"

int main()
{
    CommandLineArg arg("input.txt"); // Nx, Ny, Nz, k1, k2, nt

    Grid D(arg.nx, arg.ny, arg.nz, arg.k1, arg.k2);

    MatrixCSR M;
    M.Matrix_Portrait(D);
    Vector b(M.AI.size() - 1);

    Solver S(M, b);
    S.ConjugateGradient(arg.nt);

    return 0;
}
