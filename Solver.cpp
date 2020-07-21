#include "Solver.h"


Solver::Solver(MatrixCSR &M_in, Vector &b_in): A(M_in), b(b_in), JacobiPre()
{

}

void Solver::InitJacobiPreconditioner()
{
    // Jacobi Preconditioner is diagonal matrix: JacobiPre[diag_elem] = A[diag_elem]


    const int n = (int)A.AI.size() - 1;
    for(int i = 0; i < n; i++) { // string number
        for(int j_ = A.AI[i]; j_ < A.AI[i+1]; j_++) {
            int j = A.AJ[j_]; // column number
            if (i == j) {// diagonal element
                JacobiPre.AI.push_back(JacobiPre.AJ.size());
                JacobiPre.AJ.push_back(j); // number of column of the matrix A diagonal
                JacobiPre.A.push_back(A.A[j_]); // diagonal element of A
            }
        }
    }

    JacobiPre.AI.push_back(JacobiPre.AJ.size());
}


Vector Solver::ConjugateGradient(int nt)
{    
    BasicOperation Basic_Operation;
    int N = b.getSize();
    Vector x_k(N), z_k(N), r_k(N), p_k(N), q_k(N), Ax_k(N);
    double alpha, beta, rho_k_1 = 1.0, rho_k = 1.0;
    int count_step; // k
    bool convergence;
    double time_startSolve;

    // initialization step

    omp_set_num_threads(nt); // setting the number of threads

    A.Fill_Matrix(InitOffDiagElement); // Fill A in MatrixCSR A
    b.Fill_Vector(InitVectorElement); // Fill right part

    InitJacobiPreconditioner(); // Fill A in MatrixCSR Jacobi Preconditioner
    MatrixCSR M;
    inverse_Matrix(JacobiPre, M); // M = JacobiPre^(-1)
    count_step = 1;
    convergence = false;


    // basic operation time counting

    long long count_mult = 0, count_summ = 0, count_dot = 0;
    double start_time, time_mult_Matrix_Vector, time_summ, time_dot;

    start_time = omp_get_wtime();
    for (int i = 0; i < 100; i++)
       // Basic_Operation.mult_Matrix_Vector(A, x_k, Ax_k);
    time_mult_Matrix_Vector = omp_get_wtime() - start_time;
    time_mult_Matrix_Vector /= 100.0;
    count_mult = Basic_Operation.countMult;
    count_mult /= 100;

    start_time = omp_get_wtime();
    for (int i = 0; i < 100; i++)
        Basic_Operation.summAlpha(1.0, b, -1.0, Ax_k, r_k);
    time_summ = omp_get_wtime() - start_time;
    count_summ = Basic_Operation.countSummAlpha;

    start_time = omp_get_wtime();
    for (int i = 0; i < 100; i++)
        rho_k = Basic_Operation.dotProduct(r_k, z_k);
    time_dot = omp_get_wtime() - start_time;
    count_dot = Basic_Operation.countDot;

    Basic_Operation.resetCounters();  // resetting counters: countDot, timeDot ...
//M.printMatrixCSR_A();
    // first step of the conjugate gradient method
/*
    Basic_Operation.mult_Matrix_Vector(A, x_k, Ax_k); // Ax_0 = 0
    Basic_Operation.summAlpha(1.0, b, -1.0, Ax_k, r_k); // r_0 = b - Ax_0 = b
    Basic_Operation.mult_Matrix_Vector(M, r_k, z_k); // z_0 = M^(-1)r_0
    rho_k = Basic_Operation.dotProduct(r_k, z_k); // rho_0 = (r_0, z_0)

    Basic_Operation.resetCounters();

    time_startSolve = omp_get_wtime();

    M.printMatrixCSR_A();

    p_k = z_k; // p_0 = z_0
  */
    Basic_Operation.summAlpha(1.0, b, -1.0, Ax_k, r_k); // r_0 = b - Ax_0 = b
   // Basic_Operation.mult_Matrix_Vector(M, r_k, z_k, K); // z_0 = M^(-1)r_0
   // Basic_Operation.update(A, z_k, K);
   // rho_k = Basic_Operation.dotProduct(r_k, z_k, A.size_Own); // rho_0 = (r_0, z_0)


    Basic_Operation.resetCounters();

   // Basic_Operation.mult_Matrix_Vector(A, p_k, q_k);
   // p_k = z_k; // p_0 = z_0

    int k = 1;
    double alpha_prev = 0.0;

    while (!convergence) {
        Basic_Operation.mult_Matrix_Vector(M, r_k, z_k);
     //   Basic_Operation.update(A, z_k, K);
        //M.printMatrixCSR_A();
      //  z_k.printVector();
        alpha = Basic_Operation.dotProduct(r_k, z_k);

       // std::cout << "alpha  = " << alpha << std::endl;
        if(k == 1) {
            p_k = z_k;
        } else {
            beta = alpha/alpha_prev;
            Basic_Operation.countSolver += 1;
            Basic_Operation.summAlpha(1.0, z_k, beta, p_k, p_k);
        }

     //   std::cout << "Size = " << beta << std::endl;//p_k.printVector();
       // p_k.printVector();
        Basic_Operation.mult_Matrix_Vector(A, p_k, q_k);
     //   Basic_Operation.update(A, q_k, K);
      //  A.printMatrixCSR_A();
q_k.printVector();
        double dottmp, gamma;
        dottmp = Basic_Operation.dotProduct(p_k, q_k);

        Basic_Operation.countSolver += 1;
        gamma = alpha / dottmp;

        Basic_Operation.summAlpha(1.0, x_k, gamma, p_k, x_k);
        Basic_Operation.summAlpha(1.0, r_k, -gamma, q_k, r_k);

     //   std::cout << "gamma  = " << dottmp << std::endl;
       // x_k.printVector();

        alpha_prev = alpha;
      //  rho_k = Basic_Operation.dotProduct(r_k, r_k);

        //if (world_rank == 1)
            std::cout << "Iteration number = " << k << " " << "r_k = " << alpha << std::endl;
        if ((alpha < eps) || (k >= maxiter))
            convergence = true;
        else
            k++;
    }

/*
    while (!convergence) {
        Basic_Operation.mult_Matrix_Vector(A, p_k, q_k);
        alpha = rho_k / Basic_Operation.dotProduct(p_k, q_k);
        Basic_Operation.countSolver += 1;
        Basic_Operation.summAlpha(1.0, x_k, alpha, p_k, x_k);
        Basic_Operation.summAlpha(1.0, r_k, -alpha, q_k, r_k);

        std::cout << "Iteration number = " << count_step << " " << "r_k = " << sqrt(Basic_Operation.dotProduct(r_k, r_k)) << std::endl;
        if ((rho_k < eps) || (count_step >= maxiter))
            convergence = true;
        else
            count_step++;
        rho_k_1 = rho_k;
        Basic_Operation.mult_Matrix_Vector(M, r_k, z_k);
        rho_k = Basic_Operation.dotProduct(r_k, z_k);
        beta = rho_k / rho_k_1;
        Basic_Operation.countSolver += 1;
        Basic_Operation.summAlpha(1.0, z_k, beta, p_k, p_k);
    }
*/
    timeSolveCG = omp_get_wtime() - time_startSolve;

    // time is considered inside the solver cycle

    std::cout << std::endl;
    std::cout << "Сycle time counting " << std::endl;
    std::cout << "GFlops mult_Matrix_Vector = " << count_mult/(time_mult_Matrix_Vector * 1E9) << std::endl;
    std::cout << "GFlops summAlpha = " << count_summ/(time_summ * 1E9) << std::endl;
    std::cout << "GFlops dotProduct = " << count_dot/(time_dot * 1E9) << std::endl;
    std::cout << "GFlops Solver = " << Basic_Operation.countSolver/(timeSolveCG * 1E9) << std::endl;


    // time is considered inside functions
    std::cout << std::endl;
    std::cout << "Solver time counting " << std::endl;
    std::cout << "mult_Matrix_Vector GFlops = " << Basic_Operation.countMult/(Basic_Operation.timeMult * 1E9) << std::endl;
    std::cout << "summAlpha GFlops = " << Basic_Operation.countSummAlpha/(Basic_Operation.timeSummAlpha * 1E9) << std::endl;
    std::cout << "dotProduct GFlops = " << Basic_Operation.countDot/(Basic_Operation.timeDot * 1E9) << std::endl;
    std::cout << "Solver GFlops = " << Basic_Operation.countSolver/(timeSolveCG * 1E9) << std::endl;
    std::cout << std::endl;
double c = Basic_Operation.dotProduct(x_k, x_k);
std::cout << "AAAANSWER " << c << std::endl;
x_k.printVector();
    return x_k;
}

/*
void Solver::saveToFileData(std::string filename) const
{
    std::ofstream fout(filename);

    fout << "A.A.size() = " << A.A.size() << std::endl;
    fout << "A.AI.size() = " << A.AI.size() << std::endl;
    fout << "A.AJ.size() = " << A.AJ.size() << std::endl;

    fout << "time_spent_filling = " << time_spent_filling << std::endl;
    fout << "time_spent_solving = " << time_spent_solving << std::endl;
    fout << "time_dotProduct = " << time_dotProduct << std::endl;
    fout << "time_summAlpha = " << time_summAlpha << std::endl;
    fout << "time_mult_Matrix_Vector = " << time_mult_Matrix_Vector << std::endl;

    fout.close();
}
*/
