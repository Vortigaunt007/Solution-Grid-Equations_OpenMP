#include "MatrixCSR.h"


MatrixCSR::MatrixCSR()
{

}

MatrixCSR::MatrixCSR(std::vector<int> &AI_t, std::vector<int> &AJ_t, std::vector<double> &A_t)
{
    const int n = (int)A_t.size();
    for(int i = 0; i < n; i++)
        A.push_back(A_t[i]);

    const int m = AI_t.size();
    for(int i = 0; i < m; i++)
        AI.push_back(AI_t[i]);

    const int k = AJ_t.size();
    for(int i = 0; i < k; i++)
        AJ.push_back(AJ_t[i]);
}

void MatrixCSR::Matrix_Portrait(Grid D)
{
    std::vector<int> sort_list_of_neighbors;
    int size = 0;

    for(int k = 0; k < D.getNz(); k++)
        for(int j = 0; j < D.getNy(); j++)
            for(int i = 0; i < D.getNx(); i++) {
            //    std::cout << i << j << k << D.getIndex(i, j, k) << std::endl;
                sort_list_of_neighbors = D.Neighboring_Nodes(i, j, k); // list of neighboring nodes for a node with coordinates (i, j, k)
                size = AJ.size();
                AI.push_back(size);
                const int bound = sort_list_of_neighbors.size();
                for(int l = 0; l < bound; l++) {
                    int neighbor_index = sort_list_of_neighbors[l];
                    AJ.push_back(neighbor_index);
                    A.push_back(1.0); // random value
                }
            }
    size = AJ.size();
    AI.push_back(size);
}

void MatrixCSR::Fill_Matrix(double f(int, int))
{
    int diag_elem_index = 0;
    double elem = 0;
    double sum_of_row_elem = 0.0; // sum of row elements except diagonal


    int N = AI.size() - 1;

    for(int i = 0; i < N; i++) { // string number
        sum_of_row_elem = 0.0;
        for(int j_ = AI[i]; j_ < AI[i+1]; j_++) {
            int j = AJ[j_]; // column number
            if (i == j) { // diagonal element
                diag_elem_index = j_; // we memorize the index in matrix A so that then we push the desired value there
                A[j_] = 1.0; // random value
            }
            else {
                elem = f(i, j);
                A[j_] = elem;
                sum_of_row_elem += abs(elem);

            }

        }
        A[diag_elem_index] = InitDiagElement(sum_of_row_elem); // diagonal element
    }
/*    for(int i = 0; i < N; i++) { // string number
        for(int j_ = AI[i]; j_ < AI[i+1]; j_++) {
            int j = AJ[j_]; // column number
            std::cout << " i = " << i << " j = " << j << " " << A[j_] << std::endl;
        }

    }
*/
//std::cout << "aAAAAA " << std::endl;
    saveToFile("FillMatrixCSR.txt");
}


void inverse_Matrix(MatrixCSR &A, MatrixCSR &A_1)
{
    const int n = (int)A.AI.size();
    for(int i = 0; i < n; i++)
        A_1.AI.push_back(A.AI[i]);

    const int m = (int)A.AJ.size();
    for(int i = 0; i < m; i++)
        A_1.AJ.push_back(A.AI[i]);

    const int k = (int)A.A.size();
    for(int i = 0; i < k; i++)
        A_1.A.push_back(1.0 / A.A[i]); // diagonal matrix
}

double InitDiagElement(double a)
{
    return 1.5 * a;
}

double InitOffDiagElement(int i, int j)
{
    return cos(double(i) * double(j));
    //return sin(i + j + 1);
}

void MatrixCSR::saveToFile(std::string filename) const
{
    std::ofstream fout(filename);

    const int n = (int)A.size();
    for(int i = 0; i < n; i++)
        fout << A[i] << " ";
    fout << std::endl;

    const int m = (int)AI.size();
    for(int i = 0; i < m; i++)
        fout << AI[i] << " ";
    fout << std::endl;

    const int k = (int)AJ.size();
    for(int i = 0; i < k; i++)
        fout << AJ[i] << " ";
    fout << std::endl;

    fout.close();
}

void MatrixCSR::printMatrixCSR_A() const
{
    std::cout << "MatrixCSR A.A" << std::endl;

    const int n = (int)A.size();
    for(int i = 0; i < n; i++)
        std::cout << A[i] << " ";
    std::cout<< std::endl;
}

void MatrixCSR::printMatrixCSR_Size() const
{
    std::cout << "A.size = " << A.size() << std::endl;
    std::cout << "AI.size = " << AI.size() << std::endl;
    std::cout << "AJ.size = " << AJ.size() << std::endl;
}
