#ifndef GRID_H
#define GRID_H

#include "CommandLineArg.h"

enum cell_type
{
    PRISM,
    HEX
};

class Grid
{
    int Nx, Ny, Nz;
    int k1, k2;

public:
    Grid();
    Grid(int Nx, int Ny, int Nz, int k1, int k2);
    Grid(std::string filename);

    int getK1() const;
    int getK2() const;
    int getNx() const;
    int getNy() const;
    int getNz() const;

    int getIndex(int i, int j, int k);
    cell_type TypeCell(int i, int j, int k); // hexahedron or prism
    std::vector <int> Neighboring_Nodes(int i, int j, int k);
};

#endif // GRID_H
