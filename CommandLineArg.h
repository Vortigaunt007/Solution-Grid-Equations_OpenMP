#ifndef COMMANDLINEARG_H
#define COMMANDLINEARG_H

#include <iostream>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <omp.h>

class CommandLineArg
{
public:
    int nx, ny, nz;
    int k1, k2;
    int nt;

    CommandLineArg(int, int, int, int, int, int);
    CommandLineArg(std::string filename);
};

#endif // COMMANDLINEARG_H
