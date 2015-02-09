#ifndef __solvertype__
#define __solvertype__

enum SolverType {
    CG=0,
    CGN,
    CGS,
    CR,
    QMR,
    TCQMR,
    TFQMR,
    BICG,
    BICGSTAB,
    MINRES,
    GMRES,
    LSQR,
    JACOBI,
    SOR_FORWARD,
    SOR_BACKWARD,
    SSOR,
    RICHARDSON,
    CHEBYSHEV,
    LUMP,
    INVALID_SOLVER,
    PREONLY
};

#endif
