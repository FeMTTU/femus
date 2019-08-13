#ifndef __femus_enums_SolvertypeEnum_hpp__
#define __femus_enums_SolvertypeEnum_hpp__

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
    FGMRES,
    GMRES,
    LGMRES,
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
