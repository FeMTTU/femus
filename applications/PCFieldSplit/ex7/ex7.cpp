static char help[] = "Stokes Problem with Temperature in 2d and 3d with simplicial finite elements.\n\
We solve the Stokes problem in a rectangular\n\
domain, using a parallel unstructured mesh (DMPLEX) to discretize it.\n\n\n";

/*
TODO for Mantle Convection:
#NAME?
#NAME?
#NAME?

The isoviscous Stokes problem, which we discretize using the finite
element method on an unstructured mesh. The weak form equations are

< \nabla v, \nabla u + {\nabla u}^T > - < \nabla\cdot v, p > + < v, f > = 0
< q, \nabla\cdot v >                                                    = 0
< \nabla t, \nabla T>                                                   = q_T

Boundary Conditions:

#NAME?

#NAME?
#NAME?

Discretization:

We use a Python script to generate a tabulation of the finite element basis
functions at quadrature points, which we put in a C header file. The generic
command would be:

bin/pythonscripts/PetscGenerateFEMQuadrature.py dim order dim 1 laplacian dim order 1 1 gradient src/snes/examples/tutorials/ex62.h

We can currently generate an arbitrary order Lagrange element. The underlying
FIAT code is capable of handling more exotic elements, but these have not been
tested with this code.

Field Data:

Sieve data is organized by point, and the closure operation just stacks up the
data from each sieve point in the closure. Thus, for a P_2-P_1 Stokes element, we
have

cl{e} = {f e_0 e_1 e_2 v_0 v_1 v_2}
x     = [u_{e_0} v_{e_0} u_{e_1} v_{e_1} u_{e_2} v_{e_2} u_{v_0} v_{v_0} p_{v_0} u_{v_1} v_{v_1} p_{v_1} u_{v_2} v_{v_2} p_{v_2}]

The problem here is that we would like to loop over each field separately for
integration. Therefore, the closure visitor in DMPlexVecGetClosure() reorders
the data so that each field is contiguous

x'    = [u_{e_0} v_{e_0} u_{e_1} v_{e_1} u_{e_2} v_{e_2} u_{v_0} v_{v_0} u_{v_1} v_{v_1} u_{v_2} v_{v_2} p_{v_0} p_{v_1} p_{v_2}]

Likewise, DMPlexVecSetClosure() takes data partitioned by field, and correctly
puts it into the Sieve ordering.
*/

#include <petscdmplex.h>
#include <petscsnes.h>

int main(){
  return 1;
}

/*------------------------------------------------------------------------------
This code can be generated using 'bin/pythonscripts/PetscGenerateFEMQuadrature.py dim order dim 1 laplacian dim order 1 1 gradient dim order 1 1 identity src/snes/examples/tutorials/ex31.h'
-----------------------------------------------------------------------------*/
// #include "ex31.h"
//
// typedef enum {DIRICHLET, FREE_SLIP} BCType;
// typedef enum {RUN_FULL, RUN_TEST} RunType;
// typedef enum {FORCING_CONSTANT, FORCING_LINEAR, FORCING_CUBIC} ForcingType;
//
// typedef struct {
// DM            dm;                /* REQUIRED in order to use SNES evaluation functions */
// PetscFEM      fem;               /* REQUIRED to use DMPlexComputeResidualFEM() */
// PetscInt      debug;             /* The debugging level */
// PetscMPIInt   rank;              /* The process rank */
// PetscMPIInt   numProcs;          /* The number of processes */
// RunType       runType;           /* Whether to run tests, or solve the full problem */
// PetscBool     jacobianMF;        /* Whether to calculate the Jacobian action on the fly */
// PetscLogEvent createMeshEvent;
// PetscBool     showInitial, showSolution;
// /* Domain and mesh definition */
// PetscInt      dim;               /* The topological mesh dimension */
// PetscBool     interpolate;       /* Generate intermediate mesh elements */
// PetscReal     refinementLimit;   /* The largest allowable cell volume */
// char          partitioner[2048]; /* The graph partitioner */
// /* GPU partitioning */
// PetscInt      numBatches;        /* The number of cell batches per kernel */
// PetscInt      numBlocks;         /* The number of concurrent blocks per kernel */
// /* Element quadrature */
// PetscQuadrature q[NUM_FIELDS];
// /* Problem definition */
// void (*f0Funcs[NUM_FIELDS])(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f0[]); /* f0_u(x,y,z), and f0_p(x,y,z) */
// void (*f1Funcs[NUM_FIELDS])(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f1[]); /* f1_u(x,y,z), and f1_p(x,y,z) */
// void (*g0Funcs[NUM_FIELDS*NUM_FIELDS])(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar g0[]); /* g0_uu(x,y,z), g0_up(x,y,z), g0_pu(x,y,z), and g0_pp(x,y,z) */
// void (*g1Funcs[NUM_FIELDS*NUM_FIELDS])(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar g1[]); /* g1_uu(x,y,z), g1_up(x,y,z), g1_pu(x,y,z), and g1_pp(x,y,z) */
// void (*g2Funcs[NUM_FIELDS*NUM_FIELDS])(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar g2[]); /* g2_uu(x,y,z), g2_up(x,y,z), g2_pu(x,y,z), and g2_pp(x,y,z) */
// void (*g3Funcs[NUM_FIELDS*NUM_FIELDS])(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar g3[]); /* g3_uu(x,y,z), g3_up(x,y,z), g3_pu(x,y,z), and g3_pp(x,y,z) */
// void (*exactFuncs[NUM_BASIS_COMPONENTS_TOTAL])(const PetscReal x[], PetscScalar *u); /* The exact solution function u(x,y,z), v(x,y,z), p(x,y,z), and T(x,y,z) */
// BCType      bcType;              /* The type of boundary conditions */
// ForcingType forcingType;         /* The type of rhs */
// } AppCtx;
//
// void zero(const PetscReal coords[], PetscScalar *u)
// {
// *u = 0.0;
// }
//
// /*
// In 2D, for constant forcing,
//
// f_x = f_y = 3
//
// we use the exact solution,
//
// u = x^2 + y^2
// v = 2 x^2 - 2xy
// p = x + y - 1
// T = x + y
//
// so that
//
// -\Delta u + \nabla p + f = <-4, -4> + <1, 1> + <3, 3> = 0
// \nabla \cdot u           = 2x - 2x                    = 0
// #NAME?
// */
// void quadratic_u_2d(const PetscReal x[], PetscScalar *u)
// {
// *u = x[0]*x[0] + x[1]*x[1];
// };
//
// void quadratic_v_2d(const PetscReal x[], PetscScalar *v)
// {
// *v = 2.0*x[0]*x[0] - 2.0*x[0]*x[1];
// };
//
// void linear_p_2d(const PetscReal x[], PetscScalar *p)
// {
// *p = x[0] + x[1] - 1.0;
// };
//
// void linear_T_2d(const PetscReal x[], PetscScalar *T)
// {
// *T = x[0] + x[1];
// };
//
// /*
// In 2D, for linear forcing,
//
// f_x =  3 - 8y
// f_y = -5 + 8x
//
// we use the exact solution,
//
// u =  2 x (x-1) (1 - 2 y)
// v = -2 y (y-1) (1 - 2 x)
// p = x + y - 1
// T = x + y
//
// so that
//
// -\Delta u + \nabla p + f = <-4+8y, 4-8x> + <1, 1> + <3-8y, 8x-5> = 0
// \nabla \cdot u           = (4x-2) (1-2y) - (4y-2) (1-2x)         = 0
// #NAME?
// */
// void cubic_u_2d(const PetscReal x[], PetscScalar *u)
// {
// *u = 2.0*x[0]*(x[0]-1.0)*(1.0 - 2.0*x[1]);
// };
//
// void cubic_v_2d(const PetscReal x[], PetscScalar *v)
// {
// *v = -2.0*x[1]*(x[1]-1.0)*(1.0 - 2.0*x[0]);
// };
//
// /*
// Let \sigma = (\nabla u + \nabla u^T) = < \sigma_{ij} >, where \sigma_{ij} = \sigma_{ji}
// Then at the top and bottom (t = <1,0>),
// <\sigma_{00}, \sigma_{01}> = 0 so \sigma_{00} = A(x,y) y(1-y) \sigma_{01} = B(x,y) y(1-y)
// Using the left and right (t = <0,1>),
// <\sigma_{10}, \sigma_{11}> = 0 so \sigma_{10} = C(x,y) x(1-x) \sigma_{11} = D(x,y) x(1-x)
// Which means
// \sigma_{00} = A(x,y) y(1-y)        = 2 u_x
// \sigma_{01} = E(x,y) x(1-x) y(1-y) = u_y + v_x
// \sigma_{11} = D(x,y) x(1-x)        = 2 v_y
// Also we have
// u(x=0,1) = 0 ==> u = A'(x,y) x(1-x)
// v(y=0,1) = 0 ==> v = D'(x,y) y(1-y)
// Thus we need
// \int x - x^2 = x^2/2 - x^3/3 + C ==> 3 x^2 - 2 x^3 + 1 = 0 at x=0,1
// so that
// u =  (3 x^2 - 2 x^3 + 1) y(1-y)
// v = -(3 y^2 - 2 y^3 + 1) x(1-x)
// u_x =  6 x(1-x) y(1-y)
// v_y = -6 x(1-x) y(1-y)
// u_xx =  6 (1-2x) y(1-y)
// v_yy = -6 (1-2y) x(1-x)
//
// In 2D, for cubic forcing,
//
// f_x = -1 + 6 (1-2x) y(1-y)
// f_y = -1 - 6 (1-2y) x(1-x)
//
// we use the exact solution,
//
// u =  (3 x^2 - 2 x^3 + 1) y(1-y)
// v = -(3 y^2 - 2 y^3 + 1) x(1-x)
// p = x + y - 1
// T = x + y
//
// so that
//
// -\Delta u + \nabla p + f = <-6 (1-2x) y(1-y), 6 (1-2y) x(1-x)> + <1, 1> + <-1 + 6 (1-2x) y(1-y), -1 - 6 (1-2y) x(1-x)> = 0
// \nabla \cdot u           = 6 x(1-x) y(1-y) -6 (1-2y) x(1-x) = 0
// #NAME?
// */
// void quintic_u_2d(const PetscReal x[], PetscScalar *u)
// {
// *u = (3.0*x[0]*x[0] - 2.0*x[0]*x[0]*x[0] + 1.0)*x[1]*(1.0-x[1]);
// };
//
// void quintic_v_2d(const PetscReal x[], PetscScalar *v)
// {
// *v = -(3.0*x[1]*x[1] - 2.0*x[1]*x[1]*x[1] + 1.0)*x[0]*(1.0-x[0]);
// };
//
// void f0_u_constant(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f0[])
// {
// const PetscInt Ncomp = NUM_BASIS_COMPONENTS_0;
// PetscInt       comp;
//
// for (comp = 0; comp < Ncomp; ++comp) f0[comp] = 3.0;
// }
//
// void f0_u_linear_2d(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f0[])
// {
// f0[0] =  3.0 - 8.0*x[1];
// f0[1] = -5.0 + 8.0*x[0];
// }
//
// void f0_u_cubic_2d(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f0[])
// {
// f0[0] = -1.0 + 6.0*(1.0 - 2.0*x[0])*x[1]*(1.0 - x[1]);
// f0[1] = -1.0 - 6.0*(1.0 - 2.0*x[1])*x[0]*(1.0 - x[0]);
// }
//
// /* gradU[comp*dim+d] = {u_x, u_y, v_x, v_y} or {u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z}
// u[Ncomp]          = {p} */
// void f1_u(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f1[])
// {
// const PetscInt dim   = SPATIAL_DIM_0;
// const PetscInt Ncomp = NUM_BASIS_COMPONENTS_0;
// PetscInt       comp, d;
//
// for (comp = 0; comp < Ncomp; ++comp) {
// for (d = 0; d < dim; ++d) {
// /* f1[comp*dim+d] = 0.5*(gradU[comp*dim+d] + gradU[d*dim+comp]); */
// f1[comp*dim+d] = gradU[comp*dim+d];
// }
// f1[comp*dim+comp] -= u[Ncomp];
// }
// }
//
// /* gradU[comp*dim+d] = {u_x, u_y, v_x, v_y} or {u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z} */
// void f0_p(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f0[])
// {
// const PetscInt dim = SPATIAL_DIM_0;
// PetscInt       d;
//
// f0[0] = 0.0;
// for (d = 0; d < dim; ++d) f0[0] += gradU[d*dim+d];
// }
//
// void f1_p(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f1[])
// {
// const PetscInt dim = SPATIAL_DIM_0;
// PetscInt       d;
//
// for (d = 0; d < dim; ++d) f1[d] = 0.0;
// }
//
// void f0_T(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f0[])
// {
// f0[0] = 0.0;
// }
//
// void f1_T(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar f1[])
// {
// const PetscInt dim = SPATIAL_DIM_2;
// const PetscInt off = SPATIAL_DIM_0*NUM_BASIS_COMPONENTS_0+SPATIAL_DIM_1*NUM_BASIS_COMPONENTS_1;
// PetscInt       d;
//
// for (d = 0; d < dim; ++d) f1[d] = gradU[off+d];
// }
//
// /* < v_t, I t > */
// void g0_TT(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar g0[])
// {
// g0[0] = 1.0;
// }
//
// /* < q, \nabla\cdot v >
// NcompI = 1, NcompJ = dim */
// void g1_pu(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar g1[])
// {
// const PetscInt dim = SPATIAL_DIM_0;
// PetscInt       d;
//
// for (d = 0; d < dim; ++d) g1[d*dim+d] = 1.0; /* \frac{\partial\phi^{u_d}}{\partial x_d} */
// }
//
// /* -< \nabla\cdot v, p >
// NcompI = dim, NcompJ = 1 */
// void g2_up(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar g2[])
// {
// const PetscInt dim = SPATIAL_DIM_0;
// PetscInt       d;
//
// for (d = 0; d < dim; ++d) g2[d*dim+d] = -1.0; /* \frac{\partial\psi^{u_d}}{\partial x_d} */
// }
//
// /* < \nabla v, \nabla u + {\nabla u}^T >
// This just gives \nabla u, give the perdiagonal for the transpose */
// void g3_uu(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar g3[])
// {
// const PetscInt dim   = SPATIAL_DIM_0;
// const PetscInt Ncomp = NUM_BASIS_COMPONENTS_0;
// PetscInt       compI, d;
//
// for (compI = 0; compI < Ncomp; ++compI) {
// for (d = 0; d < dim; ++d) {
// g3[((compI*Ncomp+compI)*dim+d)*dim+d] = 1.0;
// }
// }
// }
//
// /* < \nabla t, \nabla T + {\nabla u}^T >
// This just gives \nabla T, give the perdiagonal for the transpose */
// void g3_TT(const PetscScalar u[], const PetscScalar gradU[], const PetscReal x[], PetscScalar g3[])
// {
// const PetscInt dim   = SPATIAL_DIM_2;
// const PetscInt Ncomp = NUM_BASIS_COMPONENTS_2;
// PetscInt       compI, d;
//
// for (compI = 0; compI < Ncomp; ++compI) {
// for (d = 0; d < dim; ++d) {
// g3[((compI*Ncomp+compI)*dim+d)*dim+d] = 1.0;
// }
// }
// }
//
// /*
// In 3D we use exact solution:
//
// u = x^2 + y^2
// v = y^2 + z^2
// w = x^2 + y^2 - 2(x+y)z
// p = x + y + z - 3/2
// f_x = f_y = f_z = 3
//
// so that
//
// -\Delta u + \nabla p + f = <-4, -4, -4> + <1, 1, 1> + <3, 3, 3> = 0
// \nabla \cdot u           = 2x + 2y - 2(x + y)                   = 0
// */
// void quadratic_u_3d(const PetscReal x[], PetscScalar *u)
// {
// *u = x[0]*x[0] + x[1]*x[1];
// };
//
// void quadratic_v_3d(const PetscReal x[], PetscScalar *v)
// {
// *v = x[1]*x[1] + x[2]*x[2];
// };
//
// void quadratic_w_3d(const PetscReal x[], PetscScalar *w)
// {
// *w = x[0]*x[0] + x[1]*x[1] - 2.0*(x[0] + x[1])*x[2];
// };
//
// void linear_p_3d(const PetscReal x[], PetscScalar *p)
// {
// *p = x[0] + x[1] + x[2] - 1.5;
// };
//
// PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
// {
// const char     *bcTypes[2]      = {"dirichlet", "freeslip"};
// const char     *forcingTypes[3] = {"constant", "linear", "cubic"};
// const char     *runTypes[2]     = {"full", "test"};
// PetscInt       bc, forcing, run;
//
// options->debug           = 0;
// options->runType         = RUN_FULL;
// options->dim             = 2;
// options->interpolate     = PETSC_FALSE;
// options->refinementLimit = 0.0;
// options->bcType          = DIRICHLET;
// options->forcingType     = FORCING_CONSTANT;
// options->numBatches      = 1;
// options->numBlocks       = 1;
// options->jacobianMF      = PETSC_FALSE;
// options->showInitial     = PETSC_FALSE;
// options->showSolution    = PETSC_TRUE;
//
// options->fem.quad    = (PetscQuadrature*) &options->q;
// options->fem.f0Funcs = (void (**)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[])) &options->f0Funcs;
// options->fem.f1Funcs = (void (**)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[])) &options->f1Funcs;
// options->fem.g0Funcs = (void (**)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[])) &options->g0Funcs;
// options->fem.g1Funcs = (void (**)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[])) &options->g1Funcs;
// options->fem.g2Funcs = (void (**)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[])) &options->g2Funcs;
// options->fem.g3Funcs = (void (**)(const PetscScalar[], const PetscScalar[], const PetscReal[], PetscScalar[])) &options->g3Funcs;
//
// MPI_Comm_size(comm, &options->numProcs);
// MPI_Comm_rank(comm, &options->rank);
// PetscOptionsBegin(comm, "", "Stokes Problem Options", "DMPLEX");
// PetscOptionsInt("-debug", "The debugging level", "ex31.c", options->debug, &options->debug, NULL);
// run  = options->runType;
// PetscOptionsEList("-run_type", "The run type", "ex31.c", runTypes, 2, runTypes[options->runType], &run, NULL);
//
// options->runType = (RunType) run;
//
// PetscOptionsInt("-dim", "The topological mesh dimension", "ex31.c", options->dim, &options->dim, NULL);
// PetscOptionsBool("-interpolate", "Generate intermediate mesh elements", "ex31.c", options->interpolate, &options->interpolate, NULL);
// PetscOptionsReal("-refinement_limit", "The largest allowable cell volume", "ex31.c", options->refinementLimit, &options->refinementLimit, NULL);
// PetscStrcpy(options->partitioner, "chaco");
// PetscOptionsString("-partitioner", "The graph partitioner", "pflotran.cxx", options->partitioner, options->partitioner, 2048, NULL);
// bc   = options->bcType;
// PetscOptionsEList("-bc_type","Type of boundary condition","ex31.c",bcTypes,2,bcTypes[options->bcType],&bc,NULL);
//
// options->bcType = (BCType) bc;
// forcing         = options->forcingType;
//
// PetscOptionsEList("-forcing_type","Type of forcing function","ex31.c",forcingTypes,3,forcingTypes[options->forcingType],&forcing,NULL);
//
// options->forcingType = (ForcingType) forcing;
//
// PetscOptionsInt("-gpu_batches", "The number of cell batches per kernel", "ex31.c", options->numBatches, &options->numBatches, NULL);
// PetscOptionsInt("-gpu_blocks", "The number of concurrent blocks per kernel", "ex31.c", options->numBlocks, &options->numBlocks, NULL);
// PetscOptionsBool("-jacobian_mf", "Calculate the action of the Jacobian on the fly", "ex31.c", options->jacobianMF, &options->jacobianMF, NULL);
// PetscOptionsBool("-show_initial", "Output the initial guess for verification", "ex31.c", options->showInitial, &options->showInitial, NULL);
// PetscOptionsBool("-show_solution", "Output the solution for verification", "ex31.c", options->showSolution, &options->showSolution, NULL);
// PetscOptionsEnd();
//
// PetscLogEventRegister("CreateMesh", DM_CLASSID, &options->createMeshEvent);
// return(0);
// };
//
// PetscErrorCode DMVecViewLocal(DM dm, Vec v, PetscViewer viewer)
// {
// Vec            lv;
// PetscInt       p;
// PetscMPIInt    rank, numProcs;
//
// MPI_Comm_rank(PetscObjectComm((PetscObject)dm), &rank);
// MPI_Comm_size(PetscObjectComm((PetscObject)dm), &numProcs);
// DMGetLocalVector(dm, &lv);
// DMGlobalToLocalBegin(dm, v, INSERT_VALUES, lv);
// DMGlobalToLocalEnd(dm, v, INSERT_VALUES, lv);
// PetscPrintf(PETSC_COMM_WORLD, "Local function\n");
// for (p = 0; p < numProcs; ++p) {
// if (p == rank) {VecView(lv, PETSC_VIEWER_STDOUT_SELF);}
// PetscBarrier((PetscObject) dm);
// }
// DMRestoreLocalVector(dm, &lv);
// return(0);
// }
//
// PetscErrorCode CreateMesh(MPI_Comm comm, AppCtx *user, DM *dm)
// {
// PetscInt       dim             = user->dim;
// PetscBool      interpolate     = user->interpolate;
// PetscReal      refinementLimit = user->refinementLimit;
// const char     *partitioner    = user->partitioner;
//
// PetscLogEventBegin(user->createMeshEvent,0,0,0,0);
// DMPlexCreateBoxMesh(comm, dim, interpolate, dm);
// {
// DM refinedMesh     = NULL;
// DM distributedMesh = NULL;
//
// /* Refine mesh using a volume constraint */
// DMPlexSetRefinementLimit(*dm, refinementLimit);
// DMRefine(*dm, comm, &refinedMesh);
// if (refinedMesh) {
// DMDestroy(dm);
// *dm  = refinedMesh;
// }
// /* Distribute mesh over processes */
// DMPlexDistribute(*dm, partitioner, 0, &distributedMesh);
// if (distributedMesh) {
// DMDestroy(dm);
// *dm  = distributedMesh;
// }
// }
// DMSetFromOptions(*dm);
// PetscLogEventEnd(user->createMeshEvent,0,0,0,0);
// user->dm = *dm;
// return(0);
// }
//
// PetscErrorCode PointOnBoundary_2D(const PetscScalar coords[], PetscBool onBd[])
// {
// const PetscInt  corner = 0, bottom = 1, right = 2, top = 3, left = 4;
// const PetscReal eps    = 1.0e-10;
//
// onBd[bottom] = PetscAbsScalar(coords[1])       < eps ? PETSC_TRUE : PETSC_FALSE;
// onBd[right]  = PetscAbsScalar(coords[0] - 1.0) < eps ? PETSC_TRUE : PETSC_FALSE;
// onBd[top]    = PetscAbsScalar(coords[1] - 1.0) < eps ? PETSC_TRUE : PETSC_FALSE;
// onBd[left]   = PetscAbsScalar(coords[0])       < eps ? PETSC_TRUE : PETSC_FALSE;
// onBd[corner] = onBd[bottom] + onBd[right] + onBd[top] + onBd[left] > 1 ? PETSC_TRUE : PETSC_FALSE;
// return(0);
// }
//
// PetscErrorCode CreateBoundaryPointIS_Square(DM dm, PetscInt *numBoundaries, PetscInt **numBoundaryConstraints, IS **boundaryPoints, IS **constraintIndices)
// {
// MPI_Comm       comm;
// PetscSection   coordSection;
// Vec            coordinates;
// PetscScalar    *coords;
// PetscInt       vStart, vEnd;
// IS             bcPoints;
// const PetscInt *points;
// const PetscInt corner               = 0, bottom = 1, right = 2, top = 3, left = 4;
// PetscInt       numBoundaryPoints[5] = {0, 0, 0, 0, 0}, bd, numPoints, p;
// PetscInt       *bdPoints[5], *idx;
//
// PetscObjectGetComm((PetscObject)dm,&comm);
// DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);
// /* boundary 0: corners
// boundary 1: bottom
// boundary 2: right
// boundary 3: top
// boundary 4: left
// */
// *numBoundaries = 5;
//
// PetscMalloc(*numBoundaries * sizeof(PetscInt), numBoundaryConstraints);
// PetscMalloc(*numBoundaries * sizeof(IS), boundaryPoints);
// PetscMalloc(*numBoundaries * sizeof(IS), constraintIndices);
//
// /* Set number of constraints for each boundary */
// (*numBoundaryConstraints)[corner] = 2;
// (*numBoundaryConstraints)[bottom] = 1;
// (*numBoundaryConstraints)[right]  = 1;
// (*numBoundaryConstraints)[top]    = 1;
// (*numBoundaryConstraints)[left]   = 1;
// /* Set local constraint indices for each boundary */
// PetscMalloc((*numBoundaryConstraints)[corner] * sizeof(PetscInt), &idx);
// idx[0] = 0; idx[1] = 1;
// ISCreateGeneral(comm, (*numBoundaryConstraints)[corner], idx, PETSC_OWN_POINTER, &(*constraintIndices)[corner]);
// PetscMalloc((*numBoundaryConstraints)[bottom] * sizeof(PetscInt), &idx);
// idx[0] = 1;
// ISCreateGeneral(comm, (*numBoundaryConstraints)[bottom], idx, PETSC_OWN_POINTER, &(*constraintIndices)[bottom]);
// PetscMalloc((*numBoundaryConstraints)[right] * sizeof(PetscInt), &idx);
// idx[0] = 0;
// ISCreateGeneral(comm, (*numBoundaryConstraints)[right], idx, PETSC_OWN_POINTER, &(*constraintIndices)[right]);
// PetscMalloc((*numBoundaryConstraints)[top] * sizeof(PetscInt), &idx);
// idx[0] = 1;
// ISCreateGeneral(comm, (*numBoundaryConstraints)[top], idx, PETSC_OWN_POINTER, &(*constraintIndices)[top]);
// PetscMalloc((*numBoundaryConstraints)[left] * sizeof(PetscInt), &idx);
// idx[0] = 0;
// ISCreateGeneral(comm, (*numBoundaryConstraints)[left], idx, PETSC_OWN_POINTER, &(*constraintIndices)[left]);
//
// /* Count points on each boundary */
// DMPlexGetCoordinateSection(dm, &coordSection);
// DMGetCoordinatesLocal(dm, &coordinates);
// VecGetArray(coordinates, &coords);
// DMPlexGetStratumIS(dm, "marker", 1, &bcPoints);
// ISGetLocalSize(bcPoints, &numPoints);
// ISGetIndices(bcPoints, &points);
// for (p = 0; p < numPoints; ++p) {
// PetscBool onBd[5];
// PetscInt  off, bd;
//
// if ((points[p] >= vStart) && (points[p] < vEnd)) {
// PetscSectionGetOffset(coordSection, points[p], &off);
// PointOnBoundary_2D(&coords[off], onBd);
// } else {
// PetscInt *closure = NULL;
// PetscInt closureSize, q, r;
//
// DMPlexGetTransitiveClosure(dm, points[p], PETSC_TRUE, &closureSize, &closure);
// /* Compress out non-vertices */
// for (q = 0, r = 0; q < closureSize*2; q += 2) {
// if ((closure[q] >= vStart) && (closure[q] < vEnd)) {
// closure[r] = closure[q];
// ++r;
// }
// }
// closureSize = r;
// for (q = 0; q < closureSize; ++q) {
// PetscSectionGetOffset(coordSection, closure[q], &off);
// PointOnBoundary_2D(&coords[off], onBd);
// if (!onBd[corner]) break;
// }
// DMPlexRestoreTransitiveClosure(dm, points[p], PETSC_TRUE, &closureSize, &closure);
// if (q == closureSize) SETERRQ1(comm, PETSC_ERR_PLIB, "Cannot handle face %d which has every vertex on a corner", points[p]);
// }
//
// for (bd = 0; bd < 5; ++bd) {
// if (onBd[bd]) {
// ++numBoundaryPoints[bd];
// break;
// }
// }
// }
// /* Set points on each boundary */
// for (bd = 0; bd < 5; ++bd) {
// PetscMalloc(numBoundaryPoints[bd] * sizeof(PetscInt), &bdPoints[bd]);
// numBoundaryPoints[bd] = 0;
// }
// for (p = 0; p < numPoints; ++p) {
// PetscBool onBd[5];
// PetscInt  off, bd;
//
// if ((points[p] >= vStart) && (points[p] < vEnd)) {
// PetscSectionGetOffset(coordSection, points[p], &off);
// PointOnBoundary_2D(&coords[off], onBd);
// } else {
// PetscInt *closure = NULL;
// PetscInt closureSize, q, r;
//
// DMPlexGetTransitiveClosure(dm, points[p], PETSC_TRUE, &closureSize, &closure);
// /* Compress out non-vertices */
// for (q = 0, r = 0; q < closureSize*2; q += 2) {
// if ((closure[q] >= vStart) && (closure[q] < vEnd)) {
// closure[r] = closure[q];
// ++r;
// }
// }
// closureSize = r;
// for (q = 0; q < closureSize; ++q) {
// PetscSectionGetOffset(coordSection, closure[q], &off);
// PointOnBoundary_2D(&coords[off], onBd);
// if (!onBd[corner]) break;
// }
// DMPlexRestoreTransitiveClosure(dm, points[p], PETSC_TRUE, &closureSize, &closure);
// if (q == closureSize) SETERRQ1(comm, PETSC_ERR_PLIB, "Cannot handle face %d which has every vertex on a corner", points[p]);
// }
//
// for (bd = 0; bd < 5; ++bd) {
// if (onBd[bd]) {
// bdPoints[bd][numBoundaryPoints[bd]++] = points[p];
// break;
// }
// }
// }
// VecRestoreArray(coordinates, &coords);
// ISRestoreIndices(bcPoints, &points);
// ISDestroy(&bcPoints);
// for (bd = 0; bd < 5; ++bd) {
// ISCreateGeneral(comm, numBoundaryPoints[bd], bdPoints[bd], PETSC_OWN_POINTER, &(*boundaryPoints)[bd]);
// }
// return(0);
// }
//
// PetscErrorCode CreateBoundaryPointIS_Cube(DM dm, PetscInt *numBoundaries, PetscInt **numBoundaryConstraints, IS **boundaryPoints, IS **constraintIndices)
// {
// SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Just lazy");
// return(0);
// }
//
// /* This will only work for the square/cube, but I think the interface is robust */
// PetscErrorCode CreateBoundaryPointIS(DM dm, PetscInt *numBoundaries, PetscInt **numBoundaryConstraints, IS **boundaryPoints, IS **constraintIndices)
// {
// PetscInt       dim;
//
// DMPlexGetDimension(dm, &dim);
// switch (dim) {
// case 2:
// CreateBoundaryPointIS_Square(dm, numBoundaries, numBoundaryConstraints, boundaryPoints, constraintIndices);
// break;
// case 3:
// CreateBoundaryPointIS_Cube(dm, numBoundaries, numBoundaryConstraints, boundaryPoints, constraintIndices);
// break;
// default:
// SETERRQ1(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "No boundary creatin routine for dimension %d", dim);
// }
// return(0);
// }
//
// PetscErrorCode SetupQuadrature(AppCtx *user)
// {
// user->fem.quad[0].numQuadPoints = NUM_QUADRATURE_POINTS_0;
// user->fem.quad[0].quadPoints    = points_0;
// user->fem.quad[0].quadWeights   = weights_0;
// user->fem.quad[0].numBasisFuncs = NUM_BASIS_FUNCTIONS_0;
// user->fem.quad[0].numComponents = NUM_BASIS_COMPONENTS_0;
// user->fem.quad[0].basis         = Basis_0;
// user->fem.quad[0].basisDer      = BasisDerivatives_0;
// user->fem.quad[1].numQuadPoints = NUM_QUADRATURE_POINTS_1;
// user->fem.quad[1].quadPoints    = points_1;
// user->fem.quad[1].quadWeights   = weights_1;
// user->fem.quad[1].numBasisFuncs = NUM_BASIS_FUNCTIONS_1;
// user->fem.quad[1].numComponents = NUM_BASIS_COMPONENTS_1;
// user->fem.quad[1].basis         = Basis_1;
// user->fem.quad[1].basisDer      = BasisDerivatives_1;
// user->fem.quad[2].numQuadPoints = NUM_QUADRATURE_POINTS_2;
// user->fem.quad[2].quadPoints    = points_2;
// user->fem.quad[2].quadWeights   = weights_2;
// user->fem.quad[2].numBasisFuncs = NUM_BASIS_FUNCTIONS_2;
// user->fem.quad[2].numComponents = NUM_BASIS_COMPONENTS_2;
// user->fem.quad[2].basis         = Basis_2;
// user->fem.quad[2].basisDer      = BasisDerivatives_2;
// return(0);
// }
//
// /*
// There is a problem here with uninterpolated meshes. The index in numDof[] is not dimension in this case,
// but sieve depth.
// */
// PetscErrorCode SetupSection(DM dm, AppCtx *user)
// {
// PetscSection   section;
// const PetscInt numFields           = NUM_FIELDS;
// PetscInt       dim                 = user->dim;
// PetscInt       numBC               = 0;
// PetscInt       numComp[NUM_FIELDS] = {NUM_BASIS_COMPONENTS_0, NUM_BASIS_COMPONENTS_1, NUM_BASIS_COMPONENTS_2};
// PetscInt       bcFields[2]         = {0, 2};
// IS             bcPoints[2]         = {NULL, NULL};
// PetscInt       numDof[NUM_FIELDS*(SPATIAL_DIM_0+1)];
// PetscInt       f, d;
// PetscBool      view;
//
// if (dim != SPATIAL_DIM_0) SETERRQ2(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_SIZ, "Spatial dimension %d should be %d", dim, SPATIAL_DIM_0);
// if (dim != SPATIAL_DIM_1) SETERRQ2(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_SIZ, "Spatial dimension %d should be %d", dim, SPATIAL_DIM_1);
// for (d = 0; d <= dim; ++d) {
// numDof[0*(dim+1)+d] = numDof_0[d];
// numDof[1*(dim+1)+d] = numDof_1[d];
// numDof[2*(dim+1)+d] = numDof_2[d];
// }
// for (f = 0; f < numFields; ++f) {
// for (d = 1; d < dim; ++d) {
// if ((numDof[f*(dim+1)+d] > 0) && !user->interpolate) SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "Mesh must be interpolated when unknowns are specified on edges or faces.");
// }
// }
// if (user->bcType == FREE_SLIP) {
// PetscInt numBoundaries, b;
// PetscInt *numBoundaryConstraints;
// IS       *boundaryPoints, *constraintIndices;
//
// DMPlexCreateSectionInitial(dm, dim, numFields, numComp, numDof, &section);
// /* Velocity conditions */
// CreateBoundaryPointIS(dm, &numBoundaries, &numBoundaryConstraints, &boundaryPoints, &constraintIndices);
// for (b = 0; b < numBoundaries; ++b) {
// DMPlexCreateSectionBCDof(dm, 1, &bcFields[0], &boundaryPoints[b], numBoundaryConstraints[b], section);
// }
// /* Temperature conditions */
// DMPlexGetStratumIS(dm, "marker", 1, &bcPoints[0]);
// DMPlexCreateSectionBCDof(dm, 1, &bcFields[1], &bcPoints[0], PETSC_DETERMINE, section);
// PetscSectionSetUp(section);
// for (b = 0; b < numBoundaries; ++b) {
// DMPlexCreateSectionBCIndicesField(dm, bcFields[0], boundaryPoints[b], constraintIndices[b], section);
// }
// DMPlexCreateSectionBCIndicesField(dm, bcFields[1], bcPoints[0], NULL, section);
// DMPlexCreateSectionBCIndices(dm, section);
// } else {
// if (user->bcType == DIRICHLET) {
// numBC       = 2;
// DMPlexGetStratumIS(dm, "marker", 1, &bcPoints[0]);
// bcPoints[1] = bcPoints[0];
// PetscObjectReference((PetscObject) bcPoints[1]);
// }
// DMPlexCreateSection(dm, dim, numFields, numComp, numDof, numBC, bcFields, bcPoints, &section);
// }
// PetscSectionSetFieldName(section, 0, "velocity");
// PetscSectionSetFieldName(section, 1, "pressure");
// PetscSectionSetFieldName(section, 2, "temperature");
// DMSetDefaultSection(dm, section);
// ISDestroy(&bcPoints[0]);
// ISDestroy(&bcPoints[1]);
// PetscOptionsHasName(((PetscObject) dm)->prefix, "-section_view", &view);
// if ((user->bcType == FREE_SLIP) && view) {
// PetscSection s, gs;
//
// DMGetDefaultSection(dm, &s);
// PetscSectionView(s, PETSC_VIEWER_STDOUT_WORLD);
// DMGetDefaultGlobalSection(dm, &gs);
// PetscSectionView(gs, PETSC_VIEWER_STDOUT_WORLD);
// }
// return(0);
// }
//
// PetscErrorCode SetupExactSolution(DM dm, AppCtx *user)
// {
// PetscFEM       *fem = &user->fem;
//
// switch (user->forcingType) {
// case FORCING_CONSTANT:
// if (user->bcType == FREE_SLIP) SETERRQ(PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONG, "Constant forcing is incompatible with freeslip boundary conditions");
// fem->f0Funcs[0] = f0_u_constant;
// break;
// case FORCING_LINEAR:
// switch (user->bcType) {
// case DIRICHLET:
// case FREE_SLIP:
// switch (user->dim) {
// case 2:
// fem->f0Funcs[0] = f0_u_linear_2d;
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid dimension %d", user->dim);
// }
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid boundary condition type %d", user->bcType);
// }
// break;
// case FORCING_CUBIC:
// switch (user->bcType) {
// case DIRICHLET:
// case FREE_SLIP:
// switch (user->dim) {
// case 2:
// fem->f0Funcs[0] = f0_u_cubic_2d;
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid dimension %d", user->dim);
// }
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid boundary condition type %d", user->bcType);
// }
// break;
// }
// fem->f0Funcs[1] = f0_p;
// fem->f0Funcs[2] = f0_T;
// fem->f1Funcs[0] = f1_u;
// fem->f1Funcs[1] = f1_p;
// fem->f1Funcs[2] = f1_T;
// fem->g0Funcs[0] = NULL;
// fem->g0Funcs[1] = NULL;
// fem->g0Funcs[2] = NULL;
// fem->g0Funcs[3] = NULL;
// fem->g0Funcs[4] = NULL;
// fem->g0Funcs[5] = NULL;
// fem->g0Funcs[6] = NULL;
// fem->g0Funcs[7] = NULL;
// fem->g0Funcs[8] = NULL;
// fem->g1Funcs[0] = NULL;
// fem->g1Funcs[1] = NULL;
// fem->g1Funcs[2] = NULL;
// fem->g1Funcs[3] = g1_pu;      /* < q, \nabla\cdot v > */
// fem->g1Funcs[4] = NULL;
// fem->g1Funcs[5] = NULL;
// fem->g1Funcs[6] = NULL;
// fem->g1Funcs[7] = NULL;
// fem->g1Funcs[8] = NULL;
// fem->g2Funcs[0] = NULL;
// fem->g2Funcs[1] = g2_up;      /* < \nabla\cdot v, p > */
// fem->g2Funcs[2] = NULL;
// fem->g2Funcs[3] = NULL;
// fem->g2Funcs[4] = NULL;
// fem->g2Funcs[5] = NULL;
// fem->g2Funcs[6] = NULL;
// fem->g2Funcs[7] = NULL;
// fem->g2Funcs[8] = NULL;
// fem->g3Funcs[0] = g3_uu;      /* < \nabla v, \nabla u + {\nabla u}^T > */
// fem->g3Funcs[1] = NULL;
// fem->g3Funcs[2] = NULL;
// fem->g3Funcs[3] = NULL;
// fem->g3Funcs[4] = NULL;
// fem->g3Funcs[5] = NULL;
// fem->g3Funcs[6] = NULL;
// fem->g3Funcs[7] = NULL;
// fem->g3Funcs[8] = g3_TT;      /* < \nabla t, \nabla T + {\nabla T}^T > */
// switch (user->forcingType) {
// case FORCING_CONSTANT:
// switch (user->bcType) {
// case DIRICHLET:
// switch (user->dim) {
// case 2:
// user->exactFuncs[0] = quadratic_u_2d;
// user->exactFuncs[1] = quadratic_v_2d;
// user->exactFuncs[2] = linear_p_2d;
// user->exactFuncs[3] = linear_T_2d;
// break;
// case 3:
// user->exactFuncs[0] = quadratic_u_3d;
// user->exactFuncs[1] = quadratic_v_3d;
// user->exactFuncs[2] = quadratic_w_3d;
// user->exactFuncs[3] = linear_p_3d;
// user->exactFuncs[4] = linear_T_2d;
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid dimension %d", user->dim);
// }
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid boundary condition type %d", user->bcType);
// }
// break;
// case FORCING_LINEAR:
// switch (user->bcType) {
// case DIRICHLET:
// switch (user->dim) {
// case 2:
// user->exactFuncs[0] = cubic_u_2d;
// user->exactFuncs[1] = cubic_v_2d;
// user->exactFuncs[2] = linear_p_2d;
// user->exactFuncs[3] = linear_T_2d;
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid dimension %d", user->dim);
// }
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid boundary condition type %d", user->bcType);
// }
// break;
// case FORCING_CUBIC:
// switch (user->bcType) {
// case DIRICHLET:
// case FREE_SLIP:
// switch (user->dim) {
// case 2:
// user->exactFuncs[0] = quintic_u_2d;
// user->exactFuncs[1] = quintic_v_2d;
// user->exactFuncs[2] = linear_p_2d;
// user->exactFuncs[3] = linear_T_2d;
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid dimension %d", user->dim);
// }
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid boundary condition type %d", user->bcType);
// }
// break;
// default:
// SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Invalid forcing type %d", user->forcingType);
// }
// DMPlexSetFEMIntegration(dm, FEMIntegrateResidualBatch, FEMIntegrateJacobianActionBatch, FEMIntegrateJacobianBatch);
// return(0);
// }
//
// /*
// . field - The field whose diagonal block (of the Jacobian) has this null space
// */
// PetscErrorCode CreateNullSpaces(DM dm, PetscInt field, MatNullSpace *nullSpace)
// {
// AppCtx         *user;
// Vec            nullVec, localNullVec;
// PetscSection   section;
// PetscScalar    *a;
// PetscInt       pressure = field;
// PetscInt       pStart, pEnd, p;
//
// DMGetApplicationContext(dm, (void**) &user);
// DMGetGlobalVector(dm, &nullVec);
// DMGetLocalVector(dm, &localNullVec);
// VecSet(nullVec, 0.0);
// /* Put a constant in for all pressures */
// DMGetDefaultSection(dm, &section);
// PetscSectionGetChart(section, &pStart, &pEnd);
// VecGetArray(localNullVec, &a);
// for (p = pStart; p < pEnd; ++p) {
// PetscInt fDim, off, d;
//
// PetscSectionGetFieldDof(section, p, pressure, &fDim);
// PetscSectionGetFieldOffset(section, p, pressure, &off);
// for (d = 0; d < fDim; ++d) a[off+d] = 1.0;
// }
// VecRestoreArray(localNullVec, &a);
// DMLocalToGlobalBegin(dm, localNullVec, INSERT_VALUES, nullVec);
// DMLocalToGlobalEnd(dm, localNullVec, INSERT_VALUES, nullVec);
// DMRestoreLocalVector(dm, &localNullVec);
// VecNormalize(nullVec, NULL);
// if (user->debug) {
// PetscPrintf(PetscObjectComm((PetscObject)dm), "Pressure Null Space\n");
// VecView(nullVec, PETSC_VIEWER_STDOUT_WORLD);
// }
// MatNullSpaceCreate(PetscObjectComm((PetscObject)dm), PETSC_FALSE, 1, &nullVec, nullSpace);
// DMRestoreGlobalVector(dm, &nullVec);
// return(0);
// }
//
// /*
// FormJacobianAction - Form the global Jacobian action Y = JX from the global input X
//
// Input Parameters:
// #NAME?
// #NAME?
//
// Output Parameter:
// . Y  - Local output vector
//
// Note:
// We form the residual one batch of elements at a time. This allows us to offload work onto an accelerator,
// like a GPU, or vectorize on a multicore machine.
//
// .seealso: FormJacobianActionLocal()
// */
// PetscErrorCode FormJacobianAction(Mat J, Vec X,  Vec Y)
// {
// JacActionCtx   *ctx;
// DM             dm;
// Vec            dummy, localX, localY;
// PetscInt       N, n;
//
// MatShellGetContext(J, &ctx);
// dm   = ctx->dm;
//
// /* determine whether X = localX */
// DMGetLocalVector(dm, &dummy);
// DMGetLocalVector(dm, &localX);
// DMGetLocalVector(dm, &localY);
// /* TODO: THIS dummy restore is necessary here so that the first available local vector has boundary conditions in it
// I think the right thing to do is have the user put BC into a local vector and give it to us
// */
// DMRestoreLocalVector(dm, &dummy);
// VecGetSize(X, &N);
// VecGetSize(localX, &n);
//
// if (n != N) { /* X != localX */
// VecSet(localX, 0.0);
// DMGlobalToLocalBegin(dm, X, INSERT_VALUES, localX);
// DMGlobalToLocalEnd(dm, X, INSERT_VALUES, localX);
// } else {
// DMRestoreLocalVector(dm, &localX);
// localX = X;
// }
// DMPlexComputeJacobianActionFEM(dm, J, localX, localY, ctx->user);
// if (n != N) {
// DMRestoreLocalVector(dm, &localX);
// }
// VecSet(Y, 0.0);
// DMLocalToGlobalBegin(dm, localY, ADD_VALUES, Y);
// DMLocalToGlobalEnd(dm, localY, ADD_VALUES, Y);
// DMRestoreLocalVector(dm, &localY);
// if (0) {
// Vec       r;
// PetscReal norm;
//
// VecDuplicate(X, &r);
// MatMult(ctx->J, X, r);
// VecAXPY(r, -1.0, Y);
// VecNorm(r, NORM_2, &norm);
// if (norm > 1.0e-8) {
// PetscPrintf(PETSC_COMM_WORLD, "Jacobian Action Input:\n");
// VecView(X, PETSC_VIEWER_STDOUT_WORLD);
// PetscPrintf(PETSC_COMM_WORLD, "Jacobian Action Result:\n");
// VecView(Y, PETSC_VIEWER_STDOUT_WORLD);
// PetscPrintf(PETSC_COMM_WORLD, "Difference:\n");
// VecView(r, PETSC_VIEWER_STDOUT_WORLD);
// SETERRQ1(PetscObjectComm((PetscObject)J), PETSC_ERR_ARG_WRONG, "The difference with assembled multiply is too large %g", norm);
// }
// VecDestroy(&r);
// }
// return(0);
// }
//
// int main(int argc, char **argv)
// {
// MPI_Comm       comm;
// SNES           snes;                 /* nonlinear solver */
// Vec            u,r;                  /* solution, residual vectors */
// Mat            A,J;                  /* Jacobian matrix */
// MatNullSpace   nullSpace = 0;            /* May be necessary for pressure */
// AppCtx         user;                 /* user-defined work context */
// JacActionCtx   userJ;                /* context for Jacobian MF action */
// PetscInt       its;                  /* iterations for convergence */
// PetscReal      error         = 0.0;  /* L_2 error in the solution */
// const PetscInt numComponents = NUM_BASIS_COMPONENTS_TOTAL;
//
// PetscInitialize(&argc, &argv, NULL, help);
// comm = PETSC_COMM_WORLD;
// ProcessOptions(comm, &user);
// SNESCreate(comm, &snes);
// CreateMesh(comm, &user, &user.dm);
// SNESSetDM(snes, user.dm);
// DMSetApplicationContext(user.dm, &user);
//
// SetupExactSolution(user.dm, &user);
// SetupQuadrature(&user);
// SetupSection(user.dm, &user);
//
// DMCreateGlobalVector(user.dm, &u);
// PetscObjectSetName((PetscObject) u, "solution");
// VecDuplicate(u, &r);
//
// DMCreateMatrix(user.dm, MATAIJ, &J);
// if (user.jacobianMF) {
// PetscInt M, m, N, n;
//
// MatGetSize(J, &M, &N);
// MatGetLocalSize(J, &m, &n);
// MatCreate(comm, &A);
// MatSetSizes(A, m, n, M, N);
// MatSetType(A, MATSHELL);
// MatSetUp(A);
// MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) FormJacobianAction);
//
// userJ.dm   = user.dm;
// userJ.J    = J;
// userJ.user = &user;
//
// DMCreateLocalVector(user.dm, &userJ.u);
// MatShellSetContext(A, &userJ);
// } else {
// A = J;
// }
// DMSetNullSpaceConstructor(user.dm, 1, CreateNullSpaces);
//
// DMSNESSetFunctionLocal(user.dm,  (PetscErrorCode (*)(DM,Vec,Vec,void*))DMPlexComputeResidualFEM,&user);
// DMSNESSetJacobianLocal(user.dm,  (PetscErrorCode (*)(DM,Vec,Mat,Mat,MatStructure*,void*))DMPlexComputeJacobianFEM,&user);
//
// SNESSetFromOptions(snes);
//
// {
// KSP               ksp; PC pc; Vec crd_vec;
// const PetscScalar *v;
// PetscInt          i,k,j,mlocal;
// PetscReal         *coords;
//
// SNESGetKSP(snes, &ksp);
// KSPGetPC(ksp, &pc);
// DMGetCoordinatesLocal(user.dm, &crd_vec);
// VecGetLocalSize(crd_vec,&mlocal);
// PetscMalloc(SPATIAL_DIM_0*mlocal*sizeof(*coords),&coords);
// VecGetArrayRead(crd_vec,&v);
// for (k=j=0; j<mlocal; j++) {
// for (i=0; i<SPATIAL_DIM_0; i++,k++) {
// coords[k] = PetscRealPart(v[k]);
// }
// }
// VecRestoreArrayRead(crd_vec,&v);
// PCSetCoordinates(pc, SPATIAL_DIM_0, mlocal, coords);
// PetscFree(coords);
// }
//
// DMPlexProjectFunction(user.dm, numComponents, user.exactFuncs, INSERT_ALL_VALUES, u);
// if (user.showInitial) {DMVecViewLocal(user.dm, u, PETSC_VIEWER_STDOUT_SELF);}
// if (user.runType == RUN_FULL) {
// PetscScalar (*initialGuess[numComponents])(const PetscReal x[]);
// PetscInt c;
//
// for (c = 0; c < numComponents; ++c) initialGuess[c] = zero;
// DMPlexProjectFunction(user.dm, numComponents, initialGuess, INSERT_VALUES, u);
// if (user.showInitial) {DMVecViewLocal(user.dm, u, PETSC_VIEWER_STDOUT_SELF);}
// if (user.debug) {
// PetscPrintf(comm, "Initial guess\n");
// VecView(u, PETSC_VIEWER_STDOUT_WORLD);
// }
// SNESSolve(snes, NULL, u);
// SNESGetIterationNumber(snes, &its);
// PetscPrintf(comm, "Number of SNES iterations = %D\n", its);
// DMPlexComputeL2Diff(user.dm, user.q, user.exactFuncs, u, &error);
// PetscPrintf(comm, "L_2 Error: %.3g\n", error);
// if (user.showSolution) {
// PetscPrintf(comm, "Solution\n");
// VecChop(u, 3.0e-9);
// VecView(u, PETSC_VIEWER_STDOUT_WORLD);
// }
// } else {
// PetscReal res = 0.0;
//
// /* Check discretization error */
// PetscPrintf(comm, "Initial guess\n");
// VecView(u, PETSC_VIEWER_STDOUT_WORLD);
// DMPlexComputeL2Diff(user.dm, user.q, user.exactFuncs, u, &error);
// PetscPrintf(comm, "L_2 Error: %g\n", error);
// /* Check residual */
// SNESComputeFunction(snes, u, r);
// PetscPrintf(comm, "Initial Residual\n");
// VecChop(r, 1.0e-10);
// VecView(r, PETSC_VIEWER_STDOUT_WORLD);
// VecNorm(r, NORM_2, &res);
// PetscPrintf(comm, "L_2 Residual: %g\n", res);
// /* Check Jacobian */
// {
// Vec          b;
// MatStructure flag;
// MatNullSpace nullSpace2;
// PetscBool    isNull;
//
// CreateNullSpaces(user.dm, 1, &nullSpace2);
// MatNullSpaceTest(nullSpace2, J, &isNull);
// if (!isNull) SETERRQ(comm, PETSC_ERR_PLIB, "The null space calculated for the system operator is invalid.");
// MatNullSpaceDestroy(&nullSpace2);
//
// SNESComputeJacobian(snes, u, &A, &A, &flag);
// VecDuplicate(u, &b);
// VecSet(r, 0.0);
// SNESComputeFunction(snes, r, b);
// MatMult(A, u, r);
// VecAXPY(r, 1.0, b);
// VecDestroy(&b);
// PetscPrintf(comm, "Au - b = Au + F(0)\n");
// VecChop(r, 1.0e-10);
// VecView(r, PETSC_VIEWER_STDOUT_WORLD);
// VecNorm(r, NORM_2, &res);
// PetscPrintf(comm, "Linear L_2 Residual: %g\n", res);
// }
// }
//
// if (user.runType == RUN_FULL) {
// PetscContainer c;
// PetscSection   section;
// Vec            sol;
// PetscViewer    viewer;
// const char     *name;
//
// PetscViewerCreate(comm, &viewer);
// PetscViewerSetType(viewer, PETSCVIEWERVTK);
// PetscViewerFileSetName(viewer, "ex31_sol.vtk");
// PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);
// DMGetLocalVector(user.dm, &sol);
// PetscObjectGetName((PetscObject) u, &name);
// PetscObjectSetName((PetscObject) sol, name);
// DMGlobalToLocalBegin(user.dm, u, INSERT_VALUES, sol);
// DMGlobalToLocalEnd(user.dm, u, INSERT_VALUES, sol);
// DMGetDefaultSection(user.dm, &section);
// PetscObjectReference((PetscObject) user.dm); /* Needed because viewer destroys the DM */
// PetscViewerVTKAddField(viewer, (PetscObject) user.dm, DMPlexVTKWriteAll, PETSC_VTK_POINT_FIELD, (PetscObject) sol);
// PetscObjectReference((PetscObject) sol); /* Needed because viewer destroys the Vec */
// PetscContainerCreate(comm, &c);
// PetscContainerSetPointer(c, section);
// PetscObjectCompose((PetscObject) sol, "section", (PetscObject) c);
// PetscContainerDestroy(&c);
// DMRestoreLocalVector(user.dm, &sol);
// PetscViewerDestroy(&viewer);
// }
//
// MatNullSpaceDestroy(&nullSpace);
// if (user.jacobianMF) {
// VecDestroy(&userJ.u);
// }
// if (A != J) {
// MatDestroy(&A);
// }
// MatDestroy(&J);
// VecDestroy(&u);
// VecDestroy(&r);
// SNESDestroy(&snes);
// DMDestroy(&user.dm);
// PetscFinalize();
// return 0;
// }

