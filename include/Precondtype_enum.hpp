#ifndef _precondtypeM_
#define _precondtypeM_

  enum PreconditionerTypeM {IDENTITY_PRECONDM =0,
			   JACOBI_PRECONDM,
			   BLOCK_JACOBI_PRECONDM,
			   SOR_PRECONDM,
			   SSOR_PRECONDM,
			   EISENSTAT_PRECONDM,
			   ASM_PRECONDM,
			   CHOLESKY_PRECONDM,
			   ICC_PRECONDM,
			   ILU_PRECONDM,
			   LU_PRECONDM,
			   USER_PRECONDM,
			   SHELL_PRECONDM,
                           AMG_PRECONDM,
			   INVALID_PRECONDITIONERM
  };



#endif