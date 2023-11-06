/*=========================================================================

 Program: FEMuS
 Module: NonLinearImplicitSystem
 Authors: Simone Bnà

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_MonolithicFSINonLinearImplicitSystem_hpp__
#define __femus_equations_MonolithicFSINonLinearImplicitSystem_hpp__

//----------------------------------------------------------------------------
// inclarudes :
//----------------------------------------------------------------------------
#include "NonLinearImplicitSystem.hpp"


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------

/**
 * The non linear implicit system class for MonolithicFSI NonLinearImplicitSystem
 */

class MonolithicFSINonLinearImplicitSystem : public NonLinearImplicitSystem {

public:

    /** Constructor.  Optionally initializes required data structures. */
    MonolithicFSINonLinearImplicitSystem (MultiLevelProblem& ml_probl, const std::string& name, const unsigned int number, const LinearEquationSolverType & smoother_type);

    /** Destructor */
    virtual ~MonolithicFSINonLinearImplicitSystem();

    /** The type of the parent. */
    typedef NonLinearImplicitSystem Parent;

    /** Init the system PDE structures */
    virtual void init();

    /**
     * @returns \p "MonolithicFSINonlinearImplicit".  Helps in identifying
     * the system type in an equation system file.
    */
    virtual std::string system_type () const {
        return "MonolithicFSINonlinearImplicit";
    }


    void SetElementBlockNumberFluid(unsigned const &dim_vanka_block, unsigned const &overlap=0);
    void SetElementBlockNumberSolid(unsigned const &dim_vanka_block, unsigned const &overlap=0);
    void SetElementBlockNumberPorous(unsigned const &dim_vanka_block, unsigned const &overlap=0);
    void SetElementBlockFluidAll();
    void SetElementBlockSolidAll();
    void SetElementBlockPorousAll();

protected:

    /** Create the Prolongator and Restrictor Operators for the Algebraic Multigrid Solver */
    void BuildProlongatorMatrix(unsigned gridf);

    void BuildAmrProlongatorMatrix(unsigned level);

      /** Only for FSI */
      void Build_RestrictionTranspose_OneElement_OneFEFamily_With_Pair_In_System(const LinearEquation& lspdef,
                                     const LinearEquation& lspdec,
                                     const int& ielc,
                                     SparseMatrix* Projmat,
                                     const unsigned& index_sol,
                                     const unsigned& kkindex_sol,
                                     const unsigned& index_pair_sol,
                                     const unsigned& kkindex_pair_sol,
                                     const elem_type * elem_type_in) const;
                                     

private:

};



} //end namespace femus



#endif
