#ifndef __femus_applications_00_tutorial_hpp__
#define __femus_applications_00_tutorial_hpp__


#include "Assemble_unknown.hpp"


using namespace femus;


 
const std::vector< Unknown >  systems__generate_list_of_scalar_unknowns_for_each_FE_family_lagrangian() {


    std::vector< FEFamily >     feFamily = {LAGRANGE, LAGRANGE, LAGRANGE};
    std::vector< FEOrder >       feOrder = {FIRST, SERENDIPITY, SECOND};
    std::vector< int >        time_order = {0, 0, 0};  //0 = steady, 2 = time-dependent
    std::vector< bool >   is_pde_unknown = {true, true, true};

    assert( feFamily.size() == feOrder.size());
    assert( feFamily.size() == is_pde_unknown.size());
    assert( feFamily.size() == time_order.size());

    std::vector< Unknown >  unknowns(feFamily.size());

    for (unsigned int fe = 0; fe < unknowns.size(); fe++) {

        std::ostringstream unk;
        unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
        unknowns[fe]._name           = unk.str();
        unknowns[fe]._fe_family      = feFamily[fe];
        unknowns[fe]._fe_order       = feOrder[fe];
        unknowns[fe]._time_order     = time_order[fe];
        unknowns[fe]._is_pde_unknown = is_pde_unknown[fe];

    }


    return unknowns;

}



const std::vector< Unknown >  systems__generate_list_of_scalar_unknowns_for_each_FE_family_all_fe() {


    std::vector< FEFamily >     feFamily = {LAGRANGE, LAGRANGE, LAGRANGE, DISCONTINUOUS_POLYNOMIAL, DISCONTINUOUS_POLYNOMIAL};
    std::vector< FEOrder >       feOrder = {FIRST, SERENDIPITY, SECOND, ZERO, FIRST};
    std::vector< int >        time_order = {0, 0, 0, 0, 0};  //0 = steady, 2 = time-dependent
    std::vector< bool >   is_pde_unknown = {true, true, true, true, true};

    assert( feFamily.size() == feOrder.size());
    assert( feFamily.size() == is_pde_unknown.size());
    assert( feFamily.size() == time_order.size());

    std::vector< Unknown >  unknowns(feFamily.size());

    for (unsigned int fe = 0; fe < unknowns.size(); fe++) {

        std::ostringstream unk;
        unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
        unknowns[fe]._name           = unk.str();
        unknowns[fe]._fe_family      = feFamily[fe];
        unknowns[fe]._fe_order       = feOrder[fe];
        unknowns[fe]._time_order     = time_order[fe];
        unknowns[fe]._is_pde_unknown = is_pde_unknown[fe];

    }


    return unknowns;

}





#endif
