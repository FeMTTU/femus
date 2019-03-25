// #include "FE_convergence.hpp"

#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "CurrentQuantity.hpp"


namespace femus {
    
 

    
 

  double FunctionIntegral (const uint vb, MultiLevelProblem & ml_prob, double (*pt2func)(double, const std::vector<double> ) )  {

  const uint mesh_vb = vb;
  
  
const uint Level  = ml_prob.GetMeshTwo()._NoLevels - 1;
const uint myproc = ml_prob.GetMeshTwo()._iproc;
  double time = 0.;
  
  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(Level);
  

  double integral = 0.;
  

//parallel sum
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];
  
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {
      
    CurrentElem       currelem(iel,myproc,Level,vb,NULL,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh); //element without equation
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));

    const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = ml_prob.GetMeshTwo().get_dim();
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

    currelem.SetDofobjConnCoords();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);

     
    for (uint qp = 0; qp < el_ngauss; qp++) {

       for (uint fe = 0; fe < QL; fe++)     {  
             currgp.SetPhiElDofsFEVB_g (fe,qp); 
             currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);
        }  
     
double  Jac_g=0.;
          if (vb==0)   Jac_g = currgp.JacVectVV_g(xyz);  //not xyz_refbox!      
     else if (vb==1)   Jac_g = currgp.JacVectBB_g(xyz);  //not xyz_refbox!      

   const double  wgt_g = ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);


 xyz.val_g();
double myval_g = pt2func(time,xyz._val_g); 

 
  integral += wgt_g*Jac_g*myval_g;

   
    }//gauss loop
     
    }//element loop
    
         std::cout << std::endl  << " ^^^^^^^^^^^^^^^^^L'integrale sul processore "<< myproc << " vale: " << integral << std::endl;

//     double weights_sum = 0.;
//     for (uint qp = 0; qp < el_ngauss; qp++)  weights_sum += ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
//        std::cout << std::endl << " ^^^^^^^^^^^^^^^^^ La somma dei pesi  vale: " << weights_sum << std::endl;

       double J=0.;
   #ifdef HAVE_MPI
//       MPI_Reduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );  //This one gives J only to processor 0 !
      MPI_Allreduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!
   #else
   J = integral;
   #endif
    
     std::cout << std::endl << " ^^^^^^^^^^^^^^^^^L'integrale totale vale: " << J << std::endl;

  return J;  
  
}
   
    


}  //end namespace

    
