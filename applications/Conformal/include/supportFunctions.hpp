 void GetElementNearVertexNumber(MultiLevelSolution &mlSol){

  unsigned level = mlSol.GetMLMesh()->GetNumberOfLevels() - 1;

  Solution* sol = mlSol.GetSolutionLevel (level);
  Mesh* msh = mlSol.GetMLMesh()->GetLevel (level);

  unsigned DIM = 3u;
  unsigned solIndex = mlSol.GetIndex ("ENVN");
  
  unsigned solType = mlSol.GetSolutionType (solIndex);

  sol->_Sol[solIndex]->zero();
  
  unsigned iproc = msh->processor_id();
  
  // ELEMENT LOOP: each process loops only on the elements that it owns.
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    // Number of solution element dofs.
    short unsigned ielGeom = msh->GetElementType (iel);
    unsigned nDofs  = msh->GetElementDofNumber (iel, solType);
    
    for (unsigned i = 0; i < nDofs; i++) {

      unsigned iDof = msh->GetSolutionDof (i, iel, solType);

      sol->_Sol[solIndex]->add(iDof, 1);
      
    }
  }
   
   sol->_Sol[solIndex]->close();
  
}

// Eugenio's standard FEMuS function.
void CopyDisplacement (MultiLevelSolution &mlSol,  const bool &forward) {

  unsigned level = mlSol.GetMLMesh()->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel (level);
  Mesh* msh = mlSol.GetMLMesh()->GetLevel (level);

  unsigned DIM = 3u;
  std::vector < unsigned > solDxIndex (DIM);
  solDxIndex[0] = mlSol.GetIndex ("Dx1");
  solDxIndex[1] = mlSol.GetIndex ("Dx2");
  solDxIndex[2] = mlSol.GetIndex ("Dx3");

  std::vector < unsigned > solNDxIndex (DIM);
  solNDxIndex[0] = mlSol.GetIndex ("nDx1");
  solNDxIndex[1] = mlSol.GetIndex ("nDx2");
  solNDxIndex[2] = mlSol.GetIndex ("nDx3");

  if (forward) {
    for (unsigned i = 0; i < DIM; i++) {
      * (solution->_Sol[solNDxIndex[i]]) = * (solution->_Sol[solDxIndex[i]]);
    }
  }
  else {
    for (unsigned i = 0; i < DIM; i++) {
      * (solution->_Sol[solDxIndex[i]]) = * (solution->_Sol[solNDxIndex[i]]);
    }
  }
}

double max (const double &a , const double &b) {
  return (a > b) ? a : b;
}


void ChangeTriangleConfiguration1 (const std::vector<unsigned> & ENVN, std::vector <double> &angle) {
  double scale;
  if (ENVN[0] < ENVN[1] && ENVN[0] < ENVN[2]) {
    scale = (M_PI - angle[0]) / (angle[1] + angle [2]);
    angle[1] *= scale;
    angle[2] *= scale;
  }
  else if (ENVN[0] < ENVN[1] && ENVN[0] == ENVN[2]) {
    angle[1] = M_PI - 2. * angle[0];
  }
  else if (ENVN[0] <= ENVN[1]  && ENVN[0] > ENVN[2]) {
    scale = (M_PI - angle[2]) / (angle[1] + angle [0]);
    angle[1] *= scale;
    angle[0] *= scale;
  }
  else if (ENVN[0] == ENVN[1] && ENVN[0] < ENVN[2]) {
    angle[2] = M_PI - 2. * angle[0];
  }
  else if (ENVN[0] == ENVN[1] && ENVN[0] == ENVN[2]) {
    angle[0] = angle[1] = angle[2] =  M_PI / 3.;
  }
  else if (ENVN[0] > ENVN[1] && ENVN[0] <= ENVN[2]) {
    scale = (M_PI - angle[1]) / (angle[0] + angle [2]);
    angle[0] *= scale;
    angle[2] *= scale;
  }
  else if (ENVN[0] > ENVN[1] && ENVN[0] > ENVN[2]) {
    if (ENVN[1] < ENVN[2]) {
      scale = (M_PI - angle[1]) / (angle[0] + angle [2]);
      angle[0] *= scale;
      angle[2] *= scale;
    }
    else if (ENVN[1] == ENVN[2]) {
      angle[0] = M_PI - 2. * angle[1];
    }
    else if (ENVN[1] > ENVN[2]) {
      scale = (M_PI - angle[2]) / (angle[0] + angle [1]);
      angle[0] *= scale;
      angle[1] *= scale;
    }
  }
}


void ChangeTriangleConfiguration2 (const std::vector<unsigned> & ENVN, std::vector <double> &angle) {
  unsigned type = 3; // there are 2 or 3 leading angles
  if (ENVN[0] < ENVN[1]) { // 0 leads on 1
    if (ENVN[0] < ENVN[2]) type = 0; // 0 is leading angle
    else if (ENVN[0] > ENVN[2]) type = 2; // 2 is leading angle
  }
  else if (ENVN[0] > ENVN[1]) { // 1 leads on 0
    if (ENVN[1] < ENVN[2]) type = 1; // 1 is leading angle
    else if (ENVN[1] > ENVN[2]) type = 2; // 2 is leading angle
  }
  else { // 0 equals 1
    if (ENVN[0] > ENVN[2]) type = 2; // 2 is leading angle
  }

  double scale;
  if (type == 0) {
    scale = (M_PI - angle[0]) / (angle[1] + angle [2]);
    angle[1] *= scale;
    angle[2] *= scale;
  }
  else if (type == 1) {
    scale = (M_PI - angle[1]) / (angle[0] + angle [2]);
    angle[0] *= scale;
    angle[2] *= scale;
  }
  else if (type == 2) {
    scale = (M_PI - angle[2]) / (angle[1] + angle [0]);
    angle[1] *= scale;
    angle[0] *= scale;
  }
  else {
    scale = M_PI / (angle[0] + angle[1] + angle[2]);
    angle[0] *= scale;
    angle[1] *= scale;
    angle[2] *= scale;
  }
}


