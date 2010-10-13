#include "hermes3d.h"


void OGProjection::project_internal(Tuple<Space *> spaces, WeakForm *proj_wf, scalar* target_vec, MatrixSolverType matrix_solver)
{
  _F_
  int n = spaces.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_internal().");
  for (int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_internal().", i);
  if (spaces.size() != n) error("Number of spaces must matchnumber of projected functions in project_internal().");

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = Space::assign_dofs(spaces);

  // Initialize FeProblem.
  bool is_linear = true;
  DiscreteProblem* dp = new DiscreteProblem(proj_wf, spaces, is_linear);

  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  dp->assemble(matrix, rhs);

  // Calculate the coefficient vector.
  bool solved = solver->solve();
  scalar* coeffs;
  if (solved) 
    coeffs = solver->get_solution();

  if (target_vec != NULL) 
    for (int i=0; i<ndof; i++) target_vec[i] = coeffs[i];
    
  delete solver;
  delete matrix;
  delete rhs;
  delete dp;
  delete proj_wf;
}


void OGProjection::project_global(Tuple<Space *> spaces, Tuple<MeshFunction *> source_meshfns, 
                              scalar* target_vec, MatrixSolverType matrix_solver, Tuple<ProjNormType> proj_norms)
{
  _F_
  int n = spaces.size();  

  // define temporary projection weak form
  WeakForm* proj_wf = new WeakForm(n);
  int found[100];
  for (int i = 0; i < 100; i++) found[i] = 0;
  for (int i = 0; i < n; i++) 
  {
    int norm;
    if (proj_norms == Tuple<ProjNormType>()) norm = HERMES_DEFAULT_PROJ_NORM;
    else norm = proj_norms[i];
    if (norm == HERMES_L2_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, L2projection_biform<double, scalar>, L2projection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, L2projection_liform<double, scalar>, L2projection_liform<Ord, Ord>,
                               HERMES_ANY, source_meshfns[i]);
    }
    if (norm == HERMES_H1_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, H1projection_biform<double, scalar>, H1projection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, H1projection_liform<double, scalar>, H1projection_liform<Ord, Ord>,
                               HERMES_ANY, source_meshfns[i]);
    }
    if (norm == HERMES_H1_SEMINORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, H1_semi_projection_biform<double, scalar>, H1_semi_projection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, H1_semi_projection_liform<double, scalar>, H1_semi_projection_liform<Ord, Ord>,
                               HERMES_ANY, source_meshfns[i]);
    }
    if (norm == HERMES_HCURL_NORM) 
    {
      found[i] = 1;
      proj_wf->add_matrix_form(i, i, Hcurlprojection_biform<double, scalar>, Hcurlprojection_biform<Ord, Ord>);
      proj_wf->add_vector_form(i, Hcurlprojection_liform<double, scalar>, Hcurlprojection_liform<Ord, Ord>,
                               HERMES_ANY, source_meshfns[i]);
    }
  }
  for (int i=0; i < n; i++) 
  {
    if (found[i] == 0) 
    {
      printf("index of component: %d\n", i);
      error("Wrong projection norm in project_global().");
    }
  }

  project_internal(spaces, proj_wf, target_vec, matrix_solver);
}

void OGProjection::project_global(Tuple<Space *> spaces, 
                              Tuple<Solution*> sols_src, Tuple<Solution*> sols_dest, 
                              MatrixSolverType matrix_solver, Tuple<ProjNormType> proj_norms)
{
  _F_
  scalar* target_vec = new scalar[Space::get_num_dofs(spaces)];
  Tuple<MeshFunction *> ref_slns_mf;
  for (int i = 0; i < sols_src.size(); i++) 
    ref_slns_mf.push_back(static_cast<MeshFunction*>(sols_src[i]));
  
  OGProjection::project_global(spaces, ref_slns_mf, target_vec, matrix_solver, proj_norms);
  
  Solution::vector_to_solutions(target_vec, spaces, sols_dest);
  
  delete [] target_vec;
}
