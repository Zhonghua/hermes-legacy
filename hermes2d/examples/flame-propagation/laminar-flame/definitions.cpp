class InitialSolutionTemperature : public ExactSolutionScalar
{
public:
  InitialSolutionTemperature(Mesh* mesh, double x1) : ExactSolutionScalar(mesh), x1(x1) {};

  virtual scalar value (double x, double y) const 
  {
    return (x <= x1) ? 1.0 : exp(x1 - x); 
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(10); 
  }

  // Value.
  double x1;
};

class InitialSolutionConcentration : public ExactSolutionScalar
{
public:
  InitialSolutionConcentration(Mesh* mesh, double x1) : ExactSolutionScalar(mesh), x1(x1) {};

  virtual scalar value (double x, double y) const 
  {
    return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x)); 
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(10); 
  }

  // Value.
  double x1;
};

// definition of reaction rate omega
class DXDYFilterOmega : public DXDYFilter
{
public:
  DXDYFilterOmega(Hermes::vector<MeshFunction*> solutions) : DXDYFilter(solutions) {};
protected:

  virtual void filter_fn (int n, Hermes::vector<scalar *> values, Hermes::vector<scalar *> dx, Hermes::vector<scalar *> dy, 
                          scalar* rslt, scalar* rslt_dx, scalar* rslt_dy) 
  {
    for (int i = 0; i < n; i++) 
    {
      scalar t1 = values.at(0)[i] - 1.0;
      scalar t2 = t1 * beta;
      scalar t3 = 1.0 + t1 * alpha;
      scalar t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
      scalar t5 = (beta / (t3 * t3)) * values.at(1)[i];
      rslt[i] = t4 * values.at(1)[i];
      rslt_dx[i] = t4 * (dx.at(1)[i] + dx.at(0)[i] * t5);
      rslt_dy[i] = t4 * (dy.at(1)[i] + dy.at(0)[i] * t5);
    }
  }
};

// definition of reaction rate omega_dt
class DXDYFilterOmega_dt : public DXDYFilter
{
public:
  DXDYFilterOmega_dt(Hermes::vector<MeshFunction*> solutions) : DXDYFilter(solutions) {};
protected:

  virtual void filter_fn (int n, Hermes::vector<scalar *> values, Hermes::vector<scalar *> dx, Hermes::vector<scalar *> dy, 
                          scalar* rslt, scalar* rslt_dx, scalar* rslt_dy) 
  {
    for (int i = 0; i < n; i++)
    {
      scalar t1 = values.at(0)[i] - 1.0;
      scalar t2 = t1 * beta;
      scalar t3 = 1.0 + t1 * alpha;
      scalar t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
      scalar t5 = (beta / (t3 * t3));
      rslt[i] = t4 * t5 * values.at(1)[i];
      rslt_dx[i] = 0.0;
      rslt_dy[i] = 0.0; // not important
    }
  }
};

// definition of reaction rate omega_dc
class DXDYFilterOmega_dc : public DXDYFilter
{
public:
  DXDYFilterOmega_dc(Hermes::vector<MeshFunction*> solutions) : DXDYFilter(solutions) {};
protected:

  virtual void filter_fn (int n, Hermes::vector<scalar *> values, Hermes::vector<scalar *> dx, Hermes::vector<scalar *> dy, 
                          scalar* rslt, scalar* rslt_dx, scalar* rslt_dy) 
  {
    for (int i = 0; i < n; i++)
    {
      scalar t1 = values.at(0)[i] - 1.0;
      scalar t2 = t1 * beta;
      scalar t3 = 1.0 + t1 * alpha;
      scalar t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
      rslt[i] = t4;
      rslt_dx[i] = 0.0;
      rslt_dy[i] = 0.0; // not important
    }
  }
};

// weak forms for the Newton's method

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(Hermes::vector<std::string> neumann_boundaries, Hermes::vector<std::string> newton_boundaries,
                 double time_step, double Le, double kappa, 
                 DXDYFilterOmega_dt* omega_dt, DXDYFilterOmega_dc* omega_dc, DXDYFilterOmega* omega, 
                 InitialSolutionTemperature* t_prev_time_1, InitialSolutionTemperature* t_prev_time_2, InitialSolutionTemperature* t_prev_newton,
                 InitialSolutionConcentration* c_prev_time_1, InitialSolutionConcentration* c_prev_time_2, InitialSolutionConcentration* c_prev_newton, 
                 bool JFNK = false) : WeakForm(2) 
  {

    // Jacobian forms - volumetric.
    add_matrix_form(new CustomMatrixForm_0_0(0, 0, time_step));
    add_matrix_form(new CustomMatrixForm_0_1(0, 1));
    add_matrix_form(new CustomMatrixForm_1_0(1, 0));
    add_matrix_form(new CustomMatrixForm_1_1(1, 1, time_step, Le));

    // Jacobian forms - surface.
    add_matrix_form_surf(new CustomMatrixForm_0_0_surf(0, 0, newton_boundaries, kappa));

    // Residual forms - volumetric.
    add_vector_form(new CustomVectorForm_0_0(0, time_step));
    add_vector_form(new CustomVectorForm_0_1(1));
    add_vector_form(new CustomVectorForm_1_0(0));
    add_vector_form(new CustomVectorForm_1_1(1, time_step, Le));

    add_vector_form(new CustomVectorFormVol_0(0, time_step));
    add_vector_form(new CustomVectorFormVol_1(1, time_step, Le));

    // Residual forms - surface.
    add_vector_form_surf(new CustomVectorForm_0_surf(0, newton_boundaries, kappa));
  };

  ~CustomWeakForm() {};

private:

  // Jacobian part.
  class CustomMatrixForm_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixForm_0_0(int i, int j, double time_step) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), time_step(time_step) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      Func<scalar>* domegadt = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (  1.5 * u->val[i] * v->val[i] / time_step
                        +  u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]
                        - domegadt->val[i] * u->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

    double time_step;
  };

  class CustomMatrixForm_0_0_surf : public WeakForm::MatrixFormSurf
  {
  public:
    CustomMatrixForm_0_0_surf(int i, int j, Hermes::vector<std::string> newton_boundaries, double kappa) : WeakForm::MatrixFormSurf(i, j, newton_boundaries), kappa(kappa) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (kappa * u->val[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

    double kappa;
  };

  class CustomMatrixForm_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixForm_0_1(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      Func<scalar>* domegady = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (- domegady->val[i] * u->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

  };

  class CustomMatrixForm_1_0 : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixForm_1_0(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      Func<scalar>* domegadt = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * ( domegadt->val[i] * u->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

  };

  class CustomMatrixForm_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixForm_1_1(int i, int j, double time_step, double Le) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), time_step(time_step), Le(Le) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
    {
      scalar result = 0;
      Func<scalar>* domegady = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (  1.5 * u->val[i] * v->val[i] / time_step
                          +  (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) / Le
                          + domegady->val[i] * u->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

  double time_step, Le;

  };

  // Residual part.

  class CustomVectorForm_0_0 : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorForm_0_0(int i, double time_step) : WeakForm::VectorFormVol(i), time_step(time_step) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,  
                         Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      Func<scalar>* domegadt = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (  1.5 * u_ext[0]->val[i] * v->val[i] / time_step
                        +  u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]
                        - domegadt->val[i] * u_ext[0]->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

    double time_step;
  };

  class CustomVectorForm_0_1 : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorForm_0_1(int i) : WeakForm::VectorFormVol(i) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,  
                         Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      Func<scalar>* domegady = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (- domegady->val[i] * u_ext[1]->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

  };

  class CustomVectorForm_1_0 : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorForm_1_0(int i) : WeakForm::VectorFormVol(i) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,  
                         Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      Func<scalar>* domegadt = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * ( domegadt->val[i] * u_ext[0]->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

  };

  class CustomVectorForm_1_1 : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorForm_1_1(int i, double time_step, double Le) : WeakForm::VectorFormVol(i), time_step(time_step), Le(Le) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,  
                         Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      Func<scalar>* domegady = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (  1.5 * u_ext[1]->val[i] * v->val[i] / time_step
                          +  (u_ext[1]->dx[i] * v->dx[i] + u_ext[1]->dy[i] * v->dy[i]) / Le
                          + domegady->val[i] * u_ext[1]->val[i] * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }

  double time_step, Le;

  };

  class CustomVectorFormVol_0 : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol_0(int i, double time_step) : WeakForm::VectorFormVol(i), time_step(time_step) {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* titer = u_ext[0];
      Func<scalar>* t_prev_time_1 = ext->fn[0];
      Func<scalar>* t_prev_time_2 = ext->fn[1];
      Func<scalar>* omega = ext->fn[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ( (3.0 * titer->val[i] - 4.0 * t_prev_time_1->val[i] 
                             + t_prev_time_2->val[i]) * v->val[i] / (2.0 * time_step) +
                            (titer->dx[i] * v->dx[i] + titer->dy[i] * v->dy[i]) -
                            omega->val[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the test function and solution plus two.
      return Ord(10);
    }

  private:
    double time_step;
  };

  class CustomVectorFormVol_1 : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol_1(int i, double time_step, double Le) : WeakForm::VectorFormVol(i),
    time_step(time_step), Le(Le)  {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* c_prev_newton = u_ext[1];
      Func<scalar>* c_prev_time_1 = ext->fn[0];
      Func<scalar>* c_prev_time_2 = ext->fn[1];
      Func<scalar>* omega = ext->fn[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ( (3.0 * c_prev_newton->val[i] - 4.0 * c_prev_time_1->val[i] + c_prev_time_2->val[i])
                        * v->val[i] / (2.0 * time_step) +
                        (c_prev_newton->dx[i] * v->dx[i] + c_prev_newton->dy[i] * v->dy[i]) / Le +
                        omega->val[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the test function and solution plus two.
      return Ord(10);
    }

  private:
    double time_step, Le;
  };


  class CustomVectorForm_0_surf : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorForm_0_surf(int i, Hermes::vector<std::string> newton_boundaries, double kappa) : WeakForm::VectorFormSurf(i, newton_boundaries), 
    kappa(kappa)  {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      Func<scalar>* t_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (kappa * t_prev_newton->val[i] * v->val[i]);
      return result;  
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the test function and solution plus two.
      return Ord(10);
    }

  private:
    double kappa;
  };

protected:

  double time_step; 
  double Le;
  double kappa;
  DXDYFilterOmega_dt* omega_dt;
  DXDYFilterOmega_dc* omega_dc;
  DXDYFilterOmega* omega;
  InitialSolutionTemperature* t_prev_time_1;
  InitialSolutionTemperature* t_prev_time_2;
  InitialSolutionTemperature* t_prev_newton;
  InitialSolutionConcentration* c_prev_time_1; 
  InitialSolutionConcentration* c_prev_time_2; 
  InitialSolutionConcentration* c_prev_newton; 
  bool JFNK;

};
