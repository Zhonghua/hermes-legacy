Time-Dependent Problems (09-timedep-rk)
----------------------------------

**Git reference:** Tutorial example `09-timedep-rk <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/09-timedep-rk>`_. 

Model problem
~~~~~~~~~~~~~

This example is a continuation of the example "09-timedep-basic" and it shows how 
to perform time integration with Runge-Kutta methods using arbitrary Butcher's 
tables. Currently (as of January 2011) approx. 20 tables are available by default,
as can be seen below. They are taken from various sources including J. Butcher's
book and the Wikipedia page http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods. 
If you know about some other interesting R-K methods that are missing in our database,
please let us know!

.. image:: 09/vitus1.png
   :align: center
   :width: 400
   :height: 500
   :alt: Model geometry and temperature distribution after 24 hours.

We will solve the standard heat transfer equation

.. math::
    :label: eqvit1

       c \varrho\frac{\partial T}{\partial t} - \lambda \Delta T = 0

equipped with a Dirichlet condition

.. math::

     T = T_{init}

on the bottom edge $\Gamma_{ground}$ and a Newton condition

.. math::

     \frac{\partial T}{\partial \nu} = \alpha(T_{ext}(t) - T)

on the rest of the boundary $\Gamma_{air}$. Here, $c$ is the heat capacity of the material,
$\varrho$ the material density, $\lambda$ the thermal conductivity,
$T_{init}$ the fixed temperature on the
ground (same as the initial temperature of the building), and $\alpha$
the heat transfer coefficient 
between the building and the surrounding air. The surrounding air temperature
$T_{ext}$ is time-dependent of the form

.. math::

     T_{ext}(t) = T_{init} + 10\sin(2\pi t/T_{final}),

where $T_{final}$ is 24 hours (translated into seconds).

Equation :eq:`eqvit1` is equipped with an initial condition of the
form

.. math::

     T(x,y,0) = T_{init}(x,y) \ \ \ \mbox{in} \ \Omega.

For simplicity we will use the implicit Euler method with a constant
time step $\tau$, which transforms equation :eq:`eqvit1` into

.. math::

     c \varrho\frac{T^{n+1} - T^n}{\tau} - \lambda \Delta T^{n+1} = 0.

Butcher's table
~~~~~~~~~~~~~~~

Butcher's tables type. The last number in the name always means order,
the one before last (if provided) is the number of stages::

    enum ButcherTableType
    {
       Explicit_RK_1,               // Explicit Runge-Kutta RK-1, or explicit Euler method.
       Implicit_RK_1,               // Implicit Runge-Kutta RK-1, or implicit Euler method.
       Explicit_RK_2,               // Explicit Runge-Kutta RK-2 method.
       Implicit_Crank_Nicolson_2_2, // Implicit Crank_Nicolson method.
       Implicit_SIRK_2_2,           // Implicit SIRK-2-2 method.
       Implicit_ESIRK_2_2,          // Implicit ESIRK-2-2 method.
       Implicit_SDIRK_2_2,          // Implicit SDIRK-2-2 method.
       Implicit_Lobatto_IIIA_2_2,   // Implicit Lobatto IIIA-2 method.
       Implicit_Lobatto_IIIB_2_2,   // Implicit Lobatto IIIB-2 method.
       Implicit_Lobatto_IIIC_2_2,   // Implicit Lobatto IIIB-2 method.
       Explicit_RK_3,               // Explicit Runge-Kutta RK-3 method.
       Explicit_RK_4,               // Explicit Runge-Kutta RK-4 method.
       Implicit_Lobatto_IIIA_3_4,   // Implicit Lobatto IIIA-4 method.
       Implicit_Lobatto_IIIB_3_4,   // Implicit Lobatto IIIB-4 method.
       Implicit_Lobatto_IIIC_3_4,   // Implicit Lobatto IIIB-4 method.
       Implicit_Radau_IIA_3_5,      // Implicit Radau IIA-5 method.
       Implicit_SDIRK_4_5,          // Implicit SDIRK-4-5 method.
       Implicit_DIRK_7_45_embedded  // Implicit embedded DIRK method pair of orders four in five (from the paper 
                                    // Fudziah Ismail et all: Embedded Pair of Diagonally Implicit Runge-Kutta  
                                    // Method for Solving Ordinary Differential Equations). The method has
                                    // 7 stages but the first one is explicit.
    };

Choose one Butcher's table, or define your own Butcher's table::

    ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

    ButcherTable bt(butcher_table_type);
    if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
    if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
    if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

Weak formulation
~~~~~~~~~~~~~~~~

The corresponding weak formulation is

.. math::

     \int_{\Omega} c \varrho\frac{T^{n+1}}{\tau}v + \int_{\Omega} \lambda \nabla T^{n+1}\cdot \nabla v + \int_{\Gamma_{air}} \alpha \lambda T^{n+1}v = \int_{\Omega} c \varrho\frac{T^{n}}{\tau}v + \int_{\Gamma_{air}} \alpha \lambda T_{ext}(t^{n+1})v.

In this example we use string boundary markers::

    // Boundary markers.
    const std::string BDY_GROUND = "Boundary ground";
    const std::string BDY_AIR = "Boundary air";

Boundary condition types are defined using the BCTypes class::

    // Enter boundary markers.
    BCTypes bc_types;
    bc_types.add_bc_dirichlet(BDY_GROUND);
    bc_types.add_bc_newton(BDY_AIR);

Values for Dirichlet boundary conditions are set via the BCValues class::

    // Enter Dirichlet boundary values.
    BCValues bc_values;
    bc_values.add_const(BDY_GROUND, TEMP_INIT);

Then the space for the temperature $T$ is set up::

    // Initialize an H1 space with default shepeset.
    H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
    int ndof = Space::get_num_dofs(&space);
    info("ndof = %d.", ndof);

Defining weak forms and accessing external functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bilinear and linear forms are defined as follows::

    template<typename Real, typename Scalar>
    Scalar stac_jacobian_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += -wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      }
    
      return result * LAMBDA / HEATCAP / RHO;
    }

    template<typename Real, typename Scalar>
    Scalar stac_residual_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Func<Scalar>* u_prev = u_ext[0];

      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += -wt[i] * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i]);          
      }

      return result * LAMBDA / HEATCAP / RHO;
    }

    template<typename Real, typename Scalar>
    Scalar stac_jacobian_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return - LAMBDA / HEATCAP / RHO * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    template<typename Real, typename Scalar>
    Scalar stac_residual_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Func<Scalar>* u_prev = u_ext[0];

      // This is a temporary workaround. The stage time t_n + h * c_i
      // can be accessed via u_stage_time->val[0];
      Func<Scalar>* u_stage_time = ext->fn[0]; 
  
      Scalar stage_time = u_stage_time->val[0];
      Real stage_ext_temp = temp_ext<Real>(stage_time);

      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * (stage_ext_temp - u_prev->val[i]) * v->val[i];                   
      }

      return LAMBDA / HEATCAP / RHO * ALPHA * result;
    }  

Notice how the previous time level temperature is accessed:

::

      Func<Real> *u_stage_time = ext->fn[0];
    
Setting initial condition
~~~~~~~~~~~~~~~~~~~~~~~~~ 

Next we need to initialize the previous time level solution u_prev_time with the initial condition $T_{init}$.
Besides holding the finite element solution, the Solution class
can be forced to return zero, to return a constant, or to return an arbitrary function
using the methods set_zero(), set_const() and set_exact(), respectively.
Here we simply call set_const() and supply the initial temperature::

    // Set constant initial condition.
    Solution u_prev_time(&mesh, TEMP_INIT);

Registering external functions in weak forms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The weak forms are registered as follows::

    // Initialize weak formulation.
    WeakForm wf;
    wf.add_matrix_form(callback(stac_jacobian_vol));
    wf.add_matrix_form_surf(callback(stac_jacobian_surf), BDY_AIR);
    wf.add_vector_form(callback(stac_residual_vol));
    wf.add_vector_form_surf(callback(stac_residual_surf), BDY_AIR);

As opposed to the previous example where the time-discretization was hardwired
into the weak formulation, now we only need the weak formulation of the right-hand side.

Initializing the discrete problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, the DiscreteProblem class is initialized::

    // Initialize the FE problem.
    bool is_linear = false;
    DiscreteProblem dp(&wf, &space, is_linear);

Start the time stepping
~~~~~~~~~~~~~~~~~~~~~~~

We are now ready to start the time stepping. For completeness, we show 
the entire time stepping loop below::

    // Time stepping loop:
    double current_time = 0.0; int ts = 1;
    do 
    {
      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).", 
           current_time, time_step, bt.get_size());
      bool verbose = true;
      bool is_linear = true;
      if (!rk_time_step(current_time, time_step, &bt, coeff_vec, &dp, matrix_solver, verbose, is_linear)) {
        error("Runge-Kutta time step failed, try to decrease time step size.");
      }

      // Convert coeff_vec into a new time level solution.
      Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);

      // Update time.
      current_time += time_step;

      // Show the new time level solution.
      char title[100];
      sprintf(title, "Time %3.2f, exterior temperature %3.5f", current_time, temp_ext(current_time));
      Tview.set_title(title);
      Tview.show(&u_prev_time);

      // Increase counter of time steps.
      ts++;
    } 
    while (current_time < T_FINAL);

