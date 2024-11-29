// Copyright (c) 2011-2012, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// APPLICATION INCLUDES
#include <ActiveContraction.h>
#include <BoundaryConditions.h>
#include <MechanicsModel.h>

#include <math.h>

// IBAMR INCLUDES
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/quadrature.h>
#include <libmesh/tensor_tools.h>


// STATIC VARIABLES
int ActiveContraction::act_system_num, ActiveContraction::T_system_num;
namespace  // private namespace
{
// Model constants.
static const double a = 0.35;            // dimensionless
static const double A1 = -29.0;          // dimensionless
static const double A2 = 138.0;          // dimensionless
static const double A3 = 129.0;          // dimensionless
static const double alpha_0 = 8.0;       // sec^-1
static const double alpha_1 = 30.0;      // sec^-1
static const double alpha_2 = 130.0;     // sec^-1
static const double alpha_3 = 625.0;     // sec^-1
static const double alpha_r1 = 2.0*1.0;      // sec^-1
static const double alpha_r2 = 1.75*1.0;     // sec^-1
static const double beta_0 = 4.9;        // dimensionless
static const double beta_1 = -4.0;       // dimensionless
static const double Ca_50_ref = 1.05;    // uM
static const double Ca_TRPN_max = 70.0;  // uM, LYF modification original 70
static const double gamma_trpn = 2.0;    // dimensionless
static const double k_on = 100.0;        // uM^-1 sec^-1, LYF modification original 100
static const double k_refoff = 200.0;    // sec^-1,	LYF modification original 200
static const double K_Z = 0.15;          // dimensionless
static const double n = 3.0;             // dimensionless
static const double n_r = 3.0;           // dimensionless
static const double T_ref = 56.2;        // kPa = 1000 N m^-2 = 10000 dyne cm^-2
static const double z_p = 0.85;          // dimensionless
}

// CLASS IMPLEMENTATION

void
ActiveContraction::update_active_tension_model_state_variables(
    EquationSystems* equation_systems,
    const double time,
    const double dt,
	int flag)
{
	if(flag==0)
	{
		const MeshBase& mesh = equation_systems->get_mesh();
		const int dim = mesh.mesh_dimension();
		std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, NDIM, CONSTANT);

		System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
		const DofMap& X_dof_map = X_system.get_dof_map();
	#ifdef DEBUG_CHECK_ASSERTIONS
		for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
	#endif
		//blitz::Array<std::vector<unsigned int>,1> X_dof_indices(NDIM);
		std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
	 //   for (unsigned int d = 0; d < NDIM; ++d) X_dof_indices(d).reserve(NDIM == 2 ? 9 : 27);
		std::unique_ptr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
		X_fe->attach_quadrature_rule(qrule.get());
		const std::vector<std::vector<VectorValue<double> > >& dphi_X = X_fe->get_dphi();

		System& U_system = equation_systems->get_system<System>(IBFEMethod::VELOCITY_SYSTEM_NAME);
	#ifdef DEBUG_CHECK_ASSERTIONS
		const DofMap& U_dof_map = U_system.get_dof_map();
		for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(U_dof_map.variable_type(d) == U_dof_map.variable_type(0));
		for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(U_dof_map.variable_type(d) == X_dof_map.variable_type(0));
	#endif

		System& f0_system = equation_systems->get_system<System>(MechanicsModel::f0_mv_system_num);
		const DofMap& f0_dof_map = f0_system.get_dof_map();
	#ifdef DEBUG_CHECK_ASSERTIONS
		for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(f0_dof_map.variable_type(d) == f0_dof_map.variable_type(0));
	#endif
		//blitz::Array<std::vector<unsigned int>,1> f0_dof_indices(NDIM);
		std::vector<std::vector<unsigned int> > f0_dof_indices(NDIM);

		System& T_system = equation_systems->get_system<System>(T_system_num);
		const DofMap& T_dof_map = T_system.get_dof_map();
		std::vector<unsigned int> T_dof_indices;

		System& act_system = equation_systems->get_system<System>(act_system_num);
		const DofMap& act_dof_map = act_system.get_dof_map();
		//blitz::TinyVector<std::vector<unsigned int>,NUM_ACT_VARS> act_dof_indices;
		const int NUM_ACT_VARS=6;
		std::vector<std::vector<unsigned int> > act_dof_indices(NUM_ACT_VARS);

		X_system.solution->localize(*X_system.current_local_solution);
		NumericVector<double>& X_data = *(X_system.current_local_solution);
		X_data.close();

		U_system.solution->localize(*U_system.current_local_solution);
		NumericVector<double>& U_data = *(U_system.current_local_solution);
		U_data.close();

		NumericVector<double>&  f0_data = *( f0_system.solution);
		NumericVector<double>&   T_data = *(  T_system.solution);
		NumericVector<double>& act_data = *(act_system.solution);

		TensorValue<double> FF, dFF_dt;
		//blitz::Array<double,2> X_node, U_node;
		boost::multi_array<double,2> X_node, U_node;
		const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
		const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
		
		double Ca_i_max = -1.0e300;
		double Ca_i_min = +1.0e300;
		double Ca_b_max = -1.0e300;
		double Ca_b_min = +1.0e300;
		double T_max = -1.0e300;
		for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
		{
			Elem* const elem = *el_it;

			X_fe->reinit(elem);
			for (unsigned int d = 0; d < NDIM; ++d)
			{
				X_dof_map.dof_indices(elem, X_dof_indices[d], d);
			}

			for (unsigned int d = 0; d < NDIM; ++d)
			{
				f0_dof_map.dof_indices(elem, f0_dof_indices[d], d);
			}

			T_dof_map.dof_indices(elem, T_dof_indices, 0);

			for (unsigned int d = 0; d < NUM_ACT_VARS; ++d)
			{
				act_dof_map.dof_indices(elem, act_dof_indices[d], d);
			}

	#ifdef DEBUG_CHECK_ASSERTIONS
			const unsigned int n_qp = qrule->n_points();
			TBOX_ASSERT(n_qp == 1);
	#endif
			const unsigned int qp = 0;

			get_values_for_interpolation(X_node, X_data, X_dof_indices);
			jacobian(FF,qp,X_node,dphi_X);

			get_values_for_interpolation(U_node, U_data, X_dof_indices);
			jacobian(dFF_dt,qp,U_node,dphi_X);

			VectorValue<double> f0;
			for (unsigned int d = 0; d < NDIM; ++d)
			{
				f0(d) = f0_data(f0_dof_indices[d][0]);  // piecewise constant representation
			}
			const VectorValue<double> f = FF*f0;

			const double lambda = f.norm();
			const double dlambda_dt = (0.5/lambda)*f0*((FF.transpose()*dFF_dt + dFF_dt.transpose()*FF)*f0);

			double Ca_i = act_data(act_dof_indices[CA_I_IDX][0]);
			double Ca_b = act_data(act_dof_indices[CA_B_IDX][0]);
			double Q1   = act_data(act_dof_indices[  Q1_IDX][0]);
			double Q2   = act_data(act_dof_indices[  Q2_IDX][0]);
			double Q3   = act_data(act_dof_indices[  Q3_IDX][0]);
			double z    = act_data(act_dof_indices[   Z_IDX][0]);
			
			if (elem->subdomain_id() == 21||elem->subdomain_id() == 22||elem->subdomain_id() == 23||elem->subdomain_id() == 24)
			{
				NHS_RK2_step(Ca_i, Ca_b, Q1, Q2, Q3, z, lambda, dlambda_dt, time, dt);

				act_data.set(act_dof_indices[CA_I_IDX][0], Ca_i);
				act_data.set(act_dof_indices[CA_B_IDX][0], Ca_b);
				act_data.set(act_dof_indices[  Q1_IDX][0],   Q1);
				act_data.set(act_dof_indices[  Q2_IDX][0],   Q2);
				act_data.set(act_dof_indices[  Q3_IDX][0],   Q3);
				act_data.set(act_dof_indices[   Z_IDX][0],    z);

				const double zz = max(min(lambda,1.2),0.8);
		#ifdef DEBUG_CHECK_ASSERTIONS
				assert(n   == 3.0);
				assert(n_r == 3.0);
		#endif
				const double z_p_n_r = z_p*z_p*z_p;
				const double K_Z_n_r = K_Z*K_Z*K_Z;

				const double Ca_50 = Ca_50_ref*(1.0+beta_1*(zz-1.0));
				const double Ca_TRPN_50 = Ca_TRPN_max*Ca_50/(Ca_50+(k_refoff/k_on)*(1.0-(1.0+beta_0*(zz-1.0))*0.5/gamma_trpn));

				const double Ca_TRPN_50_Ca_TRPN_max_n = (Ca_TRPN_50*Ca_TRPN_50*Ca_TRPN_50)/(Ca_TRPN_max*Ca_TRPN_max*Ca_TRPN_max);

				const double K1 = alpha_r2*(z_p_n_r/z_p)*n_r*K_Z_n_r/((z_p_n_r+K_Z_n_r)*(z_p_n_r+K_Z_n_r));
				const double K2 = alpha_r2*(z_p_n_r/(z_p_n_r+K_Z_n_r))*(1.0-n_r*K_Z_n_r/(z_p_n_r+K_Z_n_r));
				const double z_max = (alpha_0/Ca_TRPN_50_Ca_TRPN_max_n-K2)/(alpha_r1+K1+alpha_0/Ca_TRPN_50_Ca_TRPN_max_n);

				const double T_0_max = T_ref*(1.0+beta_0*(zz-1.0));
				const double T_0 = T_0_max*z/z_max;

				const double Q_sum = Q1+Q2+Q3;
				double T = (Q_sum < 0.0 ? T_0*(a*Q_sum+1.0)/(1.0-Q_sum) : T_0*(1.0+(2.0+a)*Q_sum)/(1.0+Q_sum));
				////////////////////////////////////////////////////////Here We use simplified linear + qudratic relation 
				double T_quadratic=T_ref*(1.0+beta_0*(zz-1.0));
				double t_quad=-1.0;
				if (time >=1.0)
				{
					double time_in_period = fmod(time-1.0, 0.8);
					if (elem->subdomain_id()==23 || elem->subdomain_id()==24)
					{
						if (time_in_period > 0.7 && time_in_period <= 0.8)
						{
							t_quad=time_in_period-0.7;
							T_quadratic=t_quad<0.06? T_quadratic*t_quad/0.06:T_quadratic*(0.1-t_quad)/0.04;
						}
						else
						{
							T_quadratic=0.0;
						}
						T = T_quadratic;
					}
					else if (elem->subdomain_id()==21 || elem->subdomain_id()==22)
					{
						if (time_in_period > 0.0 && time_in_period <= 0.35)
						{
							// do nothing 
							if (time_in_period > 0.34 && time_in_period <= 0.35) //here give 0.05s relaxation time to help avoid oscillation
							{
								t_quad = time_in_period-0.34;
								T = T*(0.01-t_quad)/0.01;
							}
						}
						else
						{
							T=0.0;
						}
					}
					
				}
				else
				{
					T = 0.0;
				}
				
				////////////////////////////////////////////////////////
				T = max(T, 0.0); //T should always be greater than zero
				T_data.set(T_dof_indices[0], T);
				T_max = max(T,T_max);
			}
			else
			{
				act_data.set(act_dof_indices[CA_I_IDX][0], Ca_i);
				act_data.set(act_dof_indices[CA_B_IDX][0], Ca_b);
				act_data.set(act_dof_indices[  Q1_IDX][0],   Q1);
				act_data.set(act_dof_indices[  Q2_IDX][0],   Q2);
				act_data.set(act_dof_indices[  Q3_IDX][0],   Q3);
				act_data.set(act_dof_indices[   Z_IDX][0],    z);
				const double T = 0.0;
				
				T_data.set(T_dof_indices[0], T);
				T_max = max(T,T_max);
			}
			
			
			Ca_i_max = max(Ca_i,Ca_i_max);
			Ca_i_min = min(Ca_i,Ca_i_min);
			Ca_b_max = max(Ca_b,Ca_b_max);
			Ca_b_min = min(Ca_b,Ca_b_min);
			
		}
		T_data.close();
		act_data.close();
		
		Ca_i_max = SAMRAI_MPI::maxReduction(Ca_i_max);
		Ca_i_min = SAMRAI_MPI::minReduction(Ca_i_min);
		Ca_b_max = SAMRAI_MPI::maxReduction(Ca_b_max);
		Ca_b_min = SAMRAI_MPI::minReduction(Ca_b_min);
		T_max = SAMRAI_MPI::maxReduction(T_max);
		
		
		pout<< "Ca_i_max = "<< Ca_i_max<< "\n" << "Ca_i_min="<< Ca_i_min<<"\n";
		pout<< "Ca_b_max = "<< Ca_b_max<< "\n" << "Ca_b_min="<< Ca_b_min<<"\n";
		pout<<"T_max= "<<T_max<<"\n";

	
	}
	else
	{
		// do nothing	for now
	}
    
	
    return;
}// update_active_tension_model_state_variables

void
ActiveContraction::NHS_RK2_step(
    double& Ca_i,
    double& Ca_b,
    double& Q1,
    double& Q2,
    double& Q3,
    double& z,
    const double lambda,
    const double dlambda_dt,
    const double time,
    const double dt)
{
    double Ca_i_new = Ca_i;
    double Ca_b_new = Ca_b;
    double   Q1_new = Q1;
    double   Q2_new = Q2;
    double   Q3_new = Q3;
    double    z_new = z;
    NHS_euler_step(Ca_i_new, Ca_b_new, Q1_new, Q2_new, Q3_new, z_new, lambda, dlambda_dt, time   , dt);
    NHS_euler_step(Ca_i_new, Ca_b_new, Q1_new, Q2_new, Q3_new, z_new, lambda, dlambda_dt, time+dt, dt);
    Ca_i = 0.5*(Ca_i+Ca_i_new);
    Ca_b = 0.5*(Ca_b+Ca_b_new);
    Q1   = 0.5*(Q1  +  Q1_new);
    Q2   = 0.5*(Q2  +  Q2_new);
    Q3   = 0.5*(Q3  +  Q3_new);
    z    = 0.5*(z   +   z_new);
    return;
}// NHS_RK2_step

void
ActiveContraction::NHS_euler_step(
    double& Ca_i,
    double& Ca_b,
    double& Q1,
    double& Q2,
    double& Q3,
    double& z,
    const double lambda,
    const double dlambda_dt,
    const double time,
    const double dt)
{
    // Simple Ca transient.
    static const double Ca_max = 1.0;   // uM
    static const double Ca_o = 0.01;    // uM
    static const double tau_Ca = 0.06;  // sec
	//static const double tau_Ca = 0.15;  // this is to prolong left ventricle contractio period
	
    double time_in_local_period = fmod(time-1.0, 0.8);
    double t_shift = -1.0;
	if (time > 1.0 )
	{
		t_shift = time_in_local_period;
	}
    //const double t_shift = time-1.0*BoundaryConditions::t_load;
    const double dCa_i_dt = t_shift > 0.0 ? ((Ca_max-Ca_o)/tau_Ca)*exp(1.0-t_shift/tau_Ca)*(1.0-t_shift/tau_Ca) : 0.0;//This Ca profile is no longer used. 2018/02/28
	/*
	double efficient = 0.0;
	if (t_shift > 0.0)
	{
		if (t_shift <= 0.04)
		{
			efficient = 24.0;
		}
		else if (t_shift <= 0.26)
		{
			efficient = 10.0;
		}
		else if (t_shift <= 0.4)
		{
			efficient = -10.0;
		}
		else
		{
			efficient = 0.0;
		}
	}
	const double dCa_i_dt = t_shift > 0.0 ? efficient : 0.0;//This Ca profile is no longer used. 2018/02/28
	*/
	//////////////////////////////////////////////////////////////////////////////////Here we construct a simplified Ca_ profile described with a polynomial with degree =2, -275x^2+33x+0.01 with t_peak=0.06s
	// const double dCa_i_dt = t_shift > 0.0 ? -550.0*t_shift+33.0 : 0.0;

    // The model is only valid for 0.8 <= lambda <= 1.15.
    const double zz = max(min(lambda,1.2),0.8);

    // Tropomyosin kinetics.
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(n   == 3.0);
    assert(n_r == 3.0);
#endif
    const double z_p_n_r = z_p*z_p*z_p;
    const double K_Z_n_r = K_Z*K_Z*K_Z;

    const double z_n_r = z*z*z;

    const double Ca_50 = Ca_50_ref*(1.0+beta_1*(zz-1.0));
    const double Ca_TRPN_50 = Ca_TRPN_max*Ca_50/(Ca_50+(k_refoff/k_on)*(1.0-(1.0+beta_0*(zz-1.0))*0.5/gamma_trpn));

    const double Ca_b_Ca_TRPN_50_n = (Ca_b*Ca_b*Ca_b)/(Ca_TRPN_50*Ca_TRPN_50*Ca_TRPN_50);
    const double Ca_TRPN_50_Ca_TRPN_max_n = (Ca_TRPN_50*Ca_TRPN_50*Ca_TRPN_50)/(Ca_TRPN_max*Ca_TRPN_max*Ca_TRPN_max);

    const double K1 = alpha_r2*(z_p_n_r/z_p)*n_r*K_Z_n_r/((z_p_n_r+K_Z_n_r)*(z_p_n_r+K_Z_n_r));
    const double K2 = alpha_r2*(z_p_n_r/(z_p_n_r+K_Z_n_r))*(1.0-n_r*K_Z_n_r/(z_p_n_r+K_Z_n_r));
    const double z_max = (alpha_0/Ca_TRPN_50_Ca_TRPN_max_n-K2)/(alpha_r1+K1+alpha_0/Ca_TRPN_50_Ca_TRPN_max_n);

    const double dz_dt = alpha_0*Ca_b_Ca_TRPN_50_n*(1.0-z)-alpha_r1*z-alpha_r2*z_n_r/(z_n_r+K_Z_n_r);

    // Tension development and crossbridge dynamics.
    const double T_0_max = T_ref*(1.0+beta_0*(zz-1.0));
    const double T_0 = T_0_max*z/z_max;

    const double Q_sum = Q1+Q2+Q3;
    const double T = (Q_sum < 0.0 ? T_0*(a*Q_sum+1.0)/(1.0-Q_sum) : T_0*(1.0+(2.0+a)*Q_sum)/(1.0+Q_sum));

    const double dQ1_dt = A1*dlambda_dt-alpha_1*Q1;
    const double dQ2_dt = A2*dlambda_dt-alpha_2*Q2;
    const double dQ3_dt = A3*dlambda_dt-alpha_3*Q3;

    // Troponin C-Calcium binding.
    double k_off = k_refoff*(1.0-T/(gamma_trpn*T_ref));
    k_off = max(k_off, 0.0);
    const double dCa_b_dt = k_on*Ca_i*(Ca_TRPN_max-Ca_b)-k_off*Ca_b;

    // Update time-dependent variables.
   //forcing Ca_i = 0 when reaching time BoundaryConditions::t_end_systole, then the active force T will be zero
   if (time >= 1.0 && time_in_local_period >=0.0 && time_in_local_period <=0.35)
   {
      
	   Ca_i += dt*dCa_i_dt;  Ca_i = max(0.0,Ca_i);
        Ca_b += dt*dCa_b_dt;  Ca_b = max(0.0,Ca_b);
        z    += dt*dz_dt;
		if (time_in_local_period > 0.34 && time_in_local_period <= 0.35) //here give 0.05s relaxation time to help avoid oscillation
		{
			Ca_i = Ca_i*(0.35-time_in_local_period)/0.01;
			Ca_b = Ca_b*(0.35-time_in_local_period)/0.01;
			z = z*(0.35-time_in_local_period)/0.01;
		}
   }
   else
   {
		 Ca_i = 0.0;
       Ca_b = 0.0;
       z    = 0.0;
        
   }
   Q1   += dt*dQ1_dt;
   Q2   += dt*dQ2_dt;
   Q3   += dt*dQ3_dt;

    return;
}// NHS_euler_step
