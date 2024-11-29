// Copyright (c) 2011-2013, Boyce Griffith
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
// IBAMR INCLUDES
#include <ibtk/libmesh_utilities.h>

// APPLICATION INCLUDES
#include <ActiveContraction.h>
#include <MechanicsModel.h>

// LIBMESH INCLUDES
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
// STATIC VARIABLES
//temp
std::vector<IBTK::SystemData> MechanicsModel::f0_system_data;
std::vector<IBTK::SystemData> MechanicsModel::s0_system_data;
std::vector<NumericVector<double>*> MechanicsModel::active_tension_data;
int MechanicsModel::f0_mv_system_num, MechanicsModel::s0_mv_system_num;
int MechanicsModel::tubes_f0_mv_system_num, MechanicsModel::tubes_s0_mv_system_num;
double MechanicsModel::C1_s, MechanicsModel::beta_s;
double MechanicsModel::MV_start_close_time, MechanicsModel::MV_start_cap_time, MechanicsModel::MV_start_annulus_time, MechanicsModel::MV_end_close_time;
bool  MechanicsModel::normalize_stress;
double MechanicsModel::I1_dev_max;
double MechanicsModel::I1_dev_min;
double MechanicsModel::J_dev_max;
double MechanicsModel::J_dev_min;
double MechanicsModel::I1_dil_max;
double MechanicsModel::I1_dil_min;
double MechanicsModel::J_dil_max;
double MechanicsModel::J_dil_min;
bool MechanicsModel::enable_active_tension;
double MechanicsModel::T_scale;
namespace  // private namespace
{
	//For LA 
	static const double a_LA = 0.5*5580.0;// dyne/cm^2
	static const double b_LA = 4.033;// nondimensional
    static const double a1_LA = 0.5*86520.0;// dyne/cm^2
	static const double b1_LA = 10.536;// nondimensional
	//For LV 
	static const double a_LV = 2790.0;// dyne/cm^2
	static const double b_LV = 4.033;// nondimensional
    static const double a1_LV = 43260.0;// dyne/cm^2
	static const double b1_LV = 10.536;// nondimensional
	static const double a2_LV = 8940.0;// dyne/cm^2
	static const double b2_LV = 4.284;// nondimensional
	//For MV model 
	static const double C10_ant = 1.245e3;
	static const double C01_ant = 13.6655;
    static const double k1_ant = 1.10069e5;// dyne/cm^2
	static const double k2_ant = 84.8478;// nondimensional
	static const double kappa_ant = 0.0800;
	
    static const double C10_post = 5.02e2;
	static const double C01_post = 15.0036;
    static const double k1_post = 3.0207e4;// dyne/cm^2
	static const double k2_post = 144.4848;// nondimensional
	static const double kappa_post = 0.0534;
	
	// AV and PAV
	/*
	static const double C11 = 4.74826;
	static const double C12 = 273732.0;// dyne/cm^2
    static const double C41 = 80.4291;
	static const double C42 = 119752.0;// // dyne/cm^2
	*/
	static const double C11 = 3.76;
	static const double C12 = 0.25*3.2019e6;// dyne/cm^2
    static const double C41 = 6.55;
	static const double C42 = 8.9384e6;// // dyne/cm^2
}

// CLASS IMPLEMENTATION

void
MechanicsModel::get_PK1_dev_stress_function_systems(
    std::vector<SystemData>& system_data)
{
	//temp
    system_data.resize(3);
	system_data[0].system_name = "fiber direction";
	system_data[1].system_name = "sheet direction";
	system_data[2].system_name = "active tension";
	system_data[0].vars.resize(NDIM);
	system_data[1].vars.resize(NDIM);
	system_data[2].vars.resize(1);
	for (unsigned int d = 0; d < NDIM; ++d) 
	{
		system_data[0].vars[d] = d;	
		system_data[1].vars[d] = d;	
	}
	system_data[2].vars[0] = 0;	
	f0_system_data.resize(1);
	s0_system_data.resize(1);
	f0_system_data[0]=system_data[0];
	s0_system_data[0]=system_data[1];
    return;
}// get_PK1_dev_stress_function_systems
void
MechanicsModel::get_PK1_dil_stress_function_systems(
    std::vector<SystemData>& system_data)
{
	//temp
    system_data.resize(3);
	system_data[0].system_name = "fiber direction";
	system_data[1].system_name = "sheet direction";
	system_data[2].system_name = "active tension";
	system_data[0].vars.resize(NDIM);
	system_data[1].vars.resize(NDIM);
	system_data[2].vars.resize(1);
	for (unsigned int d = 0; d < NDIM; ++d) 
	{
		system_data[0].vars[d] = d;	
		system_data[1].vars[d] = d;	
	}
	system_data[2].vars[0] = 0;	
    return;
}// get_PK1_dil_stress_function_systems

void
MechanicsModel::update_active_tension_variables(
    EquationSystems* equation_systems,
    const double time,
    const double dt)
{
	
	System& T_system = equation_systems->get_system<System>(ActiveContraction::T_system_num);
	T_system.update();
	active_tension_data.resize(1);
	active_tension_data[0] = T_system.current_local_solution.get();
	return;
}		
void
MechanicsModel::PK1_dev_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /* X */,  // current   location
    const libMesh::Point& /* s */,  // reference location
    Elem* const elem,
    const std::vector<const std::vector<double>*>& var_data,
    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
    double data_time,
    void* /* ctx */)
{
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF,NDIM);

    // Isotropic contribution.
    const double I1 = CC.tr();
    const double J = FF.det();
	const double I1_bar = I1/(std::pow(J,2.0/3));

    I1_dev_max = max(I1_dev_max, I1);
    I1_dev_min = min(I1_dev_min, I1);

    J_dev_max = max(J_dev_max, J);
    J_dev_min = min(J_dev_min, J);

    // Fiber contribution.
	//temp
    const std::vector<double>& f0_vec = *var_data[0];
    VectorValue<double> f0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        //temp
        //f0(d) = (*f0_vec)(elem->dof_number(f0_mv_system_num,d,0));
		f0(d) = f0_vec[d];	
    }
	
    const double I4f = f0*(CC*f0);
	const double Ef_ant = (I1-3.0)*kappa_ant+(I4f-1.0)*(1.0-3.0*kappa_ant);
	const double Ef_post = (I1-3.0)*kappa_post+(I4f-1.0)*(1.0-3.0*kappa_post);
	
	// Sheet contribution.
	//temp
    const std::vector<double>& s0_vec = *var_data[1];
    VectorValue<double> s0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        //temp
        //s0(d) = (*s0_vec)(elem->dof_number(s0_mv_system_num,d,0));
		s0(d) = s0_vec[d];	
    }
    const double I4s = s0*(CC*s0);
	const double Es_ant = (I1-3.0)*kappa_ant+(I4s-1.0)*(1.0-3.0*kappa_ant);
	const double Es_post = (I1-3.0)*kappa_post+(I4s-1.0)*(1.0-3.0*kappa_post);
   ///different subdomain will use different material model
   // 30: anterior
   // 40: posterior
   // 3: TV
   // 5: Housing
   // 6: AV
   // 7: PAV
   //23，24: BiAtrium
   //21，22: BiVentricle
   // 5001,5002,5003: movable annulus and housing portion
   // 31,32: Biventricle transition portion
   //100 All tube structure
	if (elem->subdomain_id() >= 3 )//Everything
	{

        if (elem->subdomain_id() == 30 ||elem->subdomain_id() == 5)//AML
        {
			PP = 2.0*C10_ant*C01_ant*exp(C01_ant*(I1_bar-3.0))*FF;
			PP = PP/(std::pow(J,2.0/3));			
			
            if (Ef_ant > 0.0)
            {
                PP += 2.0*k1_ant*kappa_ant*Ef_ant*exp(k2_ant*Ef_ant*Ef_ant)*FF+2.0*k1_ant*(1.0-3.0*kappa_ant)*Ef_ant*exp(k2_ant*Ef_ant*Ef_ant)*FF*outer_product(f0,f0);
            } 
			
			if (Es_ant > 0.0)
            {
                PP += 2.0*k1_ant*kappa_ant*Es_ant*exp(k2_ant*Es_ant*Es_ant)*FF+2.0*k1_ant*(1.0-3.0*kappa_ant)*Es_ant*exp(k2_ant*Es_ant*Es_ant)*FF*outer_product(s0,s0);
            }
        }

        else if (elem->subdomain_id() == 40||elem->subdomain_id() == 5003)//PML
        {
			PP = 2.0*C10_post*C01_post*exp(C01_post*(I1_bar-3.0))*FF;
			PP = PP/(std::pow(J,2.0/3));			
			
            if (Ef_post > 0.0)
            {
                PP += 2.0*k1_post*kappa_post*Ef_post*exp(k2_post*Ef_post*Ef_post)*FF+2.0*k1_post*(1.0-3.0*kappa_post)*Ef_post*exp(k2_post*Ef_post*Ef_post)*FF*outer_product(f0,f0);
            } 
			
			if (Es_post > 0.0)
            {
                PP += 2.0*k1_post*kappa_post*Es_post*exp(k2_post*Es_post*Es_post)*FF+2.0*k1_post*(1.0-3.0*kappa_post)*Es_post*exp(k2_post*Es_post*Es_post)*FF*outer_product(s0,s0);
            }

        }
		else if (elem->subdomain_id() == 23 || elem->subdomain_id() == 24) //biAtrium 5001 5002 are moveble portion of TV annulus, MVannulus+housing parts.
        {
           PP = a_LA*exp(b_LA*(I1_bar-3.0))*FF; 
		   PP = PP/(std::pow(J,2.0/3));
		   if (I4f > 1.0)
		   {
			   PP +=2*a1_LA*exp(b1_LA*(I4f-1)*(I4f-1))*(I4f-1)*FF*outer_product(f0,f0);
		   }

       }
	   else if (elem->subdomain_id() == 21 || elem->subdomain_id() == 22|| elem->subdomain_id() == 211)//biVentricle
        {
           PP = a_LV*exp(b_LV*(I1_bar-3.0))*FF; 
		   PP = PP/(std::pow(J,2.0/3));
		   if (I4f > 1.0)
		   {
			   PP +=2*a1_LV*exp(b1_LV*(I4f-1)*(I4f-1))*(I4f-1)*FF*outer_product(f0,f0);
		   }
		   if (I4s > 1.0)
		   {
			   PP +=2*a2_LV*exp(b2_LV*(I4s-1)*(I4s-1))*(I4s-1)*FF*outer_product(s0,s0);
		   }

       }
	   else if (elem->subdomain_id() == 3)// TV try with softer AV material
	   {
		   ////commented section is little soft for TV valve. so switch to same material as MV leaflets
		   
		  PP = 2.0*C10_ant*C01_ant*exp(C01_ant*(I1_bar-3.0))*FF; 
		  PP = PP/(std::pow(J,2.0/3));
			
            if (Ef_ant > 0.0)
            {
                PP += 2.0*k1_ant*kappa_ant*Ef_ant*exp(k2_ant*Ef_ant*Ef_ant)*FF+2.0*k1_ant*(1.0-3.0*kappa_ant)*Ef_ant*exp(k2_ant*Ef_ant*Ef_ant)*FF*outer_product(f0,f0);
            }
			PP=PP*1.0;						
	   }	   
	   else	if (elem->subdomain_id() == 100||elem->subdomain_id() == 200)// All vessel tubes
	   {
		   PP = C12*exp(C11*(I1_bar-3.0))*FF; 
		   PP = PP/(std::pow(J,2.0/3));
	   }
	   else	if (elem->subdomain_id() == 1000)// All vessel tubes
	   {
		   PP = (1.0/0.25)*100*C12*exp(C11*(I1_bar-3.0))*FF; 
		   PP = PP/(std::pow(J,2.0/3));
	   }
	   else	if (elem->subdomain_id() == 4099)// All vessel tubes
	   {
		   PP = a_LV*exp(b_LV*(I1_bar-3.0))*FF; 
		   PP = PP/(std::pow(J,2.0/3));
	   }
	   
	}
	else 
	{
		//should not reach here
		PP.zero();
	}
    
    return;
}// PK1_dev_stress_function

void
MechanicsModel::PK1_dil_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /* X */,  // current   location
    const libMesh::Point& /* s */,  // reference location
    Elem* const elem,
    const std::vector<const std::vector<double>*>& var_data,
    const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
    double  data_time ,
    void* /* ctx */)
{
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF,NDIM);
    const double I1 = CC.tr();
    const double J = FF.det();
	const double I1_bar = I1/(std::pow(J,2.0/3));

    I1_dil_max = max(I1_dil_max, I1);
    I1_dil_min = min(I1_dil_min, I1);

    J_dil_max = max(J_dil_max, J);
    J_dil_min = min(J_dil_min, J);
    
    PP.zero();
	double beta_s_NoLocking = 0.0;
	double v_stab = 0.495;
	double G = 0.0;
    //incompressibility
	if (elem->subdomain_id() >= 3 )
	{
       
	    if (elem->subdomain_id() == 30||elem->subdomain_id() == 5)
        {
			beta_s_NoLocking = 1.0E+06;
			PP += (beta_s_NoLocking*log(CC.det()))*FF_inv_trans;
			PP = PP - (I1_bar/3.0)*2.0*C10_ant*C01_ant*exp(C01_ant*(I1_bar-3.0))*FF_inv_trans;
        }
		else if (elem->subdomain_id() == 40||elem->subdomain_id() == 5003)
		{
			beta_s_NoLocking = 1.0E+06;
			PP += (beta_s_NoLocking*log(CC.det()))*FF_inv_trans;
			PP = PP - (I1_bar/3.0)*2.0*C10_post*C01_post*exp(C01_post*(I1_bar-3.0))*FF_inv_trans;
		}
		else if (elem->subdomain_id() == 23 || elem->subdomain_id() == 24)
		{
			beta_s_NoLocking = 1.0E+06;
			PP += (beta_s_NoLocking*log(CC.det()))*FF_inv_trans;
			PP = PP -(I1_bar/3.0)*a_LA*exp(b_LA*(I1_bar-3.0))*FF_inv_trans;
		}
		else if (elem->subdomain_id() == 21 || elem->subdomain_id() == 22|| elem->subdomain_id() == 211)
		{
			beta_s_NoLocking = 3.5E+06;
			PP += (beta_s_NoLocking*log(CC.det()))*FF_inv_trans;
			PP = PP -(I1_bar/3.0)*a_LV*exp(b_LV*(I1_bar-3.0))*FF_inv_trans;
		}
		else if (elem->subdomain_id() == 3)
		{
			beta_s_NoLocking = 1.0E+06;
			PP += (beta_s_NoLocking*log(CC.det()))*FF_inv_trans;
			PP = PP - (I1_bar/3.0)*1.0*2.0*C10_ant*C01_ant*exp(C01_ant*(I1_bar-3.0))*FF_inv_trans;
			
			
		}
		else if (elem->subdomain_id() == 100||elem->subdomain_id() == 200)
        {
			beta_s_NoLocking = 3.5E+06;
			PP += (beta_s_NoLocking*log(CC.det()))*FF_inv_trans;
			PP = PP -(I1_bar/3.0)*C12*exp(C11*(I1_bar-3.0))*FF_inv_trans;
        }
		else if (elem->subdomain_id() == 1000 )
        {
			beta_s_NoLocking = 3.5E+06;
			PP += (beta_s_NoLocking*log(CC.det()))*FF_inv_trans;
			PP = PP -(I1_bar/3.0)*(1.0/0.25)*100*C12*exp(C11*(I1_bar-3.0))*FF_inv_trans;
        }
		else if (elem->subdomain_id() == 4099)
        {
			beta_s_NoLocking = 1.0E+06;
			PP += (beta_s_NoLocking*log(CC.det()))*FF_inv_trans;
			PP = PP -(I1_bar/3.0)*a_LV*exp(b_LV*(I1_bar-3.0))*FF_inv_trans;
        }
	
	}
	else
	{
		//should not reach here
		PP.zero(); // do nothing
	} 
	
	if (elem->subdomain_id() == 21||elem->subdomain_id() == 22||elem->subdomain_id() == 23||elem->subdomain_id() == 24) // BiAtrium and BiVentricle
	{
		
		if (enable_active_tension)
		{
			// Fiber contribution.
			const std::vector<double>& f0_vec = *var_data[0];
			VectorValue<double> f0;
			for (unsigned int d = 0; d < NDIM; ++d)
			{
				f0(d) = f0_vec[d];	
			}
			// Active tension contribution.		
			NumericVector<double>* T_vec = active_tension_data[0];
			double T_temp = 0.0;
			if (elem->subdomain_id() == 21)
			{
				T_temp = 1.1*1.6*1.5*4.0;
			}
			else if (elem->subdomain_id() == 22)
			{
				
				T_temp = 1.1*1.8*1.5*4.0; //1.5->1.0			
			}
			else if (elem->subdomain_id() == 23)
			{
				T_temp = 1.0;//0.75->0.5
			}
			else if (elem->subdomain_id() == 24)
			{
				T_temp = 2.0;
			}
			const double T = T_temp*(*T_vec)(elem->dof_number(ActiveContraction::T_system_num,0,0))*1.0e4;  // NOTE: T is stored in kPa; must be converted to dyne/cm^2			
			if (T > 0.0)
			{
				PP += FF.det()*T*FF*outer_product(f0,f0);	
			}
			
		}
		
	}
	
	 
	
	
    return;
}// PK0_dil_stress_function
