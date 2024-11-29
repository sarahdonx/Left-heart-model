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
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>
// IBAMR INCLUDES
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
// MechanicsModel is a static class that provides data and functions required to
// implement the LV mechanics model.
class MechanicsModel
{
public:
	static  std::vector<IBTK::SystemData> f0_system_data;
	static  std::vector<IBTK::SystemData> s0_system_data;
	static 	std::vector<NumericVector<double>*> active_tension_data;
    static int f0_mv_system_num, s0_mv_system_num;
	static int tubes_f0_mv_system_num, tubes_s0_mv_system_num;
    static bool normalize_stress;
    static double C1_s, beta_s;
	static double MV_start_close_time, MV_start_cap_time, MV_start_annulus_time,MV_end_close_time;
    static double I1_dev_max, I1_dev_min, I1_dil_max, I1_dil_min;
    static double J_dev_max, J_dev_min, J_dil_max, J_dil_min;
	static bool enable_active_tension;
    static double T_scale;

    static void
    get_PK1_dev_stress_function_systems(
        std::vector<SystemData>& system_data);
	static void
    get_PK1_dil_stress_function_systems(
        std::vector<SystemData>& system_data);
	static void
    update_active_tension_variables(
        EquationSystems* equation_systems,
        double time,
        double dt);
    static void
    PK1_dev_stress_function(
        TensorValue<double>& PP,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const std::vector<const std::vector<double>*>& var_data,
		const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
        double data_time,
        void* ctx);

    static void
    PK1_dil_stress_function(
        TensorValue<double>& PP,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const std::vector<const std::vector<double>*>& var_data,
		const std::vector<const std::vector<VectorValue<double> >*>& grad_var_data,
        double data_time,
        void* ctx);

private:
    MechanicsModel();
    MechanicsModel(MechanicsModel&);
    ~MechanicsModel();
    MechanicsModel& operator=(MechanicsModel&);
};
