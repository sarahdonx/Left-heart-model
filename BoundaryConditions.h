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
#include <ibtk/libmesh_utilities.h>
// LIBMESH INCLUDES
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/quadrature.h>
#include <libmesh/point.h>

#include <string>
// BoundaryConditions is a static class that provides data and functions
// required to implement loading conditions and boundary constraints for the LV
// model.
class BoundaryConditions
{
public:
	static double P_load, t_load, kappa_s, kappa_b, P_loading;
    static BoundaryInfo* boundary_info;
	static double LV_volume; 
	static double LAA_volume; 
	static std::vector<int>  LV_endo_points_list;
    static std::vector< std::vector<double>  > LV_endo_points;
    static int LV_NoOfEndoNode;
	static std::vector<int>  RV_endo_points_list;
    static std::vector< std::vector<double>  > RV_endo_points;
    static int RV_NoOfEndoNode;
	static std::vector<int>  LA_endo_points_list;
    static std::vector< std::vector<double>  > LA_endo_points;
    static int LA_NoOfEndoNode;
	static std::vector<int>  RA_endo_points_list;
    static std::vector< std::vector<double>  > RA_endo_points;
    static int RA_NoOfEndoNode;
	static std::vector< std::vector<double> > mapping_peri;
    static double
    loading_pressure(
        double time);

    static void
    Peri_traction_function(
		double& P,
		const VectorValue<double>& n,
		const VectorValue<double>& N,
		const TensorValue<double>& FF,
		const libMesh::Point& x,
		const libMesh::Point& X,
		Elem* const elem,
		unsigned short int side,
		const vector<const vector<double>*>& var_data,
		const vector<const vector<VectorValue<double> >*>& grad_var_data,
		double time,
		void* /*ctx*/);
	static void	
	loading_force_function(
		double& P,
		const VectorValue<double>& /*n*/,
		const VectorValue<double>& /*N*/,
		const TensorValue<double>& /*FF*/,
		const libMesh::Point& /*x*/,
		const libMesh::Point& /*X*/,
		Elem* const elem,
		unsigned short int side,
		const vector<const vector<double>*>& /*var_data*/,
		const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
		double time,
		void* /*ctx*/);
		
    static void
    tether_force_function(
        VectorValue<double>& F,
		const VectorValue<double>& n,
        const VectorValue<double>& N,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const unsigned short int side,
        const vector<const vector<double>*>& var_data,
		const vector<const vector<VectorValue<double> >*>& grad_var_data,
        double time,
        void* ctx);
		
	
		
   static void
    tether_vol_force_function(
        VectorValue<double>& F,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        //const unsigned short int side,
        const vector<const vector<double>*>& var_data,
		const vector<const vector<VectorValue<double> >*>& grad_var_data,
        double time,
        void* ctx);
	static void 
    readingPeriMapping(std::vector< std::vector<double> > peri_mapping_read,
							std::vector< std::vector<double> >& mapping_peri);
	static void 
    readingPointsGeneral(MeshBase& mesh, 
                                    std::vector<int>& points_list, 
                                    std::vector< std::vector<double> >& points,
                                    int& NoOfPoints,
                                    std::string file_name);

	
	static void 
      updatePointsPositionGeneralLV(EquationSystems* equation_systems, 
                                                    std::vector<int> & points_list, 
                                     std::vector< std::vector<double> >& points,
                                                    int NoOfPoints);			
	static void 
      updatePointsPositionGeneralRV(EquationSystems* equation_systems, 
                                                    std::vector<int> & points_list, 
                                     std::vector< std::vector<double> >& points,
                                                    int NoOfPoints);
	static void 
      updatePointsPositionGeneralLA(EquationSystems* equation_systems, 
                                                    std::vector<int> & points_list, 
                                     std::vector< std::vector<double> >& points,
                                                    int NoOfPoints);
	static void 
      updatePointsPositionGeneralRA(EquationSystems* equation_systems, 
                                                    std::vector<int> & points_list, 
                                     std::vector< std::vector<double> >& points,
                                                    int NoOfPoints);													

private:
    BoundaryConditions();
    BoundaryConditions(BoundaryConditions&);
    ~BoundaryConditions();
    BoundaryConditions& operator=(BoundaryConditions&);
};
