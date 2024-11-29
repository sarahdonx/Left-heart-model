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

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petsc.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFECentroidPostProcessor.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <libmesh/centroid_partitioner.h>
// Headers for application components.
#include <ActiveContraction.h>
#include <BoundaryConditions.h>
#include <MechanicsModel.h>
#include <ModelInitialization.h>
#include <FeedbackForcer.h>
#include <VelocityBcCoefs.h>
#include <Windkessel.h>
//for linux system call
#include <stdlib.h>
#include <stdio.h>
//For IB
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/IBStrategySet.h>
#include <petscsys.h>
#include "ibamr/IBTargetPointForceSpec.h"
//
#include <libmesh/equation_systems.h>
#include <boost/multi_array.hpp>
//
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include "petscvec.h"
// Forward declaration of the main driver routine that actually runs the
// simulation.
//chordae spring force function 
inline double 
chordae_spring_force(
	const double r,
	const double* const params,
	const int lag_mastr_idx,
	const int lag_slave_idx 
)
{
	//Compute the force applied to the "master" node
	const double stf = params[0];
	const double rst = params[1];
	static const double c0_basal = 540e4;// dyne/cm^2
	static const double af_basal = 723.10e4;// dyne/cm^2
	static const double bf_basal = 22.09;// nondimensional

	static const double c0_marginal = 540e4;// dyne/cm^2
	static const double af_marginal = 100.24e4;// dyne/cm^2
	static const double bf_marginal = 33.83;// nondimensional
	double C10=0.0;
	double K1=0.0;
	double K2=0.0;
	double sigma1=0.0;
	
	double C11 = 4.74826;
	double C12 = 273732.0;// dyne/cm^2
    double C41 = 80.4291;
	double C42 = 119752.0;// // dyne/cm^2
	
	if (r > rst)
	{
		if (lag_mastr_idx >= 1449 && lag_mastr_idx <= 1562) // this is only for MV marginal chordae
		{
			C10=c0_marginal;
			K1=af_marginal;
			K2=bf_marginal;
			sigma1=2*C10*((r/rst)*(r/rst)-rst/r)+4*(r/rst)*(r/rst)*((r/rst)*(r/rst)-1)*K1*exp(K2*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1));
			sigma1=0.5625*0.0016*(rst/r)*sigma1;
	
		}
		else if (lag_mastr_idx <= 1565) // 
		{
			C10=c0_basal;
			K1=af_basal;
			K2=bf_basal;
			sigma1=2*C10*((r/rst)*(r/rst)-rst/r)+4*(r/rst)*(r/rst)*((r/rst)*(r/rst)-1)*K1*exp(K2*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1));
			sigma1=0.5625*0.0016*(rst/r)*sigma1;

		}
		else if (lag_mastr_idx >= 2871 && lag_mastr_idx <= 3083)//AV Belly is always stiffer than the free  edge region, leaflet 1 
		{
			sigma1=C12*exp(C11*((r/rst)*(r/rst)+2*rst/r-3.0))*((r/rst)*(r/rst)-rst/r)+2.0*C42*((r/rst)*(r/rst)-1)*exp(C41*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1))*(r/rst)*(r/rst);
			sigma1=sigma1*0.01;
		}
		else if (lag_mastr_idx >= 4439 && lag_mastr_idx <= 4648)//AV, leaflet 2 
		{
			sigma1=C12*exp(C11*((r/rst)*(r/rst)+2*rst/r-3.0))*((r/rst)*(r/rst)-rst/r)+2.0*C42*((r/rst)*(r/rst)-1)*exp(C41*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1))*(r/rst)*(r/rst);
			sigma1=sigma1*0.01;
		}
		else if (lag_mastr_idx >= 5980 && lag_mastr_idx <= 6448)//AV, leaflet 3 
		{
			sigma1=C12*exp(C11*((r/rst)*(r/rst)+2*rst/r-3.0))*((r/rst)*(r/rst)-rst/r)+2.0*C42*((r/rst)*(r/rst)-1)*exp(C41*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1))*(r/rst)*(r/rst);
			sigma1=sigma1*0.01;
		}
		else if (lag_mastr_idx <= 6456)//Rest of AV
		{
			sigma1=C12*exp(C11*((r/rst)*(r/rst)+2*rst/r-3.0))*((r/rst)*(r/rst)-rst/r)+2.0*C42*((r/rst)*(r/rst)-1)*exp(C41*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1))*(r/rst)*(r/rst);
			sigma1=sigma1*0.001;
		}
		else	
		{
			//shouldnt reach here;
		}
	}
	else
	{
		if (lag_mastr_idx >= 1449 && lag_mastr_idx <= 1562) // this is only for MV marginal chordae
		{
			sigma1=0.0;
	
		}
		else if (lag_mastr_idx <= 1565) // for now this is also used for TV chordae
		{
			sigma1=0.0;

		}
		else if (lag_mastr_idx >= 2871 && lag_mastr_idx <= 3083)//AV Belly is always stiffer than the free  edge region, leaflet 1 
		{
			sigma1=C12*exp(C11*((r/rst)*(r/rst)+2*rst/r-3.0))*((r/rst)*(r/rst)-rst/r)+2.0*C42*((r/rst)*(r/rst)-1)*exp(C41*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1))*(r/rst)*(r/rst);
			sigma1=sigma1*0.01;
		}
		else if (lag_mastr_idx >= 4439 && lag_mastr_idx <= 4648)//AV, leaflet 2
		{
			sigma1=C12*exp(C11*((r/rst)*(r/rst)+2*rst/r-3.0))*((r/rst)*(r/rst)-rst/r)+2.0*C42*((r/rst)*(r/rst)-1)*exp(C41*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1))*(r/rst)*(r/rst);
			sigma1=sigma1*0.01;
		}
		else if (lag_mastr_idx >= 5980 && lag_mastr_idx <= 6448)//AV, leaflet 3
		{
			sigma1=C12*exp(C11*((r/rst)*(r/rst)+2*rst/r-3.0))*((r/rst)*(r/rst)-rst/r)+2.0*C42*((r/rst)*(r/rst)-1)*exp(C41*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1))*(r/rst)*(r/rst);
			sigma1=sigma1*0.01;
		}
		else if (lag_mastr_idx <= 6456)////shouldnt reach here;
		{
			sigma1=C12*exp(C11*((r/rst)*(r/rst)+2*rst/r-3.0))*((r/rst)*(r/rst)-rst/r)+2.0*C42*((r/rst)*(r/rst)-1)*exp(C41*((r/rst)*(r/rst)-1)*((r/rst)*(r/rst)-1))*(r/rst)*(r/rst);
			sigma1=sigma1*0.001;
		}
		else	
		{
			//shouldnt reach here;
		}
	}	
	return sigma1; //linear from Joyce MV model
	//chordae spring stiffness calculation 
	// YMchord*areaChord/length  in meter if pressure is in pacasal
	// areachord = pi*r^2, if r=0.4mm, then the area will be 0.5mm^2
}
void
update_target_point_positions(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const l_data_manager,
    const double current_time,
    const double dt,
	EquationSystems* equation_systems,
	Mesh& mesh,
	std::vector<std::vector<double>> Pseudo_MVleaflet_mapping,
	std::vector<unsigned int> Pseudo_MVleaflet_mapping_firstComponent,
	std::vector<unsigned int> Pseudo_MVleaflet_mapping_secondComponent,
	std::vector<double> MV_leaflets_MappingElement_currentLocation)
{
	double time_in_period = -1.0;	
	time_in_period = fmod(current_time-1.0, 0.8);
		//////////////////////////////////

		const int level_num = hierarchy->getFinestLevelNumber();

		// Look up the Lagrangian index ranges.
		const std::pair<int,int>& plate2d_left_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, level_num);
	   // const std::pair<int,int>& plate2d_rght_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(1, level_num);

		// Get the Lagrangian mesh.
		Pointer<LMesh> l_mesh = l_data_manager->getLMesh(level_num);
		const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
		const std::vector<LNode*>& ghost_nodes = l_mesh->getGhostNodes();
		std::vector<LNode*> nodes = local_nodes;
		nodes.insert(nodes.end(), ghost_nodes.begin(), ghost_nodes.end());
		
		
		// Loop over all Lagrangian mesh nodes and update the target point
		// positions.
		
		std::vector<unsigned int>::iterator it_Pseudo;
		ofstream ofile_test;
		ofile_test.open("test.out");
		
		for (std::vector<LNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
		{
			const LNode* const node = *it;
			IBTargetPointForceSpec* const force_spec = node->getNodeDataItem<IBTargetPointForceSpec>();
			if (force_spec)
			{
				// Here we update the position of the target point.
				//
				// NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the total number of Lagrangian points)
				//        X_target     is the target position of the target point
				//        X_target(0)  is the x component of the target position
				//        X_target(1)  is the y component of the target position
				//
				// The target position is shifted to the left or right by the
				// increment dt*V
				IBTK::Point& X_target = force_spec->getTargetPointPosition();
				const int lag_idx = node->getLagrangianIndex();
				if (plate2d_left_lag_idxs.first <= lag_idx && lag_idx < plate2d_left_lag_idxs.second)
				{
					double early_sys_s = 0.00;
					double early_sys_d = 0.005;
					double early_sys_duration = early_sys_d-early_sys_s;
					double end_sys_s = 0.38;
					double end_sys_d = 0.48;
					double end_sys_duration = end_sys_d-end_sys_s;
					double XDv =  0.0;  
					double YDv =  0.0;
					double ZDv =  0.0;
					if (lag_idx==1067)
					{
						 XDv =  0.0249    ;  
						 YDv =   0.0140 ;
						 ZDv =   -0.4285;
												
					}
					else if (lag_idx==435)	
					{
						 XDv =   0.1894    ;  
						 YDv =   -0.0251;
						 ZDv =   -0.4614;
						
					}
					else if (lag_idx==64)	
					{
						 XDv =  0.01;  
						 YDv =  -0.01;
						 ZDv =  -0.3352;
						
					}
					else if (lag_idx==557)	
					{
						 XDv =  -0.1855;  
						 YDv =  -0.3734;
						 ZDv =  -0.7096;
						
					}
					else if (lag_idx==455)	
					{
						 XDv = -0.0197      ;  
						 YDv =  -0.6594;
						 ZDv =  -0.3962;
						
					}
					else if (lag_idx==767)	
					{
						 XDv =  0.0607     ; 
						 YDv = -0.5671 ;
						 ZDv =  -0.4381 ;
							
					}
					else if (lag_idx==1300)	
					{
						 XDv =  0.4941   ;  
						 YDv =  -0.4878 ;
						 ZDv =   -0.0377;
						
					}
					else if (lag_idx==179)	
					{
						 XDv =  0.5967     ;  
						 YDv =  -0.3515;
						 ZDv =   -0.0567;
						
					}
					else if (lag_idx==14)	
					{
						 XDv =  0.5087      ;  
						 YDv =  0.0015 ;
						 ZDv =   0.1268;
						
					}
					else if (lag_idx==253)	
					{
						 XDv =   0.3876     ;  
						 YDv =  -0.0851 ;
						 ZDv =   0.0290;
						
					}
					if (time_in_period >= early_sys_s && time_in_period <= early_sys_d)
					{
						X_target(0) += dt*XDv/early_sys_duration;
						X_target(1) += dt*YDv/early_sys_duration;
						X_target(2) += dt*ZDv/early_sys_duration;
						
					}
					else if (time_in_period >= end_sys_s && time_in_period<=end_sys_d)
					{
						X_target(0) -= dt*XDv/end_sys_duration;
						X_target(1) -= dt*YDv/end_sys_duration;
						X_target(2) -= dt*ZDv/end_sys_duration;
						
					}						
				}
			}
		}
			ofile_test.close();


    return;
}// update_target_point_positions
void
update_connection_point_positions(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const l_data_manager,
    const double current_time,
    const double dt,
	std::vector<std::vector<double>> Pseudo_MVleaflet_mapping,
	std::vector<unsigned int> Pseudo_MVleaflet_mapping_firstComponent,
	std::vector<unsigned int> Pseudo_MVleaflet_mapping_secondComponent,
	std::vector<double> MV_leaflets_MappingElement_currentLocation)
	{
		const int coarsest_ln = 0;
		const int finest_ln = hierarchy->getFinestLevelNumber();
		const int level_num = hierarchy->getFinestLevelNumber();
		const std::pair<int,int>& plate2d_left_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, level_num);

		// Get the Lagrangian mesh.
		Pointer<LMesh> l_mesh = l_data_manager->getLMesh(level_num);
		const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
		const std::vector<LNode*>& ghost_nodes = l_mesh->getGhostNodes();
		std::vector<LNode*> nodes = local_nodes;
		nodes.insert(nodes.end(), ghost_nodes.begin(), ghost_nodes.end());
		std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_X_current_data;
		d_X_current_data.resize(level_num + 1);
		for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
		{
			d_X_current_data[ln] = l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
		}
		
		ofstream ofile_test_connection;
		ofile_test_connection.open("test.out");
		std::vector<unsigned int>::iterator it_Pseudo;
		
		for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
		{

			boost::multi_array_ref<double, 2>& X_node = *d_X_current_data[ln]->getLocalFormVecArray();
			for (const auto& node_idx : local_nodes)
			{
				const int lag_idx = node_idx->getLagrangianIndex();
				if (plate2d_left_lag_idxs.first <= lag_idx && lag_idx < plate2d_left_lag_idxs.second)
				{
					it_Pseudo = std::find (Pseudo_MVleaflet_mapping_firstComponent.begin(), Pseudo_MVleaflet_mapping_firstComponent.end(), lag_idx);
					if (it_Pseudo != Pseudo_MVleaflet_mapping_firstComponent.end())
					{
						int location_inMapping = it_Pseudo-Pseudo_MVleaflet_mapping_firstComponent.begin();//starting from 0
						const int local_idx = node_idx->getLocalPETScIndex();
						double* const X = &X_node[local_idx][0];
						if (location_inMapping<=251)// first 252 mappings are for MV chordae-leaflet mapping
						{
							double du_x_connection,du_y_connection,du_z_connection;
							du_x_connection=MV_leaflets_MappingElement_currentLocation[location_inMapping*3]+Pseudo_MVleaflet_mapping[location_inMapping][2]-X[0];
							du_y_connection=MV_leaflets_MappingElement_currentLocation[location_inMapping*3+1]+Pseudo_MVleaflet_mapping[location_inMapping][3]-X[1];
							du_z_connection=MV_leaflets_MappingElement_currentLocation[location_inMapping*3+2]+Pseudo_MVleaflet_mapping[location_inMapping][4]-X[2];
							X[0] += du_x_connection;
							X[1] += du_y_connection;
							X[2] += du_z_connection;
						}
						else // The rest of the mapping are the AV leaflet-vessel mapping, here we use the same approach to update location.
						{				
							double du_x_connection,du_y_connection,du_z_connection;
							du_x_connection=MV_leaflets_MappingElement_currentLocation[location_inMapping*3]+Pseudo_MVleaflet_mapping[location_inMapping][2]-X[0];
							du_y_connection=MV_leaflets_MappingElement_currentLocation[location_inMapping*3+1]+Pseudo_MVleaflet_mapping[location_inMapping][3]-X[1];
							du_z_connection=MV_leaflets_MappingElement_currentLocation[location_inMapping*3+2]+Pseudo_MVleaflet_mapping[location_inMapping][4]-X[2];
							X[0] += 1.0*du_x_connection;
							X[1] += 1.0*du_y_connection;
							X[2] += 1.0*du_z_connection;				
						}
						
					}
					
				}
			}
			d_X_current_data[ln]->restoreArrays();
			
			//
			
		}
		ofile_test_connection.close();
		return;
	}
// Main driver routine that actually runs the simulation.
// Main routine that initializes the simulation libraries and calls the main
// driver routine.
int
main(
    int argc,
    char* argv[])
{
	//windkessel model needed parameters, note these also include ones for Pulmonary artery used in whole heart model. 
	std::vector<double> Q_old(2,0.0);
	std::vector<double> Q_current(2,0.0);
	std::vector<double> P_current(2,0.0);
	std::vector<double> P_predict(2,0.0);
	std::vector<double> P(4,0.0);
	std::vector<double> dQdt(2,0.0);
	double restart_flag = 0.0;
	//
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
	{    
    // Parse command line options, set some standard options from the input
    // file, initialize the restart database (if this is a restarted run), and
    // enable file logging.
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "FEMV.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // Get various standard options set in the input file.
    const bool dump_viz_data = app_initializer->dumpVizData();
    const int viz_dump_interval = app_initializer->getVizDumpInterval();
    const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();
    const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
    const string exodus_filename = app_initializer->getExodusIIFilename();
	

    const bool dump_restart_data = app_initializer->dumpRestartData();
    const int restart_dump_interval = app_initializer->getRestartDumpInterval();
    const int restart_restore_num = app_initializer->getRestartRestoreNumber();    
    const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
    const string restart_read_dirname = app_initializer->getRestartReadDirectory();

    const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
    const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
    const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
    if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
    {
        Utilities::recursiveMkdir(postproc_data_dump_dirname);
    }

    const bool dump_timer_data = app_initializer->dumpTimerData();
    const int timer_dump_interval = app_initializer->getTimerDumpInterval();

    // Configure the passive stress model.
    MechanicsModel::normalize_stress = input_db->getBool("NORMALIZE_STRESS");
	MechanicsModel::beta_s = input_db->getDouble("BETA_S");
	MechanicsModel::C1_s   = input_db->getDouble("C1_S");
    MechanicsModel::MV_start_cap_time = app_initializer->getComponentDatabase("BoundaryConditions")->getDouble("MV_start_cap_time");
	MechanicsModel::MV_start_annulus_time = app_initializer->getComponentDatabase("BoundaryConditions")->getDouble("MV_start_close_time");
    MechanicsModel::MV_start_close_time = app_initializer->getComponentDatabase("BoundaryConditions")->getDouble("MV_start_close_time");
    MechanicsModel::MV_end_close_time = app_initializer->getComponentDatabase("BoundaryConditions")->getDouble("MV_end_close_time");
    // Configure the active tension model.
    MechanicsModel::enable_active_tension = input_db->getBool("ENABLE_ACTIVE_TENSION");
    MechanicsModel::T_scale = input_db->getDouble("T_SCALE");
    // Configure the pressure loading.
    BoundaryConditions::P_load = input_db->getDouble("P_LOAD");
    BoundaryConditions::t_load = input_db->getDouble("T_LOAD");

    // Configure the penalty parameters.
    BoundaryConditions::kappa_s = input_db->getDouble("KAPPA_S");
    
    // Load and initialize the FE mesh.
    Mesh mesh(init.comm(),NDIM);
    ModelInitialization::initialize_mesh(mesh, input_db);

	
    // Create major IBAMR solver objects.
    Pointer<IBFEMethod> ibfe_method_ops
       =new IBFEMethod("IBFEMethod", app_initializer->getComponentDatabase("IBFEMethod"), &mesh,
                    app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           /*register_for_restart*/ true, restart_read_dirname, restart_restore_num); 
	//adding IB solver objects, Pseudo
    Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
    vector<Pointer<IBStrategy> > ib_ops_vec;
	ib_ops_vec.push_back(ib_method_ops);
    ib_ops_vec.push_back(ibfe_method_ops);
	Pointer<IBStrategySet> ib_ops_set = new IBStrategySet(ib_ops_vec.begin(), ib_ops_vec.end());
	//
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator("INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator("IBExplicitHierarchyIntegrator", app_initializer->getComponentDatabase("IBExplicitHierarchyIntegrator"), ib_ops_set, navier_stokes_integrator);

    // Create major SAMRAI algorithm and data objects.
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>("CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>("GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

	//creating a flowrate instance, Windkessel
	Windkessel Windkessel_model("Windkessel_model",input_db, app_initializer->getComponentDatabase("BoundaryConditions"), true);
    // Build equation systems that store various auxiliary variables, such as
    // the material axes and the values of spatially distributed active
    // contraction model state variables.
	//temp
	ibfe_method_ops->initializeFEEquationSystems();
    EquationSystems* equation_systems = ibfe_method_ops->getFEDataManager(0)->getEquationSystems();
    ModelInitialization::initialize_equation_systems(equation_systems,0);
	
	//Configure the IB solver, Pseudo
    Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer("IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
    ib_method_ops->registerLInitStrategy(ib_initializer);
    Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
    ib_force_fcn->registerSpringForceFunction(0, &chordae_spring_force);
    ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);
	
    // Configure the IBFE solver.
    pout << "Configuring the solver...\n";
	//mesh
    IBFEMethod::PK1StressFcnData PK1_dev_stress_data;
    PK1_dev_stress_data.fcn = MechanicsModel::PK1_dev_stress_function;
    MechanicsModel::get_PK1_dev_stress_function_systems(PK1_dev_stress_data.system_data);
    PK1_dev_stress_data.quad_type = QGAUSS;
    PK1_dev_stress_data.quad_order = THIRD;  // full-order integration
    ibfe_method_ops->registerPK1StressFunction(PK1_dev_stress_data,0);

    //mesh
	IBFEMethod::PK1StressFcnData PK1_dil_stress_data;
    PK1_dil_stress_data.fcn = MechanicsModel::PK1_dil_stress_function;
	MechanicsModel::get_PK1_dil_stress_function_systems(PK1_dil_stress_data.system_data);
    PK1_dil_stress_data.quad_type = QGAUSS;
    PK1_dil_stress_data.quad_order = FIRST;  // reduced-order integration
    if (MechanicsModel::normalize_stress) ibfe_method_ops->registerPK1StressFunction(PK1_dil_stress_data,0);
	
	//declare it anyway, just in case we need it in the future.
	IBFEMethod::LagSurfaceForceFcnData load_surface_force_data;
	load_surface_force_data.fcn = BoundaryConditions::tether_force_function;
    ibfe_method_ops->registerLagSurfaceForceFunction(load_surface_force_data,0);
	
	// use body force to tether structure
	IBFEMethod::LagBodyForceFcnData load_Body_force_data;
	load_Body_force_data.fcn = BoundaryConditions::tether_vol_force_function;
    ibfe_method_ops->registerLagBodyForceFunction(load_Body_force_data,0);
	

    // Set up a post-processor to reconstruct the various quantities of interest.
    Pointer<IBFEPostProcessor> ib_post_processor = new IBFECentroidPostProcessor("IBFEPostProcessor", ibfe_method_ops->getFEDataManager());

    // 1) Deformation gradient tensor FF = dX/ds.
    ib_post_processor->registerTensorVariable(
        "FF", MONOMIAL, CONSTANT, IBFEPostProcessor::FF_fcn);

    // 2) Deformed fiber and sheet axes.
    ib_post_processor->registerVectorVariable(
        "f", MONOMIAL, CONSTANT, IBFEPostProcessor::deformed_material_axis_fcn, MechanicsModel::f0_system_data);
	ib_post_processor->registerVectorVariable(
        "s", MONOMIAL, CONSTANT, IBFEPostProcessor::deformed_material_axis_fcn, MechanicsModel::s0_system_data);

    // 3) Fiber and sheet stretches.
    ib_post_processor->registerScalarVariable(
        "lambda_f", MONOMIAL, CONSTANT, IBFEPostProcessor::material_axis_stretch_fcn,MechanicsModel::f0_system_data);
	ib_post_processor->registerScalarVariable(
        "lambda_s", MONOMIAL, CONSTANT, IBFEPostProcessor::material_axis_stretch_fcn,MechanicsModel::s0_system_data);
		
    // 4) Cauchy stress sigma = (1/J) PP FF^T.
    std::pair<IBTK::TensorMeshFcnPtr,void*> PK1_dev_stress_fcn_data(&MechanicsModel::PK1_dev_stress_function,NULL);
    ib_post_processor->registerTensorVariable(
        "sigma_dev", MONOMIAL, CONSTANT,
        IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
        PK1_dev_stress_data.system_data, &PK1_dev_stress_fcn_data);

    std::pair<IBTK::TensorMeshFcnPtr,void*> PK1_dil_stress_fcn_data(&MechanicsModel::PK1_dil_stress_function,NULL);
    if (MechanicsModel::normalize_stress)
    {
        ib_post_processor->registerTensorVariable(
            "sigma_dil", MONOMIAL, CONSTANT,
            IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
            PK1_dil_stress_data.system_data, &PK1_dil_stress_fcn_data);
    }

// 5) Eulerian pressure p_f.
    Pointer<hier::Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
    Pointer<VariableContext> p_current_ctx = navier_stokes_integrator->getCurrentContext();
    HierarchyGhostCellInterpolation::InterpolationTransactionComponent p_ghostfill(
        /*data_idx*/ -1, "LINEAR_REFINE", /*use_cf_bdry_interpolation*/ false, "CONSERVATIVE_COARSEN", "LINEAR");
    FEDataManager::InterpSpec p_interp_spec("PIECEWISE_LINEAR", QGAUSS, FIFTH, /*use_adaptive_quadrature*/ false,
                                            /*point_density*/ 2.0, /*use_consistent_mass_matrix*/ true, /*use_nodal_quadrature*/ false);
    ib_post_processor->registerInterpolatedScalarEulerianVariable("p_f", LAGRANGE, FIRST, p_var, p_current_ctx,
                                                                  p_ghostfill, p_interp_spec);



    // Create Eulerian boundary condition specification objects when needed.															 																		  
	
	const string Heart_pres_data_file_name = app_initializer->getComponentDatabase("BoundaryConditions")->getStringWithDefault("Heart_Pressure_data_file_name","");
	if (Heart_pres_data_file_name.empty())
    {
        pout << "WARNING: no inlet pressure data specified!\n";
    }
    else
    {
        pout << "reading in inlet pressure boundary condition data from filename: " << Heart_pres_data_file_name << "\n\n";
    }

    double Heart_pres_data_delta_t = 0.0;
    std::vector<std::vector<double>> Heart_pres_data;
    if (!Heart_pres_data_file_name.empty())
    {
            // Expected pressure file format:
            // line 1: number of pressure records
            // line 2: the timestep between records (in seconds)
            // remaining lines: the pressures (in mmHg)
            ifstream is;
            is.open(Heart_pres_data_file_name.c_str(), ios::in);
            int num_records;
            is >> num_records;
            double delta_t;
            is >> delta_t;
            
            Heart_pres_data_delta_t = delta_t;
            Heart_pres_data.resize(num_records);
            for (int k = 0; k < num_records; ++k)
            {
                double pres_RSPV,pres_RIPV,pres_LSPV,pres_LIPV,pres_Aorta,pres_SCV,pres_ICV,pres_PA;
                is >> pres_RSPV >> pres_RIPV >> pres_LSPV >> pres_LIPV >> pres_Aorta >> pres_SCV >> pres_ICV >> pres_PA;
                pres_RSPV *= 1333.2239;  // convert from mmHg to Dyne/cm^2    
				pres_RIPV *= 1333.2239;
				pres_LSPV *= 1333.2239;
				pres_LIPV *= 1333.2239;
				pres_Aorta *= 1333.2239;
				pres_SCV *= 1333.2239;
				pres_ICV *= 1333.2239;
				pres_PA *= 1333.2239;				
				Heart_pres_data[k].resize(8);
                Heart_pres_data[k][0] = pres_RSPV ;
				Heart_pres_data[k][1] = pres_RIPV ;
				Heart_pres_data[k][2] = pres_LSPV ;
				Heart_pres_data[k][3] = pres_LIPV ;
				Heart_pres_data[k][4] = pres_Aorta ;
				Heart_pres_data[k][5] = pres_SCV ;
				Heart_pres_data[k][6] = pres_ICV ;
				Heart_pres_data[k][7] = pres_PA ;
			}
			is.close();
    }
    else
    {
            Heart_pres_data_delta_t = 0.0;
            Heart_pres_data.resize(1);
			Heart_pres_data[0].resize(1,0.0);
    }
	
	std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
	u_bc_coefs[0] = new VelocityBcCoefs("VelocityBcCoefs_0", grid_geometry, 0, Windkessel_model,Heart_pres_data_delta_t, Heart_pres_data, app_initializer->getComponentDatabase("BoundaryConditions"));
    u_bc_coefs[1] = new VelocityBcCoefs("VelocityBcCoefs_1", grid_geometry, 1, Windkessel_model,Heart_pres_data_delta_t, Heart_pres_data, app_initializer->getComponentDatabase("BoundaryConditions"));
    u_bc_coefs[2] = new VelocityBcCoefs("VelocityBcCoefs_2", grid_geometry, 2, Windkessel_model,Heart_pres_data_delta_t, Heart_pres_data, app_initializer->getComponentDatabase("BoundaryConditions"));
    navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
	
    //feedback forcer
    Pointer<FeedbackForcer> feedback_forcer = new FeedbackForcer("FeedbackForcer", patch_hierarchy);
    time_integrator->registerBodyForceFunction(feedback_forcer);
	

	
    // Set up visualization plot file writers.
    Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
	Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();//Pseudo
    if (uses_visit)
    {
		ib_initializer->registerLSiloDataWriter(silo_data_writer);//Pseudo
        time_integrator->registerVisItDataWriter(visit_data_writer);
		ib_method_ops->registerLSiloDataWriter(silo_data_writer);//Pseudo
    }
    std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

    // Initialize FE data.
    ibfe_method_ops->initializeFEData();
    ib_post_processor->initializeFEData();
	if (uses_exodus)
			{
				const bool from_restart = RestartManager::getManager()->isFromRestart();
				exodus_io->append(from_restart);
			}
    // Setup the material axes.
    ModelInitialization::initialize_material_axes(mesh, equation_systems, input_db,0);

    // Initialize hierarchy configuration and data on all patches.
    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

	//Pseudo
	const string Pseudo_MVleaflet_file_name = app_initializer->getComponentDatabase("BoundaryConditions")->getStringWithDefault("Pseudo_MVleaflet_file_name","");
	//
    // Deallocate initialization objects.
	ib_method_ops->freeLInitStrategy();//Pseudo
    ib_initializer.setNull();//Pseudo
    app_initializer.setNull();

    // Print the input database contents to the log file.
    plog << "Input database:\n";
    input_db->printClassData(plog);

    // Write out initial visualization data.
    int iteration_num = time_integrator->getIntegratorStep();
    double loop_time = time_integrator->getIntegratorTime();
	//temp
	if (MechanicsModel::enable_active_tension)
        {
			MechanicsModel::update_active_tension_variables(equation_systems, loop_time, loop_time);
		
        }	
    if (dump_viz_data)
    {
        pout << "\nWriting visualization files...\n\n";
        if (uses_visit)
        {
            time_integrator->setupPlotData();

            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);

			silo_data_writer->writePlotData(iteration_num, loop_time);//Pseudo

        }
        if (uses_exodus)
        {
            ib_post_processor->postProcessData(loop_time);
            exodus_io->write_timestep(exodus_filename, *equation_systems, iteration_num/viz_dump_interval+1, loop_time);
        }
    }
	
	
	// Setup feedback forcer object.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    feedback_forcer->d_U_data_idx = var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(), time_integrator->getCurrentContext());
	//Import mesh to label chamber cavity needed for Energy budget analysis, LV first 
	BoundaryConditions::readingPointsGeneral(mesh, BoundaryConditions::LV_endo_points_list, BoundaryConditions::LV_endo_points,BoundaryConditions::LV_NoOfEndoNode,"LV_POINTS_LIST.txt");
	BoundaryConditions::updatePointsPositionGeneralLV(equation_systems, BoundaryConditions::LV_endo_points_list, BoundaryConditions::LV_endo_points,BoundaryConditions::LV_NoOfEndoNode);
	std::vector<double> LVVector_Points;
	LVVector_Points.resize(BoundaryConditions::LV_NoOfEndoNode*3);
	if (0 == mesh.processor_id()) 
	{
		for (unsigned int i = 0; i< BoundaryConditions::LV_NoOfEndoNode; i++)
	   {
		      LVVector_Points[i*3]=BoundaryConditions::LV_endo_points[i][0];
			  LVVector_Points[i*3+1]=BoundaryConditions::LV_endo_points[i][1];
			  LVVector_Points[i*3+2]=BoundaryConditions::LV_endo_points[i][2];
	   }	
	}	
	MPI_Bcast(&LVVector_Points[0], BoundaryConditions::LV_NoOfEndoNode*3, MPI_DOUBLE, 0, init.comm().get());
	//Here read in the LV boundary facelets
	std::ifstream ifsendo("LV_Facelets_LIST.txt");
	int LVVexID,LVNoOfVex;	   
	ifsendo>>LVNoOfVex;//number of facelet is NoOfVex/3;
	std::vector<int> vv;
	std::vector<int>& LVVex_list=vv;
	LVVex_list.resize(LVNoOfVex);
	unsigned int Vex_pIndex = 0;	   
	   while (!ifsendo.eof()& Vex_pIndex<LVNoOfVex){
			ifsendo>>LVVexID;				
		    LVVex_list[Vex_pIndex]=LVVexID; 			
			Vex_pIndex = Vex_pIndex + 1;
			
		}
	//then for LA
	BoundaryConditions::readingPointsGeneral(mesh, BoundaryConditions::LA_endo_points_list, BoundaryConditions::LA_endo_points,BoundaryConditions::LA_NoOfEndoNode,"LA_POINTS_LIST.txt");
	BoundaryConditions::updatePointsPositionGeneralLA(equation_systems, BoundaryConditions::LA_endo_points_list, BoundaryConditions::LA_endo_points,BoundaryConditions::LA_NoOfEndoNode);
	std::vector<double> LAVector_Points;
	LAVector_Points.resize(BoundaryConditions::LA_NoOfEndoNode*3);
	if (0 == mesh.processor_id()) 
	{
		for (unsigned int i = 0; i< BoundaryConditions::LA_NoOfEndoNode; i++)
	   {
		      LAVector_Points[i*3]=BoundaryConditions::LA_endo_points[i][0];
			  LAVector_Points[i*3+1]=BoundaryConditions::LA_endo_points[i][1];
			  LAVector_Points[i*3+2]=BoundaryConditions::LA_endo_points[i][2];
	   }	
	}	
	MPI_Bcast(&LAVector_Points[0], BoundaryConditions::LA_NoOfEndoNode*3, MPI_DOUBLE, 0, init.comm().get());
	//Here read in the LA boundary facelets
	std::ifstream ifsendoLA("LA_Facelets_LIST.txt");
	int LAVexID,LANoOfVex;	   
	ifsendoLA>>LANoOfVex;//number of facelet is NoOfVex/3;
	std::vector<int> LAvv;
	std::vector<int>& LAVex_list=LAvv;
	LAVex_list.resize(LANoOfVex);
	unsigned int LAVex_pIndex = 0;	   
	   while (!ifsendoLA.eof()& LAVex_pIndex<LANoOfVex){
			ifsendoLA>>LAVexID;				
		    LAVex_list[LAVex_pIndex]=LAVexID; 			
			LAVex_pIndex = LAVex_pIndex + 1;
			
		}		
	//Finish Import mesh to label chamber cavity	
    // Main time step loop.
	ofstream ofile;
	ofile.open("test.out");
	ofile.close();
    pout << "Entering main time step loop...\n";
    const double loop_time_end = time_integrator->getEndTime();
    double dt = 5e-6*2.0;
	feedback_forcer->d_kappa = navier_stokes_integrator->getStokesSpecifications()->getRho()/dt;
	//////////////////////////////////Read the pseudo fibre and leaflet connection mappaing [chordae end node ID, FE mesh element ID]
	int total_mapping = 1;

	std::vector<std::vector<double>> Pseudo_MVleaflet_mapping;
	std::vector<unsigned int> Pseudo_MVleaflet_mapping_firstComponent;
	std::vector<unsigned int> Pseudo_MVleaflet_mapping_secondComponent;
	std::vector<double> MV_leaflets_MappingElement_currentLocation;
	if (!Pseudo_MVleaflet_file_name.empty())
    {
            ifstream is;
            is.open(Pseudo_MVleaflet_file_name.c_str(), ios::in);
            int mapping_number;
			int useless, useless2,useless3,useless4;
            is >> mapping_number >> useless>> useless2>> useless3>> useless4;

            
            total_mapping = mapping_number;
            Pseudo_MVleaflet_mapping.resize(total_mapping);
			Pseudo_MVleaflet_mapping_firstComponent.resize(total_mapping);
			Pseudo_MVleaflet_mapping_secondComponent.resize(total_mapping);
			MV_leaflets_MappingElement_currentLocation.resize(total_mapping*3,0.0);
            for (int k = 0; k < mapping_number; ++k)
            {
				Pseudo_MVleaflet_mapping[k].resize(5);
                int pseudo_ndid;
				int MV_eleid;
				double du_x,du_y,du_z;
                is >> pseudo_ndid >> MV_eleid>> du_x>> du_y>> du_z;
				
                Pseudo_MVleaflet_mapping[k][0] = pseudo_ndid*1.0;
				Pseudo_MVleaflet_mapping[k][1] = MV_eleid*1.0;
				Pseudo_MVleaflet_mapping[k][2] = du_x;
				Pseudo_MVleaflet_mapping[k][3] = du_y;
				Pseudo_MVleaflet_mapping[k][4] = du_z;
				Pseudo_MVleaflet_mapping_firstComponent[k] = pseudo_ndid;
				Pseudo_MVleaflet_mapping_secondComponent[k] = MV_eleid;
				
			}
			is.close();
    }
    else
    {
            Pseudo_MVleaflet_mapping.resize(1);
			Pseudo_MVleaflet_mapping_firstComponent.resize(1);
			Pseudo_MVleaflet_mapping_secondComponent.resize(1);
			Pseudo_MVleaflet_mapping[0].resize(5);
			Pseudo_MVleaflet_mapping[0][0]=1.0;
			Pseudo_MVleaflet_mapping[0][1]=1.0;
			Pseudo_MVleaflet_mapping[0][2]=1.0;
			Pseudo_MVleaflet_mapping[0][3]=1.0;
			Pseudo_MVleaflet_mapping[0][4]=1.0;
			Pseudo_MVleaflet_mapping_firstComponent[0] = 1;
			Pseudo_MVleaflet_mapping_secondComponent[0] = 1;
			MV_leaflets_MappingElement_currentLocation.resize(1*3,0.0);
    }
	
	///////////////////Finish Read the mapping 
	BoundaryConditions::boundary_info = mesh.boundary_info.get();
	int Energy_budget_analysis_flag=0;////////////////////////////////////////////////////////////////////Energy budget analysis
    while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) && time_integrator->stepsRemaining())
    {
        iteration_num = time_integrator->getIntegratorStep() ;

        pout <<                                                       endl;
        pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout << "At beginning of timestep # " << iteration_num  << endl;
        pout << "Simulation time is " << loop_time                 << endl;


        dt = time_integrator->getMaximumTimeStepSize();
		//Update Windkessel model at aorta tube boundary
		const int U_current_idx = var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(), navier_stokes_integrator->getCurrentContext());
		const int P_current_idx = var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(), navier_stokes_integrator->getCurrentContext());
		const int wgt_cc_idx = navier_stokes_integrator->getHierarchyMathOps()->getCellWeightPatchDescriptorIndex();
		const int wgt_sc_idx = navier_stokes_integrator->getHierarchyMathOps()->getSideWeightPatchDescriptorIndex();		
		Windkessel_model.advanceTimeDependentData(Q_current,P,loop_time, dt, patch_hierarchy, U_current_idx, P_current_idx, wgt_cc_idx, wgt_sc_idx,iteration_num);
		SAMRAI_MPI::sumReduction(&Q_current[0],2);	
		if (P_predict[0] <= 10.0 && restart_flag <= 1.0 && loop_time >= 0.1)// this is to for restarting step
		{
			SAMRAI_MPI::sumReduction(&P[0],4);
			P_current[0]=P[0]/P[1];
			//P_current[1]=P[2]/P[3];//this was for Pulmonary artery in whole heart model 
			Q_old[0]=Q_current[0];
			//Q_old[1]=Q_current[1];//this was for Pulmonary artery in whole heart model
			restart_flag = 1000;
		}
		if (loop_time <= 1.0)
		{
			P_current[0] = 50.0*1333.32;//ramp pressure
			//P_current[1] = 8.0*1333.32;//this was for Pulmonary artery in whole heart model
		}
		Windkessel_model.predictPressure(Q_old,Q_current,dQdt,P_current,P_predict,dt,iteration_num);
		VelocityBcCoefs::p_new = P_predict;
		P_current = P_predict;
		//Finish Update Windkessel model at aorta tube boundary
       
        MechanicsModel::I1_dev_max = -1.0e8;
        MechanicsModel::I1_dev_min = +1.0e8;
        MechanicsModel::I1_dil_max = -1.0e8;
        MechanicsModel::I1_dil_min = +1.0e8;

        MechanicsModel::J_dev_max = -1.0e8;
        MechanicsModel::J_dev_min = +1.0e8;
        MechanicsModel::J_dil_max = -1.0e8;
        MechanicsModel::J_dil_min = +1.0e8;
		
		//Update chamber cavity surface mesh location for Energy budget analysis				
		Energy_budget_analysis_flag = 0;
		if (dump_viz_data && ((iteration_num+1)%viz_dump_interval == 0))
		{
			Energy_budget_analysis_flag=1;
			//////////update LV
			BoundaryConditions::updatePointsPositionGeneralLV(equation_systems, 
										BoundaryConditions::LV_endo_points_list, 
										BoundaryConditions::LV_endo_points,
										BoundaryConditions::LV_NoOfEndoNode);	
			if (0 == mesh.processor_id()) 
			{
				for (unsigned int i = 0; i< BoundaryConditions::LV_NoOfEndoNode; i++)
			   {
					  LVVector_Points[i*3]=BoundaryConditions::LV_endo_points[i][0];
					  LVVector_Points[i*3+1]=BoundaryConditions::LV_endo_points[i][1];
					  LVVector_Points[i*3+2]=BoundaryConditions::LV_endo_points[i][2];
			   }	
			}	
			MPI_Bcast(&LVVector_Points[0], BoundaryConditions::LV_NoOfEndoNode*3, MPI_DOUBLE, 0, init.comm().get());
			MPI_Barrier(init.comm().get());
			////////update LA
			BoundaryConditions::updatePointsPositionGeneralLA(equation_systems, 
										BoundaryConditions::LA_endo_points_list, 
										BoundaryConditions::LA_endo_points,
										BoundaryConditions::LA_NoOfEndoNode);
			if (0 == mesh.processor_id()) 
			{
				for (unsigned int i = 0; i< BoundaryConditions::LA_NoOfEndoNode; i++)
			   {
					  LAVector_Points[i*3]=BoundaryConditions::LA_endo_points[i][0];
					  LAVector_Points[i*3+1]=BoundaryConditions::LA_endo_points[i][1];
					  LAVector_Points[i*3+2]=BoundaryConditions::LA_endo_points[i][2];
			   }	
			}
			
			MPI_Bcast(&LAVector_Points[0], BoundaryConditions::LA_NoOfEndoNode*3, MPI_DOUBLE, 0, init.comm().get());
			MPI_Barrier(init.comm().get());
		}
		navier_stokes_integrator->inputHeartMesh(LVVex_list,LVVector_Points,LVNoOfVex,LAVex_list,LAVector_Points,LANoOfVex,Energy_budget_analysis_flag);//
		//Finish Update chamber cavity surface mesh location
		
        time_integrator->advanceHierarchy(dt);
		//Calculate the FE element location for pseudo-FE mesh mapping		
		System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
		NumericVector<double>* X_vec = X_system.solution.get();
		NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
		X_vec->localize(*X_ghost_vec);
		DofMap& X_dof_map = X_system.get_dof_map();
		std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
		std::unique_ptr<FEBase> fe(FEBase::build(NDIM, X_dof_map.variable_type(0)));
		std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, NDIM, FIFTH);
		fe->attach_quadrature_rule(qrule.get());
		const std::vector<double>& JxW = fe->get_JxW();
		const std::vector<std::vector<VectorValue<double> > >& dphi = fe->get_dphi();
		TensorValue<double> FF;
		boost::multi_array<double, 2> X_node;
		const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
		const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
		//Prepare for iteration
		std::vector<unsigned int>::iterator it_MV;

		for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
		{
			Elem* const elem = *el_it;
			fe->reinit(elem);
			if (elem->subdomain_id() == 1000 || elem->subdomain_id() == 30|| elem->subdomain_id() == 40)// all FE mapped element should belong these three element sets
			{
				for (unsigned int d = 0; d < NDIM; ++d)
				{
					X_dof_map.dof_indices(elem, X_dof_indices[d], d);
				}
				const int n_qp = qrule->n_points();
				get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);
				
				const int elem_id = elem->id();			
				for(unsigned int ii = 0; ii < Pseudo_MVleaflet_mapping_secondComponent.size(); ++ii)
				{
					if (Pseudo_MVleaflet_mapping_secondComponent[ii] == elem_id)
					{
						MV_leaflets_MappingElement_currentLocation[ii*3]=0.25*(X_node[0][0]+X_node[1][0]+X_node[2][0]+X_node[3][0]);
						MV_leaflets_MappingElement_currentLocation[ii*3+1]=0.25*(X_node[0][1]+X_node[1][1]+X_node[2][1]+X_node[3][1]);
						MV_leaflets_MappingElement_currentLocation[ii*3+2]=0.25*(X_node[0][2]+X_node[1][2]+X_node[2][2]+X_node[3][2]);	
					}
				}
			}

		}
		std::vector<double> MV_leaflets_MappingElement_currentLocation_sum(MV_leaflets_MappingElement_currentLocation.size(),0.0);
		for (unsigned int i = 0; i < MV_leaflets_MappingElement_currentLocation.size(); ++i)
		{
			MV_leaflets_MappingElement_currentLocation_sum[i]=SAMRAI_MPI::sumReduction(MV_leaflets_MappingElement_currentLocation[i]);
		}
		
		MPI_Bcast(&MV_leaflets_MappingElement_currentLocation_sum[0], MV_leaflets_MappingElement_currentLocation_sum.size(), MPI_DOUBLE, 0, init.comm().get());
		update_target_point_positions(patch_hierarchy, ib_method_ops->getLDataManager(), loop_time, dt, equation_systems, mesh,Pseudo_MVleaflet_mapping,Pseudo_MVleaflet_mapping_firstComponent, Pseudo_MVleaflet_mapping_secondComponent,MV_leaflets_MappingElement_currentLocation_sum);//pseudo
		update_connection_point_positions(patch_hierarchy, ib_method_ops->getLDataManager(), loop_time, dt,Pseudo_MVleaflet_mapping,Pseudo_MVleaflet_mapping_firstComponent, Pseudo_MVleaflet_mapping_secondComponent,MV_leaflets_MappingElement_currentLocation_sum);//pseudo
		//Finish Calculate the FE element location for pseudo-FE mesh mapping		
		

        pout << "I1_dev_max = " << MechanicsModel::I1_dev_max << "\n"
             << "I1_dev_min = " << MechanicsModel::I1_dev_min << "\n"
             << "I1_dil_max = " << MechanicsModel::I1_dil_max << "\n"
             << "I1_dil_min = " << MechanicsModel::I1_dil_min << "\n";

        pout << "J_dev_max = " << MechanicsModel::J_dev_max << "\n"
             << "J_dev_min = " << MechanicsModel::J_dev_min << "\n"
             << "J_dil_max = " << MechanicsModel::J_dil_max << "\n"
             << "J_dil_min = " << MechanicsModel::J_dil_min << "\n";
		if (MechanicsModel::enable_active_tension)
        {
            ActiveContraction::update_active_tension_model_state_variables(equation_systems, loop_time, dt,0);
			MechanicsModel::update_active_tension_variables(equation_systems, loop_time, dt);
        }	 

        loop_time += dt;

        pout <<                                                       endl;
        pout << "At end       of timestep # " << iteration_num  << endl;
        pout << "Simulation time is " << loop_time                 << endl;
        pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        pout <<                                                       endl;

        // At specified intervals, write visualization and restart files,
        // print out timer data, and store hierarchy data for post
        // processing.
        iteration_num += 1;
		const bool last_step = !time_integrator->stepsRemaining();

        if (dump_viz_data && (iteration_num%viz_dump_interval == 0 || last_step ))
        {
            pout << "\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
				silo_data_writer->writePlotData(iteration_num, loop_time);//pseudo
            }
            if (uses_exodus)
            {
                ib_post_processor->postProcessData(loop_time);
                exodus_io->write_timestep(exodus_filename, *equation_systems, iteration_num/viz_dump_interval+1, loop_time);
	
            }
        }
        if (dump_restart_data && (iteration_num%restart_dump_interval == 0 || last_step))
        {
            pout << "\nWriting restart files...\n\n";
            RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
			ibfe_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);   
        }
        if (dump_timer_data && (iteration_num%timer_dump_interval == 0 || last_step ))
        {
            pout << "\nWriting timer data...\n\n";
            TimerManager::getManager()->print(plog);
        }
    }

    // Clean up dynamically allocated objects.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        delete u_bc_coefs[d];
    }
	
}// main_driver
    SAMRAIManager::shutdown();
    return 0;
}// main


