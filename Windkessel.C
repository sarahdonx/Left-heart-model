#include "Windkessel.h"
#include "BoundaryConditions.h"
#include <math.h>
//#include "tools.h"
/////////////////////////////// INCLUDES /////////////////////////////////////


#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>
// SAMRAI INCLUDES
//#include <CartesianGridGeometry.h>
//#include <CartesianPatchGeometry.h>
//#include <PatchLevel.h>
//#include <SideData.h>
//#include <tbox/SAMRAI_MPI.h>
#include <tbox/RestartManager.h>
//#include <tbox/Utilities.h>

#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
// BLITZ++ INCLUDES
//#include <blitz/array.h>
#include <boost/array.hpp>

// C++ STDLIB INCLUDES
#include <cassert>
#include <numeric>
//
//
#include <stdlib.h>
#include <stdio.h>
namespace
{
//conversion factors
static const double prconv = 1333.2239; //mmHg ====> dyne/cm^2
static const double R = 1.13*prconv; //resistance (dyne ml^-1 s)
static const double C = 2.57/prconv; //total arterial compliance (ml dyne^-1)
static const double Z_c = 0.071*prconv;//impedance (dyne ml^-1 s)



}


//////////////////////////////PUBLIC//////////////////////
Windkessel::Windkessel(
    const string& object_name,
    Pointer<Database> input_db,
	const tbox::Pointer<tbox::Database> input_db_bc,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_time(0.0),
      d_psrc(2,0.0),
      d_rsrc(0.0),
      d_posn(3),
	  d_P_PV_current(4,0.0),
	  d_P_PV_predict(4,0.0)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
#endif
	d_posn[0] = -12.61;
	d_posn[1] = -14.70;
	d_posn[2] = 12.69;
	d_rsrc = 0.5;
}


Windkessel::~Windkessel()
{
    return;
}// ~Windkessel


void Windkessel::advanceTimeDependentData(   
    std::vector<double> &Q_current,
	std::vector<double> &P_current,   	
	const double d_time_input,
    const double dt,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int U_idx,
    const int P_idx,
    const int wgt_cc_idx,
    const int wgt_sc_idx,
	const int iteration_num)
{
	std::vector<double> Q(4,0.0);
	std::vector<double> P(4,0.0);
	for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
		{
			Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
			for (PatchLevel<NDIM>::Iterator p(level); p; p++)
				{
					Pointer<Patch<NDIM> > patch = level->getPatch(p());
					Pointer<SideData<NDIM,double> > U_data = patch->getPatchData(U_idx);
					Pointer<SideData<NDIM,double> > wgt_sc_data = patch->getPatchData(wgt_sc_idx);
					Pointer<CellData<NDIM,double> > P_data = patch->getPatchData(P_idx);
					Pointer<CellData<NDIM,double> > wgt_cc_data = patch->getPatchData(wgt_cc_idx);
					const Box<NDIM>& patch_box = patch->getBox();
					const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
					Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
					const double* const XLower = pgeom->getXLower();
					const double* const dx = pgeom->getDx();
					Pointer<CartesianGridGeometry<NDIM> > grid_geometry = hierarchy->getGridGeometry();
					const double* const dx_coarsest = grid_geometry->getDx();
					double dx_target[NDIM];
					const int target_ln = hierarchy->getFinestLevelNumber();
					const IntVector<NDIM>& target_ratio = hierarchy->getPatchLevel(target_ln)->getRatio();
					for (int d = 0; d < NDIM; ++d)
						{
							dx_target[d] = dx_coarsest[d] / double(target_ratio(d));
						}

					const double L[NDIM] = { max(2.0*dx_target[0] , dx_coarsest[0]),
											 max(2.0*dx_target[1] , dx_coarsest[1]),
											 max(2.0*dx_target[2] , dx_coarsest[2]) };

					const int offset[NDIM] = { int(L[0]/dx[0])-1 ,
											   int(L[1]/dx[1])-1 ,
											   int(L[2]/dx[2])-1 };

					const IntVector<NDIM>& ratio = pgeom->getRatio();
					const Box<NDIM> domain_box = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0],ratio);
					const SAMRAI::hier::Index<NDIM>& domain_lower = domain_box.lower();
					const SAMRAI::hier::Index<NDIM>& domain_upper = domain_box.upper();
					for (int bdry_normal_axis = 0; bdry_normal_axis < NDIM; ++bdry_normal_axis)
						{
							for (int side = 0; side <= 1; ++side)
							{
								if (pgeom->getTouchesRegularBoundary(bdry_normal_axis,side))
								{
									Box<NDIM> bdry_box = domain_box;
									if (side == 0)
									{
										bdry_box.upper()(bdry_normal_axis) = domain_lower(bdry_normal_axis);//+offset[bdry_normal_axis];
									}
									else
									{
										bdry_box.lower()(bdry_normal_axis) = domain_upper(bdry_normal_axis);//-offset[bdry_normal_axis];
									}
									for (int component = 0; component < NDIM; ++component)
									{
										if (component == bdry_normal_axis)
										{
											if (bdry_normal_axis == 2) //here we compute the flow integration
											{
												//for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box*patch_box,component)); b; b++)
												for (Box<NDIM>::Iterator b(bdry_box*patch_box); b; b++)
												{
													if (side == 1)//compute at the upper tube circle,Aorta [-0.777,-21.290,148.407]
													{
														const SAMRAI::hier::Index<NDIM>& i = b();
														const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Upper);
														double X[NDIM];
														for (int d = 0; d < NDIM; ++d)
														{
															X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+(d == component ? 1.0 : 0.5));
														}
														const double r =  sqrt((X[0]+0.81)*(X[0]+0.81) + (X[1]+21.32)*(X[1]+21.32));
														if (r <= 2.5)
														{
															Q[0] +=  (*U_data)(i_s)*dx[0]*dx[1];
															Q[1] +=  dx[0]*dx[1];			
														}
														if (r <= 1.0)// this P is only used in restart steps
														{
															P[0] +=  (*P_data)(i)*dx[0]*dx[1]*dx[2];
															P[1] +=  dx[0]*dx[1]*dx[2];
														}
													}
												}
												
											}
										}
									}
								}
							}
						}										   
				}
		}
	Q_current[0] =	Q[0];
	P_current[0] = P[0];
	P_current[1] = P[1];
	Q_current[1] =	Q[2];
	P_current[2] = P[2];
	P_current[3] = P[3];
	return;
	}
	
void Windkessel::predictPressure(   
    std::vector<double> &Q_old,
	std::vector<double> &Q_current,
	std::vector<double> &dQdt,
	std::vector<double> &P_current,
	std::vector<double> &P_predict,	
	const double dt,	
	const int iteration_num)
	{
		double fq = 0.0;
		double aorta_R = 0.6*R; 
		double aorta_C = 1.33*C; 
		double aorta_ZC = 0.4*Z_c; 
		dQdt[0] = (Q_current[0]-Q_old[0])/dt;
		Q_old[0] = Q_current[0];		
		fq = Q_current[0]*(1+aorta_ZC/aorta_R)+aorta_ZC*aorta_C*dQdt[0];
		P_predict[0] = (P_current[0]+dt*fq/aorta_C)/(1+dt/(aorta_R*aorta_C));
		return;
	}
	
void
Windkessel::putToDatabase(
    Pointer<Database> db)
{
    //db->putDouble("d_time",d_time);
    //db->putDoubleArray("d_psrc",&d_psrc[0],2);
	//db->putDoubleArray("d_posn",&d_posn[0],3);
	//db->putDoubleArray("d_Q_PV_current",&d_Q_PV_current[0],4);
	//db->putDoubleArray("d_Q_PV_predict",&d_Q_PV_predict[0],4);
    //db->putDoubleArray("d_rsrc",d_rsrc);
    return;
}// putToDatabase

