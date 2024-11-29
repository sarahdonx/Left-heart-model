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
#include <FeedbackForcer.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <SideData.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <numeric>

// CLASS IMPLEMENTATION

namespace
{
inline double
delta_4(
    const double& r)
{
    const double r1 = abs(r);
    const double r2 = r*r;
    if (r1 < 1.0)
    {
        return 0.125*(3.0 - 2.0*r1 + sqrt( 1.0+ 4.0*r1-4.0*r2));
    }
    else if (r1 < 2.0)
    {
        return 0.125*(5.0 - 2.0*r1 - sqrt(-7.0+12.0*r1-4.0*r2));
    }
    else
    {
        return 0.0;
    }
}// delta_4
}

double FeedbackForcer::MV_start_close_time;
double FeedbackForcer::MV_end_close_time;

FeedbackForcer::FeedbackForcer(
    const string& object_name,
    const Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
    : d_U_data_idx(-1),
      d_object_name(object_name),
      d_patch_hierarchy(patch_hierarchy)
{
    return;
}// FeedbackForcer

FeedbackForcer::~FeedbackForcer()
{
    // intentionally blank
    return;
}// ~FeedbackForcer

void
FeedbackForcer::setDataOnPatch(
    const int data_idx,
    Pointer<hier::Variable<NDIM> >  /*var*/,
    Pointer<Patch<NDIM> > patch,
    const double data_time,
    const bool initial_time,
    Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    Pointer<SideData<NDIM,double> > F_data = patch->getPatchData(data_idx);
    F_data->fillAll(0.0);

    if (initial_time) return;

    Pointer<SideData<NDIM,double> > U_data = patch->getPatchData(d_U_data_idx);

    const Box<NDIM>& patch_box = patch->getBox();
    const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_patch_hierarchy->getGridGeometry();

    const double* const dx_coarsest = grid_geometry->getDx();
    double dx_target[NDIM];
    const int target_ln = d_patch_hierarchy->getFinestLevelNumber();
    const IntVector<NDIM>& target_ratio = d_patch_hierarchy->getPatchLevel(target_ln)->getRatio();
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
	double time_in_period = -1;
	if (data_time >= 1.0)
	{
		time_in_period = fmod(data_time-1.0, 0.8);
	}
	

    // WARNING: The following code will not compute the correct feedback force
    // near corners except in the case that the upper/lower z boundary has both
    // no-penetration and no-slip boundary conditions near edges/corners in the
    // computational domain.
    for (int bdry_normal_axis = 0; bdry_normal_axis < NDIM; ++bdry_normal_axis)
    {
        for (int side = 0; side <= 1; ++side)
        {
            if (pgeom->getTouchesRegularBoundary(bdry_normal_axis,side))
            {
                Box<NDIM> bdry_box = domain_box;
                if (side == 0)
                {
                    bdry_box.upper()(bdry_normal_axis) = domain_lower(bdry_normal_axis)+offset[bdry_normal_axis];
                }
                else
                {
                    bdry_box.lower()(bdry_normal_axis) = domain_upper(bdry_normal_axis)-offset[bdry_normal_axis];
                }

                for (int component = 0; component < NDIM; ++component)
                {
                    if (component != bdry_normal_axis)
                    {
                        // Penalize deviations from the specified tangential
                        // velocity boundary conditions.
                        for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box*patch_box,component)); b; b++)
                        {
                            const SideIndex<NDIM> i_s(b(), component, SideIndex<NDIM>::Lower);
                            (*F_data)(i_s) = d_kappa*(0.0 - (*U_data)(i_s));
                        }
                    }
					else if (bdry_normal_axis == 0)
					{
						// Penalize deviations from the specified normal
                        // velocity boundary conditions.
                        for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box*patch_box,component)); b; b++)
                        {
                            const SAMRAI::hier::Index<NDIM>& i = b();
                            const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                            double X[NDIM];
                            for (int d = 0; d < NDIM; ++d)
                            {
                                X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+(d == component ? 0.0 : 0.5));
                            }

                            // Compute the radius from the centerline of the
                            // valve tester to the current grid point.
							if (side ==0)
							{
								const double r_RIPV = sqrt((X[1]+12.896)*(X[1]+12.896) + (X[2]-144.008)*(X[2]-144.008));//thickness 0.15 cm
								const double r_RSPV = sqrt((X[1]+17.159)*(X[1]+17.159) + (X[2]-145.39)*(X[2]-145.39));//thickness 0.15 cm
								const double R = 1.58;
								const double h = dx[0];
								double radial_grader = 0.0;
								if (r_RIPV >= (R+0.6) && r_RSPV >= (R+0.6))
								{
									if (r_RIPV <= r_RSPV)
									{
										radial_grader = delta_4((R+0.6-r_RIPV)/h)/delta_4(0.0);
									}
									else
									{
										radial_grader = delta_4((R+0.6-r_RSPV)/h)/delta_4(0.0);
									}
								}
								else if (r_RIPV <= R)
								{
									radial_grader=abs((*U_data)(i_s))/4000.0;
									//radial_grader=0.0;
									if (time_in_period>= 0.0 && time_in_period <= 0.8)	
									{
										if ((*U_data)(i_s)<=0)
										{
											radial_grader=1.0;
										}										
									}
								}								
								else if (r_RSPV <= R)
								{
									radial_grader=abs((*U_data)(i_s))/4000.0;
									//radial_grader=0.0;
									if (time_in_period>= 0.0 && time_in_period <= 0.8)	
									{
										if ((*U_data)(i_s)<=0)
										{
											radial_grader=1.0;
										}										
									}
														
								}
								else
								{
									radial_grader = 1.0;
								}	
																
								(*F_data)(i_s) = radial_grader*d_kappa*(0.0 - (*U_data)(i_s));
							}
							
 
                        }//for (Box<NDIM>)
					}
				else if (bdry_normal_axis == 1)
				{
					// Penalize deviations from the specified normal
                        // velocity boundary conditions.
                        for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box*patch_box,component)); b; b++)
                        {
                            const SAMRAI::hier::Index<NDIM>& i = b();
                            const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                            double X[NDIM];
                            for (int d = 0; d < NDIM; ++d)
                            {
                                X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+(d == component ? 0.0 : 0.5));
                            }

                            // Compute the radius from the centerline of the
                            // valve tester to the current grid point.
							if (side ==1)
							{
								const double r_LIPV = sqrt((X[0]-5.782)*(X[0]-5.782) + (X[2]-144.527)*(X[2]-144.527));// thickness 0.15 cm
								
								const double R = 1.6;
								const double h = dx[0];
								double radial_grader = 0.0;
								if (r_LIPV >= (R+0.6))
								{
									radial_grader = delta_4((R+0.6-r_LIPV)/h)/delta_4(0.0);
								}
								else if (r_LIPV <= R)
								{
									radial_grader=abs((*U_data)(i_s))/4000.0;
									//radial_grader=0.0;
									if (time_in_period>= 0.0 && time_in_period <= 0.8)	
									{
										if ((*U_data)(i_s)>=0)
										{
											radial_grader=1.0;
										}										
									}									
								}		
								else
								{
									radial_grader = 1.0;	
								}
					
								(*F_data)(i_s) = radial_grader*d_kappa*(0.0 - (*U_data)(i_s));
							}
          
                        }//for (Box<NDIM>)
				}
				else if (bdry_normal_axis == 2)
				{
					// Penalize deviations from the specified normal
                        // velocity boundary conditions.
                        for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box*patch_box,component)); b; b++)
                        {
                            const SAMRAI::hier::Index<NDIM>& i = b();
                            const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                            double X[NDIM];
                            for (int d = 0; d < NDIM; ++d)
                            {
                                X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+(d == component ? 0.0 : 0.5));
                            }

                            // Compute the radius from the centerline of the
                            // valve tester to the current grid point.
							if (side ==0)
							{
								//do nothing
							}
							else if (side ==1)
							{
								
								const double r_LSPV = sqrt((X[0]-6.6)*(X[0]-6.6) + (X[1]+14.662)*(X[1]+14.662));
								const double r_Aorta = sqrt((X[0]+0.81)*(X[0]+0.81) + (X[1]+21.32)*(X[1]+21.32));// thickness 0.15 cm
								const double h = dx[0];
								double radial_grader = 0.0;
								double re_LSPV = r_LSPV-(1.6+0.6);
								double re_Aorta = r_Aorta-(2.58+0.7);
								if (r_LSPV >= (1.6+0.6) && r_Aorta >= (2.58+0.7))
								{								
									radial_grader = 1.0;
								}
								else if (r_LSPV <= 1.6)
								{
									radial_grader=abs((*U_data)(i_s))/4000.0;
									//radial_grader=0.0;
									if (time_in_period>= 0.0 && time_in_period <= 0.8)	
									{
										if ((*U_data)(i_s)>=0)
										{
											radial_grader=1.0;
										}										
									}									
								}	
								else if (r_Aorta <= 2.58)
								{
									radial_grader = 0.0;
									if (time_in_period>= 0.285 && time_in_period <= 0.8)	
									{
										radial_grader=abs((*U_data)(i_s))/2000.0;								
									}
									else if (data_time <= 1.0)
									{
										radial_grader=abs((*U_data)(i_s))/2000.0;
									}
									
									if ((*U_data)(i_s)<=0.0)
									{
										radial_grader=abs((*U_data)(i_s))/4000.0;
									}
																		
								}									
								else
								{
									radial_grader = 1.0;
								}																
								(*F_data)(i_s) = radial_grader*d_kappa*(0.0 - (*U_data)(i_s));
							}
          
                        }//for (Box<NDIM>)
				}
                }// for (int component = 0; component < NDIM; ++component)
            } //if (pgeom->getTouchesRegularBoundary(bdry_normal_axis,side))
        }//for (int side = 0; side <= 1; ++side)
    }//for (int bdry_normal_axis = 0; bdry_normal_axis < NDIM; ++bdry_normal_axis)
    return;
}// setDataOnPatch
