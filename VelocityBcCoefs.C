// Filename: VelocityBcCoefs.C
// Last modified: <11.Jun.2009 16:15:16 griffith@box230.cims.nyu.edu>
// Created on 04 May 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "VelocityBcCoefs.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <tbox/Utilities.h>

// NAMESPACE
using namespace IBTK;
double VelocityBcCoefs::P_dynamic; 
std::vector<double> VelocityBcCoefs::p_new; 							   
/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcCoefs::VelocityBcCoefs(
    const string& object_name,
    const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry,
    const int component_normal_axis,
	const Windkessel& Windkessel_model,
    const double& Heart_pres_data_delta_t,
    const vector<vector<double>>& Heart_pres_data,
    const tbox::Pointer<tbox::Database> input_db)
    : d_object_name(object_name),
      d_grid_geometry(grid_geometry),
      d_component_normal_axis(component_normal_axis),
	  d_Windkessel_model(Windkessel_model),
      d_radius(input_db->getDouble("radius")),
      d_init_time(input_db->getDouble("init_time")),
      d_Heart_pres_data_delta_t(Heart_pres_data_delta_t),
      d_Heart_pres_data(Heart_pres_data)
{
    // intentionally blank
    return;
}// VelocityBcCoefs

VelocityBcCoefs::~VelocityBcCoefs()
{
    // intentionally blank
    return;
}// ~VelocityBcCoefs

void
VelocityBcCoefs::setBcCoefs(
    tbox::Pointer<pdat::ArrayData<NDIM,double> >& acoef_data,
    tbox::Pointer<pdat::ArrayData<NDIM,double> >& bcoef_data,
    tbox::Pointer<pdat::ArrayData<NDIM,double> >& gcoef_data,
    const tbox::Pointer<hier::Variable<NDIM> >& variable,
    const hier::Patch<NDIM>& patch,
    const hier::BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
    // Compute P_inlet at fill_time.
	const double Mytime= fill_time-d_init_time;

	double RSPV_pressure_inlet = 0.0;
	double RIPV_pressure_inlet = 0.0;
	double LSPV_pressure_inlet = 0.0;
	double LIPV_pressure_inlet = 0.0;
	double Aorta_inlet = 0.0;
	double SCV_inlet = 0.0;
	double ICV_inlet = 0.0;
	double PA_inlet = 0.0;
    if (fill_time <= d_init_time)
    {
		RSPV_pressure_inlet = (fill_time/d_init_time)*d_Heart_pres_data[0][0];
		RIPV_pressure_inlet = (fill_time/d_init_time)*d_Heart_pres_data[0][1];
		LSPV_pressure_inlet = (fill_time/d_init_time)*d_Heart_pres_data[0][2];
		LIPV_pressure_inlet = (fill_time/d_init_time)*d_Heart_pres_data[0][3];
		Aorta_inlet = (fill_time/d_init_time)*d_Heart_pres_data[0][4];
		SCV_inlet = (fill_time/d_init_time)*d_Heart_pres_data[0][5];
		ICV_inlet = (fill_time/d_init_time)*d_Heart_pres_data[0][6];
		PA_inlet = (fill_time/d_init_time)*d_Heart_pres_data[0][7];

		//P_inlet = d_inpres_data[0];
    }
    else
    {
        const double t = fill_time-d_init_time;
        const double dt = d_Heart_pres_data_delta_t;
        const int k = int(floor(t/dt));
        if (k < 0)
        {
			RSPV_pressure_inlet = d_Heart_pres_data[0][0];
			RIPV_pressure_inlet = d_Heart_pres_data[0][1];
			LSPV_pressure_inlet = d_Heart_pres_data[0][2];
			LIPV_pressure_inlet = d_Heart_pres_data[0][3];
			Aorta_inlet = 		  d_Heart_pres_data[0][4];
			SCV_inlet = 		  d_Heart_pres_data[0][5];
			ICV_inlet = 		  d_Heart_pres_data[0][6];
			PA_inlet = 			  d_Heart_pres_data[0][7];			
				
        }                                                                                                                                                                                                                
        else if (k+1 >= int(d_Heart_pres_data.size()))                                                                                                                                                                       
        {                                                                                                                                                                                                                
			RSPV_pressure_inlet = d_Heart_pres_data[d_Heart_pres_data.size()-1][0];
			RIPV_pressure_inlet = d_Heart_pres_data[d_Heart_pres_data.size()-1][1];
			LSPV_pressure_inlet = d_Heart_pres_data[d_Heart_pres_data.size()-1][2];
			LIPV_pressure_inlet = d_Heart_pres_data[d_Heart_pres_data.size()-1][3];
			Aorta_inlet = 		  d_Heart_pres_data[d_Heart_pres_data.size()-1][4];
			SCV_inlet = 		  d_Heart_pres_data[d_Heart_pres_data.size()-1][5];
			ICV_inlet = 		  d_Heart_pres_data[d_Heart_pres_data.size()-1][6];
			PA_inlet = 			  d_Heart_pres_data[d_Heart_pres_data.size()-1][7];					
        }                                                                                                                                                                                                                
        else                                                                                                                                                                                                             
        {                                                                                                                                                                                                                
            const vector<double> f1 = d_Heart_pres_data[k  ];                                                                                                                                                                        
            const vector<double> f2 = d_Heart_pres_data[k+1];                                                                                                                                                                        
            RSPV_pressure_inlet = (f2[0]-f1[0])/dt*t+(f1[0]*dt-f2[0]*double(k)*dt+f1[0]*double(k)*dt)/dt;  
			RIPV_pressure_inlet = (f2[1]-f1[1])/dt*t+(f1[1]*dt-f2[1]*double(k)*dt+f1[1]*double(k)*dt)/dt; 
			LSPV_pressure_inlet = (f2[2]-f1[2])/dt*t+(f1[2]*dt-f2[2]*double(k)*dt+f1[2]*double(k)*dt)/dt; 
			LIPV_pressure_inlet = (f2[3]-f1[3])/dt*t+(f1[3]*dt-f2[3]*double(k)*dt+f1[3]*double(k)*dt)/dt; 	
			Aorta_inlet = 		(f2[4]-f1[4])/dt*t+(f1[4]*dt-f2[4]*double(k)*dt+f1[4]*double(k)*dt)/dt; 	
			SCV_inlet = 		(f2[5]-f1[5])/dt*t+(f1[5]*dt-f2[5]*double(k)*dt+f1[5]*double(k)*dt)/dt; 	
			ICV_inlet = 		(f2[6]-f1[6])/dt*t+(f1[6]*dt-f2[6]*double(k)*dt+f1[6]*double(k)*dt)/dt; 	
			PA_inlet = 			(f2[7]-f1[7])/dt*t+(f1[7]*dt-f2[7]*double(k)*dt+f1[7]*double(k)*dt)/dt; 				
			 			
        }                                                                                                                                                                                                                
    }                                                                                                                                                                                                                    

    // Log BC data.
    static double t_print = -1.0;
    if (fill_time > t_print && fill_time < 100.0) // XXXX: hack!
    {
        tbox::plog << "PVs_pressure_inlet  = " << RSPV_pressure_inlet * 0.00075006158 << " mmHg\n";
		tbox::plog << "Aorta_inlet  = " << Aorta_inlet * 0.00075006158 << " mmHg\n";
        t_print = fill_time;
    }

    // Patch geometry data.
    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const XLower = pgeom->getXLower();
    const double* const grid_XLower = d_grid_geometry->getXLower();
    const double* const grid_XUpper = d_grid_geometry->getXUpper();
    const double grid_XCenter[NDIM] = { 0.5*(grid_XUpper[0]+grid_XLower[0]) , 0.5*(grid_XUpper[1]+grid_XLower[1]) , 0.5*(grid_XUpper[2]+grid_XLower[2]) };
    const double* const dx = pgeom->getDx();
    double X[NDIM];

    // Boundary box data.
    const int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis =  location_index / 2;
    const hier::Box<NDIM>& bc_coef_box = acoef_data->getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
#endif

    // Prescribe boundary condition values.
    for (hier::Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const hier::Index<NDIM>& i = bc();
        for (int d = 0; d < NDIM; ++d)
        {
            if (d != bdry_normal_axis)
            {
                X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
            }
            else
            {
                X[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d)));
            }
        }

        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i,0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i,0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i,0) : dummy);

        if (d_component_normal_axis != bdry_normal_axis)  // tangential velocity boundary conditions
        {
            // All tangential velocities are clamped to zero.
            a = 1.0;
            b = 0.0;
            g = 0.0;
        }
        else                                              // normal     velocity boundary conditions
        {
            if      (location_index == 0) // RIPV [-4.8 -12.896 144.008] and RSPV [-4.8 -17.159 145.39]
            {
				const double r_RIPV = sqrt((X[1]+12.896)*(X[1]+12.896) + (X[2]-144.008)*(X[2]-144.008));
                const bool inside_RIPV = r_RIPV < 1.58;
				const double r_RSPV = sqrt((X[1]+17.159)*(X[1]+17.159) + (X[2]-145.39)*(X[2]-145.39));
                const bool inside_RSPV = r_RSPV < 1.58;
				if (inside_RIPV)
				{
					a = 0.0;
					b = 1.0; 
					g = -RIPV_pressure_inlet;
				}
				else if (inside_RSPV)
				{
					a = 0.0;
					b = 1.0; 
					g = -RSPV_pressure_inlet;
				}
				else
				{
					a = 0.0;
					b = 1.0; 
					g = 0.0;
				}

            }
			else if (location_index == 1)//None
			{
				a = 0.0;
                b = 1.0; 
				g = 0.0;	

			}
			else if (location_index == 2)//None
			{
				a = 0.0;
                b = 1.0; 
				g = 0.0;
			}
			else if (location_index == 3)//LIPV[5.782 -10.42 144.527]
			{
				const double r_LIPV = sqrt((X[0]-5.782)*(X[0]-5.782) + (X[2]-144.527)*(X[2]-144.527));
                const bool inside_LIPV = r_LIPV < 1.6;
				if (inside_LIPV)
				{
					a = 0.0;
					b = 1.0; 
					g = -LIPV_pressure_inlet;
				}
				else
				{
					a = 0.0;
					b = 1.0; 
					g = 0.0;
				}
			}
			else if (location_index == 4)
			{
					a = 0.0;
					b = 1.0; 
					g = 0.0;
			}
			else if (location_index == 5)//LSPV[6.6 -14.662 149.45] 1.6, SCV[-0.467 -17.055 149.45] 1.55, Aorta [-0.81 -21.32 149.45] 2.55, PAtrunk[3.529 -18.166 149.45] 2.1
			{
				const double r_LSPV = sqrt((X[0]-6.6)*(X[0]-6.6) + (X[1]+14.662)*(X[1]+14.662));
                const bool inside_LSPV = r_LSPV < 1.6;
				const double r_Aorta = sqrt((X[0]+0.81)*(X[0]+0.81) + (X[1]+21.32)*(X[1]+21.32));
                const bool inside_Aorta = r_Aorta < 2.58;
				if (inside_Aorta)
				{
					a = 0.0;
					b = 1.0; 
					g = -Aorta_inlet;
					if (Mytime >= 1.0)
					{
						g = -p_new[0];
					}
				}
				else if (inside_LSPV)
				{
					a = 0.0;
					b = 1.0; 
					g = -LSPV_pressure_inlet;
				}
				else
				{
					a = 0.0;
					b = 1.0; 
					g = 0.0;
				}
			}
        }
    }
    return;
}// setBcCoefs

hier::IntVector<NDIM>
VelocityBcCoefs::numberOfExtensionsFillable() const
{
    return 128;
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
