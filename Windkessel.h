#ifndef included_Windkessel
#define included_Windkessel
//#include "tools.h"

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>
// SAMRAI INCLUDES
#include <PatchHierarchy.h>

#include <vector>
//#include <blitz/tinyvec.h>
#include <boost/array.hpp>
//for linux system call
#include <stdlib.h>
#include <stdio.h>
// Windkessel is a static class that provides data and functions
// required to implement flow rate calculating.
class Windkessel
{
public:

	string d_object_name;

	bool d_registered_for_restart;

	/*!
     * \brief Windkessel model data.
     */
    double d_time, d_rsrc;
    std::vector<double> d_psrc, d_posn, d_P_PV_current,d_P_PV_predict;

	//double d_t_end_diastole;

	/*!
     * \brief The level of the patch hierarchy on which the Lagrangian
     * structures that interface the boundary are located.
     */
    //int d_bdry_interface_level_number;

	/*!
     * \brief Upstream pressure waveform data.
     */
    //vector<double> d_upstream_pres_waveform;
    //double d_upstream_pres_waveform_dt;


	/*!
     * \brief Constructor
     */
    Windkessel(
        const string& object_name,
        Pointer<Database> input_db,
		const tbox::Pointer<tbox::Database> input_db_bc,
        bool register_for_restart=true);

    /*!
     * \brief Destructor.
     */
    virtual
    ~Windkessel();

	void
    putToDatabase(
        Pointer<Database> db);


    void
    advanceTimeDependentData(
    std::vector<double> &Q_current,
	std::vector<double> &P_current, 
	const double d_time_input,
    const double dt,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int U_idx,
    const int P_idx,
    const int wgt_cc_idx,
    const int wgt_sc_idx,
	const int iteration_num
	);
	
	void 
	predictPressure(   
    std::vector<double> &Q_old,
	std::vector<double> &Q_current,
	std::vector<double> &dQdt,
	std::vector<double> &P_current,
	std::vector<double> &P_predict,	
	const double dt,	
	const int iteration_num);

private:

	Windkessel();

    Windkessel(Windkessel&);

    Windkessel& operator=(Windkessel&);
	

	void
    getFromRestart();
};

#endif
