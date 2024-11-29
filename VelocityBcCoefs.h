#ifndef included_VelocityBcCoefs
#define included_VelocityBcCoefs

// Filename: VelocityBcCoefs.h
// Last modified: <11.Jun.2009 15:55:30 griffith@box230.cims.nyu.edu>
// Created on 04 May 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include "Windkessel.h"
// C++ STDLIB
#include <vector>

// NAMESPACE
using namespace SAMRAI;
using namespace std;

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class VelocityBcCoefs is an implementation of the strategy class
 * solv::RobinBcCoefStrategy that is used to specify velocity boundary
 * conditions that are appropriate for the pressure boundary conditions defined
 * by class PressureBcCoefs.
 */
class VelocityBcCoefs
    : public solv::RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor
     */
	 static double P_dynamic; 
	 static std::vector<double> p_new;
    VelocityBcCoefs(
        const string& object_name,
        const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry,
        const int component_normal_axis,
		const Windkessel& Windkessel_model,
        const double& Heart_pres_data_delta_t,
        const vector<vector<double>>& Heart_pres_data,
        const tbox::Pointer<tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual
    ~VelocityBcCoefs();

    /*!
     * \name Implementation of solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     *
     * \note In the original solv::RobinBcCoefStrategy interface, it was assumed
     * that \f$ b = (1-a) \f$.  In the new interface, \f$a\f$ and \f$b\f$ are
     * independent.
     *
     * \see solv::RobinBcCoefStrategy::setBcCoefs()
     *
     * \param acoef_data  Boundary coefficient data.
     *        The array will have been defined to include index range
     *        for corresponding to the boundary box \a bdry_box and
     *        appropriate for the alignment of the given variable.  If
     *        this is a null pointer, then the calling function is not
     *        interested in a, and you can disregard it.
     * \param bcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it
     *        is to be filled with the b coefficient.
     * \param gcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it
     *        is to be filled with the g coefficient.
     * \param variable    Variable to set the coefficients for.
     *        If implemented for multiple variables, this parameter
     *        can be used to determine which variable's coefficients
     *        are being sought.
     * \param patch       Patch requiring bc coefficients.
     * \param bdry_box    Boundary box showing where on the boundary the coefficient data is needed.
     * \param fill_time   Solution time corresponding to filling, for use when coefficients are time-dependent.
     */
    virtual void
    setBcCoefs(
        tbox::Pointer<pdat::ArrayData<NDIM,double> >& acoef_data,
        tbox::Pointer<pdat::ArrayData<NDIM,double> >& bcoef_data,
        tbox::Pointer<pdat::ArrayData<NDIM,double> >& gcoef_data,
        const tbox::Pointer<hier::Variable<NDIM> >& variable,
        const hier::Patch<NDIM>& patch,
        const hier::BoundaryBox<NDIM>& bdry_box,
        double fill_time=0.0) const;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     *
     * The "extension" used here is the number of cells that a boundary box
     * extends past the patch in the direction parallel to the boundary.
     *
     * Note that the inability to fill the sufficient number of cells past the
     * edge or corner of the patch may preclude the child class from being used
     * in data refinement operations that require the extra data, such as linear
     * refinement.
     *
     * The boundary box that setBcCoefs() is required to fill should not extend
     * past the limits returned by this function.
     */
    virtual hier::IntVector<NDIM>
    numberOfExtensionsFillable() const;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VelocityBcCoefs();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VelocityBcCoefs(
        const VelocityBcCoefs& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VelocityBcCoefs&
    operator=(
        const VelocityBcCoefs& that);

    /*
     * The object name.
     */
    const string d_object_name;

    /*
     * The grid geometry object.
     */
    const tbox::Pointer<geom::CartesianGridGeometry<NDIM> > d_grid_geometry;

    /*
     * The velocity component to be set by this object.
     */
    const int d_component_normal_axis;

    /*
     * The radius of the tube.
     */
    const double d_radius;

    /*
     * List of inlet pressures and times.
     */
    const double d_init_time;
    const double d_Heart_pres_data_delta_t;
    const vector<vector<double>>& d_Heart_pres_data;
	const Windkessel& d_Windkessel_model;


    /*
     * The (constant) pressure at the tube outlet.
     */
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <VelocityBcCoefs.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VelocityBcCoefs
