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

#ifndef included_FeedbackForcer
#define included_FeedbackForcer

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>

// SAMRAI INCLUDES
#include <PatchHierarchy.h>

// NAMESPACE
#include <ibamr/app_namespaces.h>

// FeedbackForcer implements a simple feedback-forcing approach to imposing
// Dirichlet boundary conditions weakly at domain boundaries.
class FeedbackForcer
    : public CartGridFunction
{
public:
	static double MV_start_close_time;
    static double MV_end_close_time; 
    /*!
     * \brief Constructor
     */
    FeedbackForcer(
        const string& object_name,
        const Pointer<PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Destructor.
     */
    ~FeedbackForcer();

    /*!
     * Velocity patch data descriptor index.
     */
    int d_U_data_idx;

    /*!
     * Penalty parameter.
     */
    double d_kappa;

    /*!
     * \name Implementation of CartGridFunction interface.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    bool
    isTimeDependent() const { return true; }

    /*!
     * \brief Set data on the specified patch interior.
     */
    void
    setDataOnPatch(
        int data_idx,
        Pointer<hier::Variable<NDIM> > var,  // needed to get around namespace clash between SAMRAI and libMesh
        Pointer<Patch<NDIM> > patch,
        double data_time,
        bool initial_time=false,
        Pointer<PatchLevel<NDIM> > patch_level=Pointer<PatchLevel<NDIM> >(NULL));

    //\}

private:
    FeedbackForcer();
    FeedbackForcer(const FeedbackForcer&);
    FeedbackForcer& operator=(const FeedbackForcer&);

    /*
     * The object name.
     */
    string d_object_name;

    /*
     * The patch hierarchy object.
     */
    Pointer<PatchHierarchy<NDIM> > d_patch_hierarchy;
};

#endif //#ifndef included_FeedbackForcer
