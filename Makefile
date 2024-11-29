######################################################################
IBAMR_SRC_DIR=/xlwork3/2390378d/IBAMR/IBAMR
IBAMR_BUILD_DIR=/xlwork3/2390378d/IBAMR/ibamr-objs-opt
######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc
######################################################################
## Build the application.
SOURCES = Windkessel.C       Windkessel.h       \
FeedbackForcer.C       FeedbackForcer.h       \
		  VelocityBcCoefs.C      VelocityBcCoefs.h      \
          BoundaryConditions.C   BoundaryConditions.h   \
          MechanicsModel.C       MechanicsModel.h       \
          ModelInitialization.C  ModelInitialization.h  \
		  ActiveContraction.C  ActiveContraction.h  \
INSStaggeredHierarchyIntegrator.C  INSStaggeredHierarchyIntegrator.h  \
INSHierarchyIntegrator.C  INSHierarchyIntegrator.h  \
          main.C                                       
OBJS    = main.o FeedbackForcer.o VelocityBcCoefs.o BoundaryConditions.o MechanicsModel.o ModelInitialization.o ActiveContraction.o Windkessel.o INSHierarchyIntegrator.o INSStaggeredHierarchyIntegrator.o
PDIM    = 3

main3d: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(IBAMR_LIB_3D) $(IBTK_LIB_3D)  $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) -o main3d

clean:
	$(RM) main3d
	$(RM) *.o *.lo *.objs *.ii *.int.c
	$(RM) -r .libs

#depend:
#	makedepend -s "# THE FOLLOWING IS AUTOMATICALLY GENERATED BY MAKEDEPEND" \
#	-- -DNDIM=$(PDIM) $(CPPFLAGS) -- $(SOURCES)

######################################################################
# THE FOLLOWING IS AUTOMATICALLY GENERATED BY MAKEDEPEND