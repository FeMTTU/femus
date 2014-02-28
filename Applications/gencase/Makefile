############## GENCASE MAKEFILE ########################

########### Consistency of boundary conditions #######
# What if the physics tells you that the boundary conditions
# for a model are something,
# and then after you do the weak formulation it comes out 
# that those boundary conditions cannot be retrieved 
# in the boundary integrals?
# it means that those boundary conditions are MATHEMATICALLY WRONG...
# What if a bc is MATHEMATICALLY WRONG but PHYSICALLY CORRECT?

# Is the weak formulation the best criterion for checking
# the correctness of given boundary conditions from physics?
# Well, that is the point where integration by parts HIGHLIGHTS BOUNDARY INTEGRALS...

# Well, if you want to check the correctness 
# the ultimate criterion is VALIDATION with experiments...

##################
# When you have a PDE, the fact that it is LINEAR or NONLINEAR 
# is something that comes out from the STRONG formulation,
# the WEAK Formulation in PDE world only comes LATER,
# if you have to prove stability results,
# or do finite element discretization,
# or say something about the space of solutions...
# (what if the solution of a PDE, be it linear or nonlinear,
# does not belong to a VECTOR SPACE?
# then basically you could not do anything, because a linear combination
# would not be in the set...
# So, even if the equation is nonlinear, the space is a LINEAR SPACE... at least... :)
# )




######## moving to the new makefile
# Allora, prima c'era il make common di libmesh che ci pensava lui
# Ora il make common di libmesh e' stato spostato da qualche altra parte
# pero' se includo quello nuovo
#sembra non funzionare.
# Sembra che il nuovo sistema di libmesh non funzioni con il make common,
# ma in qualche modo invece dovrebbe essere cosi'...

############# CONFIGURE PHILOSOPHY ##########
# the philosophy is that you HARD CODE the system in TEXT FILES,
# so that NEXT TIME you dont have to set anything in the shell environment...
# In realta' alcune shell environment di LIBRERIE ESTERNE mi sa che devi tenertele...
# E' possibile che funzioni senza PETSC_... e senza MPI_DIR?
# (Chiaramente PATH c'e' sempre, e LD_LIBRARY_PATH)
# mi sa che di questi devo avere cura...

# L'importante e' che SE CAMBI QUALSIASI COSA nel SISTEMA
# (compilatore, librerie esterne)
# devi ricompilare le libmesh.......
# 
# Io avevo fatto in modo che il CAMBIO di SISTEMA
# sia SENZA SCRITTURA di FILE, totalmente indipendente





######## flags for external libraries ########
# every application must decide separately for the libraries it needs
# but gencase is a SPECIAL Application: it is the only one that 
# DOES NOT RECOMPILE the LIBRARY, but it uses the 
# SAME LIBRARY as the CORRESPONDING APPLICATION!
# So, like i pick the same configurations as 'fm_myapp',
# I should also pick the same MAKEFILE CONFIGURATIONS as 'fm_myapp'...
# While this is possible with the conf files of the application,
# it is not so easy for the makefile also

# Concerning gencase, we have to say this: we do not want to 
# RECOMPILE the LIB because of the gencase!
# We want to compile the lib for our App, and then gencase
# uses the lib AS COMPILED BY OUR APP, without modifying it!
# so, the fact of enabling libmesh only for gencase DOESNT MEAN
# that you have to recompile ALL of LIBFEMUS !

# in order to do this, you have to compile
# FIRST YOUR APP,
# then GENCASE!!!
# otherwise you may risk recompiling the library, which is what i dont want!
# in order to do that, i might say that gencase does not depend on the library,
# but only on its objects

fm-enable-mpi     = yes
fm-enable-hdf5    = yes
fm-enable-petsc   = yes
fm-enable-libmesh = yes
###############################################


###############################################################################
FM_TARGET 	   := ./gencase_$(FM_FEMUS_METHOD)

all:: $(FM_TARGET)

##############################The library
include ../../contrib/Make.libfemus

#### FEMuS specific APP #######

FM_APPSRCFILES  := $(wildcard ./*$(FM_SRCEXT))
FM_APPOBJECTS   := $(patsubst %$(FM_SRCEXT), %_$(FM_OBJ_SUFFIX), $(FM_APPSRCFILES))



$(FM_TARGET) : $(FM_APPOBJECTS)
	@$(MAKE) libfemus
	@echo "Linking the application objects of "$@" with libfemus..."
	$(FEMUS_CXX) $(FEMUS_CXXFLAGS) -o $@  $^ $(FEMUS_LIBS) 



.PHONY: clean clobber distclean



clean:
	rm -f $(FM_APPOBJECTS) *~

clobber:
	@$(MAKE) clean
	@rm -f $(FM_TARGET) out_*

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o


# include the dependency list
include .depend


# Dependencies

.depend:
	@$(FM_PERL) $(FEMUS_DIR)/contrib/make_dependencies.pl \
         $(FEMUS_INCLUDE) "-S\$$(FM_OBJ_SUFFIX)" $(FM_APPSRCFILES) > .depend
	@echo "Updated .depend"