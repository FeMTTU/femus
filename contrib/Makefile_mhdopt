#############################################################
############## FEMUS MAIN STATIC MAKEFILE ###################
#############################################################

######## flags for External Libraries for this Application ########
fm-enable-mpi     = yes
fm-enable-hdf5    = yes
fm-enable-petsc   = yes
fm-enable-libmesh = no
###################################################################


###############################################################################
# main="basename $PWD" from shell

FM_TARGET 	   :=  ./main_$(FM_FEMUS_METHOD)


all:: $(FM_TARGET)

### library
include ../../contrib/Make.libfemus


#### FEMuS specific APP #######
FM_APPSRCFILES  := $(wildcard ./*$(FM_SRCEXT))
FM_APPOBJECTS   := $(patsubst %$(FM_SRCEXT), %_$(FM_OBJ_SUFFIX), $(FM_APPSRCFILES))


$(FM_TARGET) :  $(FM_APPOBJECTS)
	@$(MAKE) libfemus
	@echo "Linking the application objects of "$@" with libfemus..."
	@echo "========="
	@echo $(FEMUS_LIBS)  | tr " " "\n"
	@echo "========="
	@echo $(FEMUS_LDFLAGS)  | tr " " "\n"
	@echo "========="
	@echo $(FEMUS_CPPFLAGS)  | tr " " "\n"
	@echo "========="
	@echo $(libmesh_LIBS)  | tr " " "\n"
	@echo "========="
	$(FEMUS_CXX) $(FEMUS_CXXFLAGS) -o $@ $^  $(FEMUS_LIBS)



.PHONY: clean clobber distclean

# Useful rules.
# in this rule, clean both the App objects, and the associated App library
clean:
	rm -f $(FM_APPOBJECTS) *~

clobber:
	@$(MAKE) clean
	@rm -f $(FM_TARGET)

distclean:
	@$(MAKE) clobber
	@rm -f src/*.o *.o *.g.o *.pg.o



# include the dependency list
include .depend


# Dependencies

.depend:
	@$(FM_PERL) $(FEMUS_DIR)/contrib/make_dependencies.pl -I. \
         $(FEMUS_INCLUDE) $(foreach i, $(wildcard $(FEMUS_DIR)/include/*), -I$(i)) "-S\$$(FM_OBJ_SUFFIX)"\
         $(FM_APPSRCFILES) > .depend
	@echo "Updated .depend"

###############################################################################
