
PACK_DIR := "$(FEMUS_DIR)/packs"

all: 
	@echo "******* ---> Go to the ./applications folder of FEMuS library <--- *********"

myclean:
	rm -r -f *~
	rm -f *.o src/*.o 
	rm -r output/*
	rm -r input/*
	rm -f  output/* data_in/*.in data_in/*.h5 data_in/*.xmf
	rm -f contrib/gencase/*.o
	rm -f applications/gencase/*.o
	rm -f applications/main/*.o
	rm -f applications/gencase/gencase-opt applications/gencase/gencase-dbg 
	rm -f   applications/main/main   applications/main/main-opt   applications/main/main-dbg
	rm -f gencase gencase-opt gencase-dbg

pack:
	rm -r -f *~ 
	rm -f *.o src/*.o
	tar czvf $(PACK_DIR)/femus.tar.gz ../femus --exclude=data_in --exclude=output --exclude=packs \
	--exclude=contrib/hdf5-1.8.5-patch1-linux-x86_64-shared   --exclude=contrib/ioformat --exclude=contrib/vof  --exclude=contrib/oldies --exclude=fem --exclude=lib --exclude=.git --exclude=*/*/*.o --exclude=*/*/*/*kdev4 --exclude=*/*/*~ --exclude=*/*/*.old \
         --exclude=*.kdevelop*  --exclude=tests  --exclude=.git

packall:
	rm -r -f *~ 
	rm -f *.o src/*.o
	tar czvf  $(PACK_DIR)/femusall.tar.gz ../femus --exclude=data_in --exclude=output --exclude=packs


include Make.libfemus
