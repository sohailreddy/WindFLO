FCOMPILER = gfortran-10
CCOMPILER = g++-10
LGFORTRAN = #-L/usr/local/lib/gcc/10/ -I/usr/local/lib/gcc/10/


LIB_COMPILER = ar


ifndef OS
	OS = OSX
endif


EXENAME = WindFLO
LIBNAME = libWindFLO
RELEASEDIR = ../release/


OBJECTS = octree.o terrain.o WindFLO_f.o WindFLO.o wakeModels.o sobol.o ambientWindModels.o turbine.o \
			costModels.o convexHull.o userDefined.o CubicSpline.o API.o main.o

LIBOBJECTS := $(filter-out main.o,$(OBJECTS))





ifeq ($(OS),OSX)
	LIBFLAG = -std=c++11 -shared -fPIC -O2 $(LGFORTRAN) -lgfortran -Wl,-all_load
endif		
ifeq ($(OS),LINUX)
	LIBFLAG = -std=c++11 -shared -fPIC -O2 -Wl,--whole-archive	  
endif
ifeq ($(OS),WINDOWS)
	LIBFLAG = -std=c++11 -shared -fPIC -O2 $(LGFORTRAN) -lgfortran -Wl,-all_load
endif		

FFLAG = -x f95-cpp-input -g -fcheck=all -fbacktrace -fdec-math -fPIC 
CFLAG = -fPIC -std=c++11 -I./Include/


OPER = $(OBJECTS)
OPER += LINKALL
OPER += CREATE_LIB
#OPER += clean


BUILD:
	make COMPILE;
	make LINKALL;
	make CREATE_LIB;


COMPILE: $(OBJECTS)

octree.o: octree.f90
	$(FCOMPILER) $(FFLAG) -lstdc++ -c octree.f90
terrain.o: terrain.f90
	$(FCOMPILER) $(FFLAG) -lstdc++ -c terrain.f90	
WindFLO_f.o: WindFLO_f.f90
	$(FCOMPILER) $(FFLAG) -lstdc++ -c WindFLO_f.f90
WindFLO.o: WindFLO.cpp
	$(CCOMPILER) $(CFLAG) -c WindFLO.cpp
turbine.o: turbine.cpp
	$(CCOMPILER) $(CFLAG) -c turbine.cpp		
ambientWindModels.o: ambientWindModels.cpp
	$(CCOMPILER) $(CFLAG) -c ambientWindModels.cpp	
wakeModels.o: wakeModels.cpp
	$(CCOMPILER) $(CFLAG) -c wakeModels.cpp
userDefined.o: userDefined.cpp
	$(CCOMPILER) $(CFLAG) -c userDefined.cpp	
costModels.o: costModels.cpp
	$(CCOMPILER) $(CFLAG) -c costModels.cpp	
convexHull.o: convexHull.cpp
	$(CCOMPILER) $(CFLAG) -c convexHull.cpp	
sobol.o: sobol.cpp
	$(CCOMPILER) $(CFLAG) -c sobol.cpp	
API.o: API.cpp
	$(CCOMPILER) $(CFLAG) -c API.cpp	
CubicSpline.o: CubicSpline.cpp
	$(CCOMPILER) $(CFLAG) -c CubicSpline.cpp	
main.o: main.cpp
	$(CCOMPILER) $(CFLAG) -c main.cpp



###
#	LINK OBJECTS INTO EXECUTABLE
##           
LINKALL:	           
ifeq ($(OS),OSX)
	$(CCOMPILER) $(OBJECTS) $(LGFORTRAN) -g -lgfortran -o ../release/$(EXENAME)
endif	

ifeq ($(OS),LINUX)
	$(CCOMPILER) $(OBJECTS) $(LGFORTRAN) -g -lgfortran -o ../release/$(EXENAME)
endif		

ifeq ($(OS),WINDOWS)
	$(CCOMPILER) $(OBJECTS) $(LGFORTRAN) -g -lgfortran -o ../release/$(EXENAME)
endif		


###
#	CREATE LIBRARIES
##           
CREATE_LIB: 
	@$(LIB_COMPILER) rc $(LIBNAME).a $(LIBOBJECTS); cp $(LIBNAME).a $(RELEASEDIR)
	@$(CCOMPILER) -o $(LIBNAME).so $(CFLAG) $(LIBFLAG) $(LIBNAME).a
	cp $(LIBNAME).so $(RELEASEDIR)
	cp $(LIBNAME).so $(RELEASEDIR)../API/





#	@$(LIB_COMPILER) rc $(LIBNAME).a $(LIBOBJECTS)
#	@$(CCOMPILER) -o $(LIBNAME).so $(CFLAG) $(LIBFLAG) $(LIBNAME).a	
#	@$(CCOMPILER) -o $(LIBNAME).so -std=c++11 $(LGFORTRAN) -lgfortran -shared -fPIC -O2 $(CFLAG) userDefined.o -Wl,-all_load  $(LIBNAME).a
#	@$(CCOMPILER) -o $(LIBNAME).so -std=c++11 -fPIC -O2 $(CFLAG) -Wl,--whole-archive $(LIBNAME).a -Wl,--no-whole-archive userDefined.o $(LGFORTRAN) -lgfortran -shared			           
#	@$(CCOMPILER) -std=c++11 $(LGFORTRAN) -lgfortran -shared -fPIC -O2 $(CFLAG) -Wl,-all_load  $(LIBNAME).a userDefined.o -o $(LIBNAME).so 			
           
clean :
	rm 	$(OBJECTS) $(LIBNAME).* *.mod
	