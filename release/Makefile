FCOMPILER = gfortran
CCOMPILER = g++
LIB_COMPILER = ar



ifndef OS
	OS = OSX
endif

ifdef MAIN
	OBJECTS = $(MAIN).o
endif

EXENAME = WindFLO_$(OS)
LIBNAME = libWindFLO_$(OS)
API = ../API/



OBJECTS += userDefined.o

LGFORTRAN = -L/usr/local/lib/gcc/9/


ifeq ($(OS),OSX)
	LIBFLAG = -std=c++11 -shared -fPIC -O2 $(LGFORTRAN) -lgfortran -Wl,-all_load
endif		
ifeq ($(OS),LINUX)
	LIBFLAG = -std=c++11 -shared -fPIC -O2 -Wl,--whole-archive	  
endif

FFLAG += -x f95-cpp-input -g -fcheck=all -fbacktrace -fdec-math  -fPIC 
CFLAG += -fPIC -std=c++11 -I$(API)


OPER = $(OBJECTS)
OPER += LINKALL
OPER += CREATE_LIB


BUILD:
	make COMPILE;
ifdef MAIN	
	make LINKALL;
endif
	make CREATE_LIB;


COMPILE: $(OBJECTS)

userDefined.o: userDefined.cpp
	$(CCOMPILER) $(CFLAG) -c userDefined.cpp	
$(MAIN).o: $(MAIN).cpp
	$(CCOMPILER) $(CFLAG) -c $(MAIN).cpp	




###
#	LINK OBJECTS INTO EXECUTABLE
##           
LINKALL:	           
	g++ $(OBJECTS) $(LIBNAME).a $(LGFORTRAN) -g -lgfortran -o $(EXENAME)



###
#	CREATE LIBRARIES
##           
CREATE_LIB: 
	@$(LIB_COMPILER) rc $(LIBNAME).a userDefined.o

ifeq ($(OS),OSX)
	@$(CCOMPILER) -o $(LIBNAME).so $(CFLAG) $(LIBFLAG) $(LIBNAME).a
endif		

ifeq ($(OS),LINUX)
	@$(CCOMPILER) -o $(LIBNAME).so $(CFLAG) $(LIBFLAG) $(LIBNAME).a	$(LGFORTRAN) -lgfortran -Wl,--no-whole-archive
endif



clean :
	rm 	*.o $(LIBNAME).so
	