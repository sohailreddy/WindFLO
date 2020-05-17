CCOMPILER = g++

#PATH TO LIBGFORTRAN
LGFORTRAN = -L/usr/local/lib/gcc/9/


#PATH TO STATIC LIBRARY
LWINDFLO = ../../release/libWindFLO.a
#PATH TO API
API = ../../API/
#PATH TO OPTIMIZERS
OPTIMIZERS = ../../Optimizers/

EXENAME = example2

ifndef MAIN
	MAIN = main
endif
ifndef EXENAME
	EXENAME = $(MAIN)
endif


OBJECTS = $(MAIN).o GA.o

CFLAG += -std=c++11 -I$(API)

BUILD:
	make COMPILE;
	make LINKALL;

COMPILE: $(OBJECTS)

$(MAIN).o: $(MAIN).cpp
	$(CCOMPILER) $(CFLAG) -c $(MAIN).cpp	
GA.o: $(OPTIMIZERS)GA.cpp
	$(CCOMPILER) $(CFLAG) -c $(OPTIMIZERS)GA.cpp	


LINKALL:
	g++ $(OBJECTS) $(LWINDFLO) $(LGFORTRAN) -g -lgfortran -o $(EXENAME)
	           	    	           
clean :
	rm *.o $(EXENAME)
