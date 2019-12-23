FCOMPILER = gfortran-8
CCOMPILER = g++-8

EXENAME = WindFLO

OBJECTS = userDefined.o

LGFORTRAN = -L/usr/local/lib/gcc/4.8/

FFLAG += -x f95-cpp-input -g
CFLAG += -std=c++11

BUILD:
	make COMPILE;
	make LINKALL;


COMPILE: $(OBJECTS)

userDefined.o: userDefined.cpp
	$(CCOMPILER) $(CFLAG) -c userDefined.cpp	

LINKALL:
	g++ *.o WindFLO.a $(LGFORTRAN) -g -lgfortran -o $(EXENAME)
	           	    	           
clean :
	rm 	$(OBJECTS)
	
