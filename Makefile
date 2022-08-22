targets = pframe wtpframe quadratic testmpeigen testmp
FLAGS = -std=c++11 -O3 -Wpedantic -msse2
LIBS = -lmpfr -lgmp
INCLUDES = -Iinclude -Iinclude/eigen3 -I.

build/wtpframe: wtpframe.cpp parseinput.cpp ./wtpframe.h solvers/LBFGS.h
	if [ ! -d build ]; then mkdir build; fi
	g++ $(FLAGS) $(INCLUDES) wtpframe.cpp parseinput.cpp -o build/wtpframe $(LIBS)

pframe: pframe.cpp ./pframe.h solvers/LBFGS.h solvers/descent.h
	g++ $(FLAGS) $(INCLUDES)  pframe.cpp -o build/pframe $(LIBS)

quadratic: example-quadratic.cpp solvers/descent.h
	g++ $(FLAGS) $(INCLUDES) example-quadratic.cpp -o build/quadratic $(LIBS)

testmpeigen: testmpeigen.cpp
	g++ $(FLAGS) $(INCLUDES) testmpeigen.cpp -o build/testmpeigen $(LIBS) 

testmp: testmp.cpp
	g++ $(FLAGS) testmp.cpp -o build/testmp $(LIBS)

include/mpreal.h:
	if [ ! -d include ]; then mkdir include; fi
	wget "https://raw.githubusercontent.com/advanpix/mpreal/master/mpreal.h"
	mv mpreal.h include/

include/eigen3:
	wget "https://gitlab.com/libeigen/eigen/-/archive/3.2.10/eigen-3.2.10.zip" 
	unzip eigen-3.2.10.zip  
	rm eigen-3.2.10.zip
	if [ ! -d include/eigen3 ]; then mkdir -p include/eigen3; fi
	mv eigen-3.2.10 eigen3
	mv eigen3 include/

.PHONY: clean all
clean:
	-rm -rf build include

all: include/eigen3 include/mpreal.h build/wtpframe
