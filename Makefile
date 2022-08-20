targets = pframe wtpframe quadratic testmpeigen testmp
FLAGS = -std=c++11 -O3 -Wpedantic -msse2
LIBS = -lmpfr -lgmp
INCLUDES = -I. -I/usr/include/eigen3 

wtpframe: wtpframe.cpp parseinput.cpp ./wtpframe.h ./solvers/LBFGS.h
	g++ $(FLAGS) $(INCLUDES) wtpframe.cpp parseinput.cpp -o wtpframe $(LIBS)

pframe: pframe.cpp ./include/pframe.h ./LBFGS.h ./descent.h
	g++ $(FLAGS) $(INCLUDES)  pframe.cpp -o pframe $(LIBS)

quadratic: example-quadratic.cpp ./solvers/descent.h
	g++ $(FLAGS) $(INCLUDES) example-quadratic.cpp -o quadratic $(LIBS)

testmpeigen: testmpeigen.cpp
	g++ $(FLAGS) $(INCLUDES) testmpeigen.cpp -o testmpeigen $(LIBS) 

testmp: testmp.cpp
	g++ $(FLAGS) testmp.cpp -o testmp $(LIBS)

getmpfr:
	wget "https://raw.githubusercontent.com/advanpix/mpreal/master/mpreal.h"
	mv mpreal.h include/

geteigen:
	wget -o mpfrcpp.zip "https://gitlab.com/libeigen/eigen/-/archive/3.2.10/eigen-3.2.10.zip" 
	unzip eigen-3.2.10.zip  
	rm eigen-3.2.10.zip
	mv eigen-3.2.10 include/eigen3

# sine: sine.cpp
# 	g++ -O2 -Wno-deprecated-declarations sine.cpp -o sine -lmpfr -lgmp 

.PHONY: clean
clean:
	-rm -f *.o 
	-rm $(targets)
