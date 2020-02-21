targets = pframe wtpframe quadratic testmpeigen testmp
FLAGS = -std=c++11 -O3 -Wpedantic
LIBS = -lmpfr -lgmp
INCLUDES = -Iinclude -I/usr/include/eigen3

wtpframe: wtpframe.cpp ./include/wtpframe.h ./include/solvers/LBFGS.h
	g++ $(FLAGS) $(INCLUDES) wtpframe.cpp -o wtpframe $(LIBS)

pframe: pframe.cpp ./include/pframe.h ./include/LBFGS.h ./include/descent.h
	g++ $(FLAGS) $(INCLUDES)  pframe.cpp -o pframe $(LIBS)

quadratic: example-quadratic.cpp ./include/solvers/descent.h
	g++ $(FLAGS) $(INCLUDES) example-quadratic.cpp -o quadratic $(LIBS)

testmpeigen: testmpeigen.cpp
	g++ $(FLAGS) $(INCLUDES) testmpeigen.cpp -o testmpeigen $(LIBS) 

testmp: testmp.cpp
	g++ $(FLAGS) testmp.cpp -o testmp $(LIBS)

getmpfr:
	wget -o mpfrcpp.zip "http://www.holoborodko.com/pavel/downloads/mpfrc++-3.6.2.zip" 
	unzip -p mpfrc++-3.6.2.zip mpreal.h > mpreal.h
	mv mpreal.h include/

# sine: sine.cpp
# 	g++ -O2 -Wno-deprecated-declarations sine.cpp -o sine -lmpfr -lgmp 

.PHONY: clean
clean:
	-rm -f *.o 
	-rm $(targets)
