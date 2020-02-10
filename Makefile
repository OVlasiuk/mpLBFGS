targets = pframe wtpframe quadratic testmpeigen testmp
FLAGS = -std=c++11 -O3 -Wno-deprecated-declarations
LIBS = -lmpfr -lgmp

wtpframe: wtpframe.cpp wtpframe.h ./include/LBFGS.h
	g++ $(FLAGS) -Iinclude -I/usr/include/eigen3 wtpframe.cpp -o wtpframe $(LIBS)

pframe: pframe.cpp pframe.h ./include/LBFGS.h ./include/descent.h
	g++ $(FLAGS) -Iinclude -I/usr/include/eigen3  pframe.cpp -o pframe $(LIBS)

quadratic: example-quadratic.cpp ./include/LBFGS.h
	g++ $(FLAGS) -Iinclude -I/usr/include/eigen3 example-quadratic.cpp -o quadratic $(LIBS)

testmpeigen: testmpeigen.cpp
	g++ $(FLAGS) -Iinclude -I/usr/include/eigen3 testmpeigen.cpp -o testmpeigen $(LIBS) 

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
