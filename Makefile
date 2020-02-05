targets = pframe wtpframe quadratic testmpeigen testmp

wtpframe: wtpframe.cpp wtpframe.h ./include/LBFGS.h
	g++ -O3 -fopenmp -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations wtpframe.cpp -o wtpframe -lmpfr -lgmp

pframe: pframe.cpp pframe.h ./include/LBFGS.h ./include/descent.h
	g++ -O3 -fopenmp -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations pframe.cpp -o pframe -lmpfr -lgmp

quadratic: example-quadratic.cpp ./include/LBFGS.h
	g++ -O3 -fopenmp -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations example-quadratic.cpp -o quadratic -lmpfr -lgmp

testmpeigen: testmpeigen.cpp
	g++ -O2 -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations testmpeigen.cpp -o testmpeigen -lmpfr -lgmp 

testmp: testmp.cpp
	g++ -O2 -Wno-deprecated-declarations testmp.cpp -o testmp -lmpfr -lgmp 

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
