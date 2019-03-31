targets = pframe wtpframe quadratic testmpeigen testmp

pframe: pframe.cpp pframe.h ./include/LBFGS.h ./include/descent.h
	g++ -O3 -fopenmp -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations pframe.cpp -o pframe -lmpfr -lgmp

wtpframe: wtpframe.cpp wtpframe.h ./include/LBFGS.h
	g++ -O3 -fopenmp -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations wtpframe.cpp -o wtpframe -lmpfr -lgmp

quadratic: example-quadratic.cpp ./include/LBFGS.h
	g++ -O3 -fopenmp -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations example-quadratic.cpp -o quadratic -lmpfr -lgmp

testmpeigen: testmpeigen.cpp
	g++ -O2 -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations testmpeigen.cpp -o testmpeigen -lmpfr -lgmp 

testmp: testmp.cpp
	g++ -O2 -Wno-deprecated-declarations testmp.cpp -o testmp -lmpfr -lgmp 

# sine: sine.cpp
# 	g++ -O2 -Wno-deprecated-declarations sine.cpp -o sine -lmpfr -lgmp 

.PHONY: clean
clean:
	-rm -f *.o 
	-rm $(targets)
