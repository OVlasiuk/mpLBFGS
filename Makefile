targets = pframe testmpeigen

pframe: pframe.cpp
	g++ -O3 -fopenmp -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations pframe.cpp -o pframe -lmpfr -lgmp

testmpeigen: testmpeigen.cpp
	g++ -O2 -Iinclude -I/usr/include/eigen3 -Wno-deprecated-declarations testmpeigen.cpp -o testmpeigen -lmpfr -lgmp 

testmp: testmp.cpp
	g++ -O2 -Wno-deprecated-declarations -Iinclude testmp.cpp -o testmp -lmpfr -lgmp 

.PHONY: clean
clean:
	-rm -f *.o 
	-rm $(targets)
