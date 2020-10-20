test: testsrc/test.o cpp/MT.o cpp/conf.o testcpp/conf_test.o cpp/particles.o testcpp/particles_test.o cpp/cells.o testcpp/cells_test.o cpp/jamming.o testcpp/jamming_test.o
	g++ -o $@ $^

%.o: %.cpp
	g++ -o $@ -c $<
