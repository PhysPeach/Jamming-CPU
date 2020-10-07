test: testsrc/test.o
	g++ -o $@ $^

%.o: %.cpp
	g++ -o $@ -c $<
