test: testsrc/test.o cpp/conf.o testcpp/conf_test.o
	g++ -o $@ $^

%.o: %.cpp
	g++ -o $@ -c $<
