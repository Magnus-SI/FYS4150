all:	compile	execute

compile:
	c++	-o	main.out	main.cpp ising2D.cpp
	c++	-o	tests.out	tests-main.cpp	test-functions.cpp ising2D.cpp

execute:
	./tests.out
	./main.out
