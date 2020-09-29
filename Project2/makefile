all:	compile	execute

compile:
	c++	-o	main.out	main.cpp Jacobi_rotation.cpp	-larmadillo
	c++ -o tests.out tests-main.cpp test-functions.cpp Jacobi_rotation.cpp -larmadillo

execute:
	./main.out
	./tests.out
