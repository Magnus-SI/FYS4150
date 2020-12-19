all:	compile	execute

compile:
	c++	-o	main.out	main.cpp diffusion.cpp
	c++	-o	isomain.out	isomain.cpp	diffusion.cpp
	c++	-o	tests.out	tests-main.cpp	test-functions.cpp diffusion.cpp

execute:
	./tests.out
	./main.out
	./isomain.out
