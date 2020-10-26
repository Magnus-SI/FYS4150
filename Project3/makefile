all:	compile	execute

compile:
	python3 make_initial.py
	c++ -o tests.out tests-main.cpp test-functions.cpp solar_system.cpp
	c++	-o	main.out	main.cpp solar_system.cpp


execute:
	./tests.out
	./main.out
