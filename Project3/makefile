all:	compile	execute

compile:
	python3 init.py
	c++	-o	main.out	main.cpp solar_system.cpp
	c++ -o tests.out tests-main.cpp test-functions.cpp solar_system.cpp

execute:
	./main.out
	./tests.out
