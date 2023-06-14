all: build

build: main.o functions.o
	g++ main.o functions.o -o build.out

# # # # #
main.o: main.cpp
	g++ -c main.cpp

functions.o: functions.cpp
	g++ -c functions.cpp

##################################################
#------------------- TESTS ----------------------#
##################################################
tests: tests.o functions.o
	g++ tests.o functions.o -o tests.out

tests.o: tests/tests.cpp
	g++ -c tests/tests.cpp


##################################################
#------------------ UTILITY ---------------------#
##################################################
clean:
	rm -f *.o *.gch

purge:
	rm -f *.o *.gch *.out; rm -rf build/