.PHONY : all
all : bin
	g++ *.o -o a.out
	echo "-----------"
	./a.out
.PHONY : bin
bin:
	g++ -c runge.cpp
	g++ -c main.cpp
.PHONY : clean
clean:
	rm -f *.o *.out *.csv *.png