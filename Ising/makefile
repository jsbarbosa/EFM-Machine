FILE = ising_generator.c

animate : animation.py results/*.dat
	python3 animation.py

run : results.dat
	python3 plot.py

results.dat : a.out
	./a.out

a.out : $(FILE) headers.h
	gcc -O2 $(FILE) -lm -Wall -g

clean :
	rm a.out
