default: prog

SRC = RKVektor.f95

%.f95: 
	f95 -o RKVektor $<

run:
	./RKVektor > werteliste.dat

prog: 
	f95 RKVektor.f95 -o RKVektor

plot: run
	gnuplot -e "set terminal png; plot 'werteliste.dat' using 1:2 with lines" > graph.png

view: plot
	gwenview graph.png