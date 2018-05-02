all: kruskal, boruvka

kruskal: kruskal-io.c kruskal.cpp
	mpic++ -Wall -O3 kruskal-io.c -lm -o kruskal

boruvka: boruvka-io.c boruvka.cpp
	mpic++ -Wall -O3 boruvka-io.c -lm -o boruvka

clean:
	rm kruskal
	rm boruvka
