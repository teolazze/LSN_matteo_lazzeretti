CC = g++

CFLAGS = -Wall -O3 --std=c++11

esercizio_2.1.exe : esercizio_2.1.o ../../rangen/random.o
	$(CC) esercizio_2.1.o ../../rangen/random.o -o esercizio_2.1.exe

esercizio_2.1.o : esercizio_2.1.cpp
	$(CC) -c esercizio_2.1.cpp -o esercizio_2.1.o $(CFLAGS)

../../rangen/random.o : ../../rangen/random.cpp
	$(CC) -c ../../rangen/random.cpp -o ../../rangen/random.o $(CFLAGS)

clean:
	rm -f esercizio_2.1.o ../../rangen/random.o esercizio_2.1.exe