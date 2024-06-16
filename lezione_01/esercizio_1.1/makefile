
CC = g++

CFLAGS = -Wall -O3 --std=c++11

esercizio_1.1.exe : esercizio_1.1.o ../../rangen/random.o
	$(CC) esercizio_1.1.o ../../rangen/random.o -o esercizio_1.1.exe

esercizio_1.1.o : esercizio_1.1.cpp
	$(CC) -c esercizio_1.1.cpp -o esercizio_1.1.o $(CFLAGS)

../../rangen/random.o : ../../rangen/random.cpp
	$(CC) -c ../../rangen/random.cpp -o ../../rangen/random.o $(CFLAGS)

clean:
	rm -f esercizio_1.1.o ../../rangen/random.o esercizio_1.1.exe
