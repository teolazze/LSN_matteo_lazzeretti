
CC = g++

CFLAGS = -Wall -O3 --std=c++11

esercizio_1.3.exe : esercizio_1.3.o ../../rangen/random.o
	$(CC) esercizio_1.3.o ../../rangen/random.o -o esercizio_1.3.exe

esercizio_1.3.o : esercizio_1.3.cpp
	$(CC) -c esercizio_1.3.cpp -o esercizio_1.3.o $(CFLAGS)

../../rangen/random.o : ../../rangen/random.cpp
	$(CC) -c ../../rangen/random.cpp -o ../../rangen/random.o $(CFLAGS)

clean:
	rm -f esercizio_1.3.o ../../rangen/random.o esercizio_1.3.exe