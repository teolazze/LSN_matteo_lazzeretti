CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo

all: esercizio_9.1.exe

esercizio_9.1.exe : esercizio_9.1.o system.o ../../rangen/random.o
	$(CC) esercizio_9.1.o system.o ../../rangen/random.o -o esercizio_9.1.exe $(AFLAGS)

esercizio_9.1.o : esercizio_9.1.cpp system.h ../../rangen/random.h
	$(CC) -c esercizio_9.1.cpp -o esercizio_9.1.o $(CFLAGS)

system.o : system.cpp system.h ../../rangen/random.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)

random.o : random.cpp ../../rangen/random.h
	$(CC) -c random.cpp -o ../../rangen/random.o $(CFLAGS)

clean:
	rm -f *.o esercizio_9.1.exe
