CC = g++
CFLAGS = -O3 --std=c++11

all: esercizio_8.2.exe

esercizio_8.2.exe : esercizio_8.2.o ../system.o ../../rangen/random.o
	$(CC) esercizio_8.2.o ../system.o ../../rangen/random.o -o esercizio_8.2.exe

esercizio_8.2.o : esercizio_8.2.cpp ../system.h ../../rangen/random.h
	$(CC) -c esercizio_8.2.cpp -o esercizio_8.2.o $(CFLAGS)

system.o : ../system.cpp ../system.h ../../rangen/random.h
	$(CC) -c ../system.cpp -o ../system.o $(CFLAGS)

random.o : random.cpp ../../rangen/random.h
	$(CC) -c random.cpp -o ../../rangen/random.o $(CFLAGS)

clean:
	rm -f *.o esercizio_8.2.exe
