CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: simulator.exe
simulator.exe : system.o NSL_SIMULATOR.o particle.o random.o
	$(CC) system.o particle.o random.o NSL_SIMULATOR.o -o simulator.exe $(AFLAGS)
NSL_SIMULATOR.o : NSL_SIMULATOR.cpp system.h particle.h
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CFLAGS)
system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o simulator.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*

move1_e:
	rm -r ../../lezione_04/OUTPUT1/Equilibration
	cp -r ../OUTPUT ../../lezione_04/OUTPUT1/
	mv ../../lezione_04/OUTPUT1/OUTPUT ../../lezione_04/OUTPUT1/Equilibration

move1:
	rm -r ../../lezione_04/OUTPUT1
	cp -r ../OUTPUT ../../lezione_04
	mv ../../lezione_04/OUTPUT ../../lezione_04/OUTPUT1

move2_e:
	rm -r ../../lezione_04/OUTPUT2/Equilibration
	cp -r ../OUTPUT ../../lezione_04/OUTPUT2/
	mv ../../lezione_04/OUTPUT2/OUTPUT ../../lezione_04/OUTPUT2/Equilibration

move2:
	rm -r ../../lezione_04/OUTPUT2
	cp -r ../OUTPUT ../../lezione_04
	mv ../../lezione_04/OUTPUT ../../lezione_04/OUTPUT2

move3_e:
	rm -r ../../lezione_04/OUTPUT3/Equilibration
	cp -r ../OUTPUT ../../lezione_04/OUTPUT3/
	mv ../../lezione_04/OUTPUT3/OUTPUT ../../lezione_04/OUTPUT3/Equilibration

move3:
	rm -r ../../lezione_04/OUTPUT3
	cp -r ../OUTPUT ../../lezione_04
	mv ../../lezione_04/OUTPUT ../../lezione_04/OUTPUT3

move_lezione6:
	rm -r ../../lezione_06/OUTPUT
	cp -r ../OUTPUT ../../lezione_06

move_es7.2_solid:
	rm ../../lezione_07/esercizio7.2/potential_energy_solid.dat
	cp ../OUTPUT/potential_energy.dat ../../lezione_07/esercizio7.2
	mv ../../lezione_07/esercizio7.2/potential_energy.dat ../../lezione_07/esercizio7.2/potential_energy_solid.dat

move_es7.2_liquid:
	rm ../../lezione_07/esercizio7.2/potential_energy_liquid.dat
	cp ../OUTPUT/potential_energy.dat ../../lezione_07/esercizio7.2
	mv ../../lezione_07/esercizio7.2/potential_energy.dat ../../lezione_07/esercizio7.2/potential_energy_liquid.dat

move_es7.2_gas:
	rm ../../lezione_07/esercizio7.2/potential_energy_gas.dat
	cp ../OUTPUT/potential_energy.dat ../../lezione_07/esercizio7.2
	mv ../../lezione_07/esercizio7.2/potential_energy.dat ../../lezione_07/esercizio7.2/potential_energy_gas.dat

move_es7.4_solid_0:
	rm -r ../../lezione_07/esercizio7.4/solid/OUTPUT0
	cp -r ../OUTPUT/ ../../lezione_07/esercizio7.4/solid
	mv ../../lezione_07/esercizio7.4/solid/OUTPUT ../../lezione_07/esercizio7.4/solid/OUTPUT0

move_es7.4_solid_1:
	rm -r ../../lezione_07/esercizio7.4/solid/OUTPUT1
	cp -r ../OUTPUT/ ../../lezione_07/esercizio7.4/solid
	mv ../../lezione_07/esercizio7.4/solid/OUTPUT ../../lezione_07/esercizio7.4/solid/OUTPUT1

move_es7.4_liquid_0:
	rm -r ../../lezione_07/esercizio7.4/liquid/OUTPUT0
	cp -r ../OUTPUT/ ../../lezione_07/esercizio7.4/liquid
	mv ../../lezione_07/esercizio7.4/liquid/OUTPUT ../../lezione_07/esercizio7.4/liquid/OUTPUT0

move_es7.4_liquid_1:
	rm -r ../../lezione_07/esercizio7.4/liquid/OUTPUT1
	cp -r ../OUTPUT/ ../../lezione_07/esercizio7.4/liquid
	mv ../../lezione_07/esercizio7.4/liquid/OUTPUT ../../lezione_07/esercizio7.4/liquid/OUTPUT1

move_es7.4_gas_0:
	rm -r ../../lezione_07/esercizio7.4/gas/OUTPUT0
	cp -r ../OUTPUT/ ../../lezione_07/esercizio7.4/gas
	mv ../../lezione_07/esercizio7.4/gas/OUTPUT ../../lezione_07/esercizio7.4/gas/OUTPUT0

move_es7.4_gas_1:
	rm -r ../../lezione_07/esercizio7.4/gas/OUTPUT1
	cp -r ../OUTPUT/ ../../lezione_07/esercizio7.4/gas
	mv ../../lezione_07/esercizio7.4/gas/OUTPUT ../../lezione_07/esercizio7.4/gas/OUTPUT1