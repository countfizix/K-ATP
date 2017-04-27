all: model fixed

atp.o : atp.h atp.c
	gcc  -c atp.c

driver.o : atp.h driver.c atp.c
	gcc  -c driver.c 

out.o : atp.h out.c
	gcc  -c out.c

radau5.o: atp.h radau5.f
	gfortran -fopenmp -fno-align-commons -L/usr/local/g2c/lib  -c radau5.f

decsol.o: atp.h decsol.f
	gfortran  -fno-align-commons -L/usr/local/g2c/lib  -c decsol.f

model:  radau5.o decsol.o out.o driver.o atp.o atp.h cblock.h
	gfortran -fno-align-commons -W -Wall  -o model radau5.o decsol.o  out.o  driver.o -lm

fixed_finder.o: atp.h fixed_finder.c
	gcc -c fixed_finder.c -lm

fixed:  fixed_finder.o fixed_finder.c atp.h
	gcc -W -Wall -o stable fixed_finder.o -lm

clean:
	$(RM) *.o
