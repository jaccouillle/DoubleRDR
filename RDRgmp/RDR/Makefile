CC = gcc
CFLAGS = -Wall -g 
LIB = -lgmp
OBJ = RDRgmp.o

all: essaiGMP.exe

essaiGMP.exe: essaiGMP.c $(OBJ)
	$(CC) $(CFLAGS) essaiGMP.c $(OBJ) -o essaiGMP.exe $(LIB)

RDRgmp.o: RDRgmp.c

clean:
	rm *.o
	rm *~

