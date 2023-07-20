CC=gcc
COPTS=-g -Wall -std=c99
ALL=LUtest

all: $(ALL)

JUNK= *.o *~ *.dSYM LUtest

clean:
	-rm -rf $(JUNK)

run: LUtest
	./LUtest

LUtest: LUtest.o LUfact.o
	$(CC) $(COPTS) $^ -lm -o $@

LUfact.o: LUfact.c LUfact.h
LUtest.o: LUtest.c LUfact.h

.c.o:
	$(CC) -c $(COPTS) $<