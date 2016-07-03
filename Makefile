all: undust3 1628

undust3: undust3.cpp utils.o pnm.o
1628: 1628.c
test_pnm: test_pnm.cpp pnm.o
pnm.o:   pnm.cpp pnm.h
utils.o: utils.cpp utils.h