LDLIBS=-lm
all: 1628  undust3 test_uncorr test_shift

undust3: undust3.c pnm.c pnm.h
test_uncorr: test_uncorr.c pnm.c pnm.h
test_shift:  test_shift.c pnm.c pnm.h
