ifeq ($(os),win)
	CC=x86_64-w64-mingw32-gcc
	OUT=cusum.exe
else
	CC=gcc
	OUT=cusum
endif
	
CFLAGS=-D_GNU_SOURCE -O3 -Wall -Wextra -lm --static
DEPS=bessel.h detector.h io.h stepfit.h lmmin_int64.h utils.h
ODIR=obj
_OBJ=main.o bessel.o detector.o io.o lmmin_int64.o stepfit.o utils.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
LIBS=-lm


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
	
$(OUT): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(OUT) $(ODIR)/*.o *~ core $(INCDIR)/*~ 
