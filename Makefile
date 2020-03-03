O = o
E =
CC = gcc
OUT = cusum$E
CFLAGS = -D_GNU_SOURCE -O3 -Wall -Wextra --static
DEPS = bessel.h detector.h io.h stepfit.h lmmin_int64.h utils.h
ODIR = obj
_OBJ = main.$O bessel.$O detector.$O io.$O lmmin_int64.$O stepfit.$O utils.$O
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
LIBS = -lm

$(ODIR)/%.$O: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

$(OUT): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean win cleanwin

clean:
	rm -f $(OUT) $(ODIR)/*.$O *~ core $(INCDIR)/*~
	rm -f $(OUT).exe w$(ODIR)/*.$O *~ core $(INCDIR)/*~

win:
	$(MAKE) CC=x86_64-w64-mingw32-gcc E=.exe O=obj ODIR=wobj
	
