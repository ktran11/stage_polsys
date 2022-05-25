CFLAGS = -std=c11 -Wall -Wextra
CPPFLAGS = -I../include/
LDFLAGS =

all: malpol

malpol: malpol.o
	$(CC) $(CFLAGS) -o $@ $^

malpol.o: malpol.c ../include/malpol.o
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

clean:
	@rm -f *~ *.o malpol

.PHONY: all clean
