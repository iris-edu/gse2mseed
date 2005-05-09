
# Standard compiler parameters
CFLAGS = -O -I../libmseed

# Options specific for GCC
GCC = gcc
GCCFLAGS = -g -O2 -Wall -I../libmseed

LDFLAGS = -L../libmseed
LDLIBS = -lmseed -lm

OBJS = cm6.o gse2mseed.o

all: gse2mseed

gse2mseed: $(OBJS)
	$(CC) $(CFLAGS) -o ../$@ $(OBJS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJS) ../gse2mseed

cc:
	@$(MAKE) "CC=$(CC)" "CFLAGS=$(CFLAGS)"

gcc:
	@$(MAKE) "CC=$(GCC)" "CFLAGS=$(GCCFLAGS)"

gcc64:
	@$(MAKE) "CC=$(GCC)" "CFLAGS=-m64 $(GCCFLAGS)"

debug:
	$(MAKE) "CFLAGS=-g $(CFLAGS)"

gccdebug:
	$(MAKE) "CC=$(GCC)" "CFLAGS=-g $(GCCFLAGS)"

install:
	@echo
	@echo "No install target, copy the executable(s) yourself"
	@echo
