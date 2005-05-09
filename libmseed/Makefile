
# Standard compiler parameters
CFLAGS = -O 

# Options specific for GCC
GCC = gcc
GCCFLAGS = -O2 -Wall

LDFLAGS =
LDLIBS =

OBJS = gswap.o lookup.o genutils.o fileutils.o traceutils.o msrutils.o unpackdata.o unpack.o packdata.o pack.o

LIB_A = libmseed.a
LIB_SO = libmseed.so.1.0
LIB_SO_ALIAS = libmseed.so

all: static

static: $(LIB_A)

shared: $(LIB_SO)

$(LIB_A): $(OBJS)
	rm -f $(LIB_A)
	ar -csq $(LIB_A) $(OBJS)

$(LIB_SO): $(OBJS)
	rm -f $(LIB_SO) $(LIB_SO_ALIAS)
	$(CC) -shared -Wl,-soname -Wl,$(LIB_SO_ALIAS) -o $(LIB_SO) $(OBJS)

clean:
	rm -f $(OBJS) $(LIB_A) $(LIB_SO) $(LIB_SO_ALIAS)

cc:
	@$(MAKE) "CC=$(CC)" "CFLAGS=$(CFLAGS)"

gcc:
	@$(MAKE) "CC=$(GCC)" "CFLAGS=$(GCCFLAGS)"

gcc64:
	@$(MAKE) "CC=$(GCC)" "CFLAGS=-m64 $(GCCFLAGS)"

debug:
	$(MAKE) "CC=$(CC)" "CFLAGS=-g $(CFLAGS)"

gccdebug:
	$(MAKE) "CC=$(GCC)" "CFLAGS=-g $(GCCFLAGS)"

install:
	@echo
	@echo "No install target, copy the library and header yourself"
	@echo
