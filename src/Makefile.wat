#
# Wmake File for gse2mseed - For Watcom's wmake
# Use 'wmake -f Makefile.wat'

.BEFORE
	@set INCLUDE=.;$(%watcom)\H;$(%watcom)\H\NT
	@set LIB=.;$(%watcom)\LIB386

cc     = wcc386
cflags = -zq
lflags = OPT quiet OPT map LIBRARY ..\libmseed\libmseed.lib
cvars  = $+$(cvars)$- -DWIN32

BIN = ..\gse2mseed.exe

INCS = -I..\libmseed

all: $(BIN)

$(BIN):	gse2mseed.obj cm6.obj
	wlink $(lflags) name $(BIN) file {gse2mseed.obj cm6.obj}

# Source dependencies:
gse2mseed.obj:	gse2mseed.c cm6.h
cm6.obj:	cm6.c cm6.h

# How to compile sources:
.c.obj:
	$(cc) $(cflags) $(cvars) $(INCS) $[@ -fo=$@

# Clean-up directives:
clean:	.SYMBOLIC
	del *.obj *.map $(BIN)
