CC = gcc
CXX = g++
GMPHOME = /opt/local

LDFLAGS =  -L$(GMPHOME)/lib -lmpfr -lgmp 
INCLIST = -I$(GMPHOME)/include

BIN = remez


CFLAGS=  -O3
CXXFLAGS=  -O3
DFLAGS +=  -DGMP

#  targets
all:$(BIN)

.SUFFIXES:
.SUFFIXES:  .o .C .c

CSRC :=$(wildcard *.c)
CCSRC :=$(wildcard *.C)

COBJ=$(CSRC:.c=.o)
CCOBJ=$(CCSRC:.C=.o)

OBJS_SRC = $(SOBJ) $(CCOBJ) $(COBJ)
OBJS := $(notdir $(OBJS_SRC))

$(BIN):  $(OBJS) 
	$(CXX) $(OBJS) $(LDFLAGS) -o $(BIN)

.c.o:
	$(CC) -o $@ $(CFLAGS) $(DFLAGS) -c $(INCLIST) $<
.C.o:
	$(CXX) -o $@ $(CXXFLAGS) $(DFLAGS) -c $(INCLIST) $<

clean:
	rm -f *.dat *.o  *~ *\#* *.stackdump $(BIN)
