# 
# Small makefile for Saisie
#

CCC=g++
CFLAGS=-Wall -I/urs/local/include
CFLAGS=
LDFLAGS= -lGL -lGLEW -lGLU -lglut  -larmadillo

DEBUBFLAG=-g

SRCDIR = src
BINDIR = bin

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%.o,$(SOURCES))

EXEC = $(BINDIR)/main

all: $(EXEC)

$(EXEC): $(OBJECTS)
		$(CCC) -g -o $@ $^ $(LDFLAGS)

$(BINDIR)/%.o: $(SRCDIR)/%.cpp
		@mkdir -p $(BINDIR)
		$(CCC) -g -o $@ -c $< $(CFLAGS)

clean:
		rm -rf $(BINDIR)

mrproper: clean
		rm -rf $(EXEC)
