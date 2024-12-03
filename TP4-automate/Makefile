SYSCONF_LINK = g++
CPPFLAGS     = -O3
CPPFLAGS     =
LDFLAGS      =
LIBS         =

SRCDIR = src
BINDIR = bin

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%.o,$(SOURCES))

EXEC = $(BINDIR)/postscript

all: $(EXEC)

$(EXEC): $(OBJECTS)
		$(SYSCONF_LINK) -g -o $@ $^ $(LDFLAGS)

$(BINDIR)/%.o: $(SRCDIR)/%.cpp
		@mkdir -p $(BINDIR)
		$(SYSCONF_LINK) -g -o $@ -c $< $(CFLAGS)

clean:
		rm -rf $(BINDIR)

mrproper: clean
		rm -rf $(EXEC)