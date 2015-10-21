prefix=/Users/holgerthies/iRRAM/installed
exec_prefix=/Users/holgerthies/iRRAM/installed

CC := clang++ # This is the main compiler
SRCDIR := src
BUILDDIR := build
TESTDIR := test
 
SRCEXT := cc
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g -Wall -std=c++11 -Xlinker -rpath -Xlinker /Users/holgerthies/iRRAM/installed/lib
LIB := -L/Users/holgerthies/iRRAM/installed/lib -liRRAM -lmpfr -lgmp -lm -lpthread

INC := -I include -I/Users/holgerthies/iRRAM/installed/include

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) $(LIB) -o $@ $<"; $(CC) $(CFLAGS) $(INC) $(LIB) -o $@ $<

$(TESTDIR)/%.o: $(TESTDIR)/%.$(SRCEXT)
	@mkdir -p $(TESTDIR)
	@echo " $(CC) $(CFLAGS) $(INC) $(LIB) -o $@ $< build/ANALYTIC.o"; $(CC) $(CFLAGS) $(INC) $(LIB) -o $@ $< src/ANALYTIC.cc src/combinatorics.cc

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

