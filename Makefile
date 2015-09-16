CC := g++ # This is the main compiler
SRCDIR := src
BUILDDIR := build
TESTDIR := test
 
SRCEXT := cc
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g -Wall -std=c++11
LIB := -pthread -lmongoclient -L lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt
INC := -I include


$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -o $@ $<"; $(CC) $(CFLAGS) $(INC) -o $@ $<

$(TESTDIR)/%.o: $(TESTDIR)/%.$(SRCEXT)
	@mkdir -p $(TESTDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -o $@ $<"; $(CC) $(CFLAGS) $(INC) -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

