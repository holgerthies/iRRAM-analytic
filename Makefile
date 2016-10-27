prefix=/Users/holgerthies/iRRAM/installed
exec_prefix=/Users/holgerthies/iRRAM/installed

CC := clang++ # This is the main compiler
BUILDDIR := build
TESTDIR := test
 
SRCEXT := cc
CFLAGS := -g -Wall -std=c++14 -Xlinker -rpath -Xlinker /Users/holgerthies/iRRAM/installed/lib
LIB := -L/Users/holgerthies/iRRAM/installed/lib -liRRAM -lmpfr -lgmp -lm -lpthread

INC := -I include -I/Users/holgerthies/iRRAM/installed/include


$(TESTDIR)/%.o: $(TESTDIR)/%.$(SRCEXT)
	@mkdir -p $(TESTDIR)
	$(CC) $(CFLAGS) $(INC) $(LIB) -lboost_unit_test_framework -o $@ $< 

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

