#Copyright 2014 Oleguer Sagarra. All rights reserved. Code under License GPLv3.

CC   = clang                  # compiler, can be set to gcc also

# Compiler options:
INCLUDES = -I/opt/local/include/ -I/usr/local/include/ -I./include # include paths
CFLAGS   = -O2 -Wall $(INCLUDES)  #flags of the compiler

# Linker options:
LIBS     = -lm -lgsl -lgslcblas # libraries
LDFLAGS  = -L/opt/local/lib/ -L/usr/local/lib/ $(LIBS)  # library directories

# Files:
CFILES   = main.c ula_gen_funcs.c ula_stat_funcs.c ula_read_funcs.c \
		   ula_w_graph_funcs.c ula_w_graph_dist_funcs.c \
		   ula_null_models.c others_null_models.c     

# Actual files are in src directory:
SRCDIR = src
CSRCS    = $(addprefix $(SRCDIR)/, $(CFILES))

# Objects will be created in $(OBJDIR)/*.o
OBJDIR = src/obj
OBJECTS  = $(addprefix $(OBJDIR)/, $(CFILES:%.c=%.o))

all: $(OBJDIR) MultiEdgeGen

MultiEdgeGen: $(OBJDIR) $(OBJECTS)
	@echo Linking MultiEdgeGen
	@$(CC) $(OBJECTS) -o MultiEdgeGen $(LDFLAGS)


$(OBJDIR):
	@echo Create object directory: $@
	@mkdir -p $(OBJDIR)


$(OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.c
	@echo Compiling $<
	@$(CC) $(CFLAGS) -c $< -o $@

clean:
	@rm -rf $(OBJDIR) *~ src/*~ MultiEdgeGen

# By default make understands that each rule refers to a file/directory.
# By saying that "all" and "clean" rules are PHONY we are telling make
# that it should not expect a "./all" nor "./clean" file to be created
.PHONY: all clean

