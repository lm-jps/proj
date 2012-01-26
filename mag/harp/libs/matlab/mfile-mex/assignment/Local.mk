
LIBSUSED += -lmatch
CMEXLDFLAGS += -L$(OUTDIR)/$(CSDIR)
MEX2C_LIBS += -lmatch
LDFLAGS += -L.

# CC=/usr/local/gcc43/bin/gcc

# removed upon "make clean"
CLEAN_FILES += wmatch.o glib.o libmatch.a

LIBMATCH	:= $(OUTDIR)/$(CSDIR)/libmatch.a

LIBMATCH_OBJ	:= $(addprefix $(OUTDIR)/$(CSDIR)/, wmatch.o glib.o)
$(LIBMATCH):	$(LIBMATCH_OBJ)
	ar crus $@ $^

$(LIBMATCH_OBJ) :	$(notdir $(LIBMATCH_OBJ:%.o=%.c))
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

assignment.$(MEXEXT): $(LIBMATCH)

assignment.cli: $(LIBMATCH)

$(OUTDIR)/$(CSDIR)/wmatch.o: match.defs wmatch.c pairs.c pointer.c term.c unpairs.c readgraph.c
	$(CMEX) -outdir $(OUTDIR)/$(CSDIR) -c wmatch.c

$(OUTDIR)/$(CSDIR)/glib.o: match.defs glib.c
	$(CMEX) -outdir $(OUTDIR)/$(CSDIR) -c glib.c

