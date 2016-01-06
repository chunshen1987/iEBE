# ===========================================================================
#  Makefile osc2u                                    Chun Shen Mar. 19, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile		[uqmd.f]
##			TYPE	=	operating system	['uname']
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  


# ----------------------------------------- 


FC := gfortran
FFLAGS= $(shell ../crank/configure_compiler_fortran_Flags.sh)

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(FFLAGS)
SYSTEMFILES     =       $(SRCGNU)



# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	osc2u.e
endif

SRC		=	osc2u.f pdg2ityp.f vni_procev.f blockres.f dectim.f \
			gnuranf.f

INC		= 	ucoms.f comres.f

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.f:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(FC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj*

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------

./osc2u.f:	ucoms.f 
./procev.f:	ucoms.f
./pdg2ityp.f:	comres.f
