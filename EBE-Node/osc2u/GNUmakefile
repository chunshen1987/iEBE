# ===========================================================================
#  Makefile ubap                                    Henning Weber 10.09.1996
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

ifeq "$(TYPE)" ""
   TYPE := $(shell uname)
endif


ifeq "$(TYPE)" "AIX" 

FC      	= 	xlf
LD              =       xlf
RM              =       rm
O		=	.o
#FFLAGS          =       -O  -qextname  -qcheck -qextchk -bloadmap:loadmap.out
#LDFLAGS         =       -O  -qextname  -qcheck -qextchk -bloadmap:loadmap.out
FFLAGS     	=       -O2 -qextname 
LDFLAGS    	=       -O2 -qextname
SYSTEMFILES	=	$(SRCAIX)

endif

ifeq "$(TYPE)" "Darwin"

FC              =       gfortran
LD		=	gfortran
RM		=	rm
O               =       .o
FFLAGS          =       -O -Wall -Wsurprising
LDFLAGS         =       -O -Wall -Wsurprising
SYSTEMFILES     =       $(SRCGNU)

endif

ifeq "$(TYPE)" "GNU"

FC              =       g77
LD		=	g77
RM              =       rm
O               =       .o
FFLAGS          =       -g  
LDFLAGS         =       -g
SYSTEMFILES     =       $(SRCGNU)

endif

ifeq "$(TYPE)" "Linux" 

FC              =       g77
LD              =       g77
RM              =       rm 
O               =       .o
FFLAGS          =       -O -Wall
LDFLAGS         =       -O -Wall
SYSTEMFILES     =       $(SRCGNU)

endif


ifeq "$(TYPE)" "HP-UX"

FC              =       f77
LD		=	f77
RM		=	rm
O               =       .o
FFLAGS          =       -w
LDFLAGS         =       -g
SYSTEMFILES     =       $(SRCGNU)

endif

ifeq "$(TYPE)" "IRIX64"

FC              =       f77
LD		=	f77
RM		=	rm
O               =       .o
FFLAGS          =       -g -n32 -trapuv -C 
LDFLAGS         =       -g -n32
SYSTEMFILES     =       $(SRCSGI)

endif

ifeq "$(TYPE)" "PURE"

FC              =       f77
LD              =       purify f77
RM              =       rm 
O               =       .o 
FFLAGS          =       -w 
LDFLAGS         =       -g 
SYSTEMFILES     =       $(SRCGNU)

endif

ifeq "$(TYPE)" "ALPHA"

FC              =       f77
LD              =       f77
RM              =       rm 
O               =       .o 
FFLAGS          =       -V -g  -C -align dcommons   -check overflow
LDFLAGS         =       -V -g  -C -align dcommons   -check overflow
SYSTEMFILES     =       $(SRCALPHA)

endif

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	osc2u
endif

SRC		=	osc2u.f pdg2ityp.f vni_procev.f blockres.f dectim.f \
			gnuranf.f

INC		= 	ucoms.f comres.f

# -------------------------------------------------

OBJDIR		=	obj_$(TYPE)
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN).$(TYPE)
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
		$(LD) $(LDFLAGS) $(OBJECTS) -o $(TARGET)
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS) loadmap.out

distclean:	
		-rm $(TARGET) loadmap.out
		-rm -r obj_*

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------

./osc2u.f:	ucoms.f 
./procev.f:	ucoms.f
./pdg2ityp.f:	comres.f
