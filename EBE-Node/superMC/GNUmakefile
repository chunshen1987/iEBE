# ===========================================================================
#  Makefile superMC                                    Chun Shen Mar. 19, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := $(shell which icpc)
LD := $(shell which icpc)
CFLAGS= -O3 -fast -Wall $(shell gsl-config --cflags --libs)
ifeq "$(CC)" ""
  CC := $(shell which g++)
  LD := $(shell which g++)
  CFLAGS= -O3 -Wall $(shell gsl-config --cflags --libs)
endif

RM		=	rm -f
O               =       .o
LDFLAGS         =       -O3 -Wall $(shell gsl-config --cflags --libs)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	superMC.e
endif

SRC		=	main.cpp Bases.cpp MCnucl.cpp GlueDensity.cpp MakeDensity.cpp \
			KLNModel.cpp OverLap.cpp Largex.cpp Regge96.cpp rcBKfunc.cpp \
			MathBasics.cpp ParameterReader.cpp arsenal.cpp EOS.cpp \
			GaussianNucleonsCal.cpp NBD.cpp RandomVariable.cpp \
			TableFunction.cpp Table.cpp \

INC		= 	ArraySaver.h Bases.h CollisionPair.h EOS.h GaussianNucleonsCal.h \
			GlueDensity.h Integral.h KLNModel.h KLNfunc.h Largex.h \
			MCnucl.h MakeDenstiy.h MathBasics.h NBD.h Overlap.h \
			ParamDefs.h ParameterReader.h Participant.h Particle.h \
                  RandomVariable.h Regge96.h Stopwatch.h Table.h TableFunction.h \
                  UnintegPartonDist.h arsenal.h rcBKfunc.h \


# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
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
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------

