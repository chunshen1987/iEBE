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

CC := $(shell ../../crank/configure_compiler_C++.sh)
CFLAGS= $(shell ../../crank/configure_compiler_C++_Flags.sh) $(shell gsl-config --cflags)

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS) $(shell gsl-config --libs)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	superMC.e
endif

SRC		=	main.cpp Bases.cpp MCnucl.cpp GlueDensity.cpp MakeDensity.cpp \
			KLNModel.cpp OverLap.cpp Largex.cpp Regge96.cpp rcBKfunc.cpp \
			MathBasics.cpp ParameterReader.cpp arsenal.cpp EOS.cpp \
			GaussianNucleonsCal.cpp NBD.cpp RandomVariable.cpp \
			TableFunction.cpp Table.cpp HulthenFunc.cpp \
                  GaussianDistribution.cpp Particle.cpp

INC		= 	ArraySaver.h Bases.h CollisionPair.h EOS.h GaussianNucleonsCal.h \
			GlueDensity.h Integral.h KLNModel.h KLNfunc.h Largex.h \
			MCnucl.h MakeDenstiy.h MathBasics.h NBD.h OverLap.h \
			ParamDefs.h ParameterReader.h Participant.h Particle.h \
                  RandomVariable.h Regge96.h Stopwatch.h Table.h TableFunction.h \
                  UnintegPartonDist.h arsenal.h rcBKfunc.h Spectator.h \
                  HulthenFunc.h GaussianDistribution.h


# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	../

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
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
./EOS.cpp: arsenal.h
./GaussianNucleonsCal.cpp: arsenal.h OverLap.h ParameterReader.h
./GlueDensity.cpp: MCnucl.h Participant.h CollisionPair.h
./KLNModel.cpp: Regge96.h kkp.cpp Bases.h ParamDefs.h UnintegPartonDist.h OverLap.h
./Largex.cpp: KLNModel.h rcBKfunc.h
./MCnucl.cpp: Regge96.h arsenal.h GaussianNucleonsCal.h ParameterReader.h NBD.h Particle.h OverLap.h KLNModel.h Participant.h GlueDensity.h CollisionPair.h Largex.h HulthenFunc.h
./MakeDensity.cpp: Regge96.h KLNfunc.h rcBKfunc.h ParamDefs.h MCnucl.h EOS.h ParameterReader.h ArraySaver.h Stopwatch.h OverLap.h UnintegPartonDist.h
./NBD.cpp: arsenal.h RandomVariable.h
./OverLap.cpp: MathBasics.h
./ParameterReader.cpp: arsenal.h
./Participant.h: Particle.h
./RandomVariable.cpp: arsenal.h TableFunction.h
./Table.cpp: arsenal.h
./TableFunction.cpp: arsenal.h Table.h 
./main.cpp: MakeDensity.h ParamDefs.h ParameterReader.h
./rcBKfunc.cpp: ParamDefs.h UnintegPartonDist.h
./HulthenFunc.cpp: HulthenFunc.h RandomVariable.h
./GaussianDistribution.cpp : GaussianDistribution.h RandomVariable.h
./Particle.cpp: Particle.h GaussianDistribution.h
