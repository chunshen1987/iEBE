# ===========================================================================
#  Makefile freestreaming                             Chun Shen Mar. 19, 2014
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS= -O3 $(shell gsl-config --cflags)

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS) $(shell gsl-config --libs)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	lm.e
endif

SRC		=	LdMatching.cpp Freestreaming.cpp CellData.cpp lmtest.cpp \
                  EOS.cpp arsenal.cpp mistools.cpp gauss_quadrature.cpp \
                  ParameterReader.cpp

INC		= 	LdMatching.h CellData.h Freestreaming.h Table.h \
                  gauss_quadrature.h EOS.h LdMatching.h arsenal.h mistools.h \
                  ParameterReader.h


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
./LdMatching.cpp : LdMatching.h Freestreaming.h CellData.h EOS.h ParameterReader.h gauss_quadrature.h mistools.h arsenal.h 
./Freestreaming.cpp : Freestreaming.h gauss_quadrature.h mistools.h arsenal.h
./CellData.cpp : CellData.h
./lmtest.cpp : LdMatching.h Freestreaming.h gauss_quadrature.h ParameterReader.h
./EOS.cpp : EOS.h arsenal.h
./arsenal.cpp : arsenal.h
./mistools.cpp : mistools.h
./gauss_quadrature.cpp : gauss_quadrature.h
./ParameterReader.cpp : arsenal.h ParameterReader.h
