# ===========================================================================
#  Makefile photon_emission                           Chun Shen May 5, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

#CC := $(shell ../crank/configure_compiler_C++.sh)
#CFLAGS= $(shell ../crank/configure_compiler_C++_Flags.sh)
CC := g++
CFLAGS = -g

RM		=	rm -f
O               =       .o
HDF5LD          = -I/usr/local/hdf5/include
HDF5FLAGS       = -L/usr/local/hdf5/lib /usr/local/hdf5/lib/libhdf5_hl_cpp.a /usr/local/hdf5/lib/libhdf5_cpp.a /usr/local/hdf5/lib/libhdf5_hl.a /usr/local/hdf5/lib/libhdf5.a -lz -lm
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	hydro_photonEmission.e
endif

SRC		=	main.cpp Hydroinfo_h5.cpp Arsenal.cpp ThermalPhoton.cpp \
                  Table2D.cpp PhotonEmission.cpp ParameterReader.cpp tensor_trans.cpp \
                  gauss_quadrature.cpp

INC		= 	Hydroinfo_h5.h Arsenal.h ThermalPhoton.h Table2D.h Stopwatch.h \
                  PhotonEmission.h ParameterReader.h tensor_trans.h gauss_quadrature.h


# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $(HDF5LD) -c $< -o $@

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
		$(CC) $(LDFLAGS) $(HDF5FLAGS) $(OBJECTS) -o $(TARGET)
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
main.cpp : Hydroinfo_h5.h Stopwatch.h Arsenal.h ParameterReader.h PhotonEmission.h gauss_quadrature.h
Hydroinfo_h5.cpp : Hydroinfo_h5.h
ParameterReader.cpp : ParameterReader.h Arsenal.h
PhotonEmission.cpp : PhotonEmission.h Hydroinfo_h5.h ThermalPhoton.h tensor_trans.h ParameterReader.h
tensor_trans.cpp : tensor_trans.h
ThermalPhoton.cpp : Table2D.h ThermalPhoton.h ParameterReader.h Arsenal.h gauss_quadrature.h
Table2D.cpp : Arsenal.h Table2D.h
gauss_quardrature.cpp : gauss_quadrature.h
Arsenal.cpp : Arsenal.h

