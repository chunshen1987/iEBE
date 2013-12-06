# ===========================================================================
#  Makefile VISHNew-Jet                               Chun Shen Apr. 23, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

FC := h5fc
FFLAGS= -O3 -cpp -lm -lz -fno-align-commons

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(FFLAGS) 
LIBSHDF         =       -L/usr/local/hdf5/lib /usr/local/hdf5/lib/libhdf5hl_fortran.a /usr/local/hdf5/lib/libhdf5_hl.a /usr/local/hdf5/lib/libhdf5_fortran.a /usr/local/hdf5/lib/libhdf5.a 
LIBSHDF         =  
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	VISHNew.e
endif

SRC		=	VISH2p1V1.10.0.for PhyBdary-1.10.for InputEOS-1.3.for \
                  OSCARoutput.for Arsenal-0.8.for Initialization-1.03.for \
                  InputFun-1.29RC6.for Jetoutputh5.for cornelius2.f90

INC		= 	

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.for
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

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
		$(FC) $(LDFLAGS) $(LIBSHDF) $(OBJECTS) -o $(TARGET)
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------

