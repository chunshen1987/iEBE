# ===========================================================================
#  Makefile urQMD                                    Chun Shen Mar. 19, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile		[uqmd.f]
##			TYPE	=	operating system	['uname']
##                      LINK    =       linking additional libraries
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  


# ----------------------------------------- 

FC := gfortran
FFLAGS = -O3 -cpp

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(FFLAGS)
SYSTEMFILES     =       $(SRCGNU)


# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	urqmd.e
endif

SRC = 1fluid.f addpart.f angdis.f anndec.f bessel.f blockres.f boxprg.f \
      cascinit.f coload.f dectim.f delpart.f detbal.f dwidth.f error.f \
      getmass.f getspin.f hepnam.f hepchg.f hepcmp.f init.f input.f iso.f \
      ityp2pdg.f jdecay2.f make22.f numrec.f output.f paulibl.f proppot.f \
      pythia6409.f saveinfo.f scatter.f siglookup.f string.f tabinit.f \
      uhmerge.f upmerge.f whichres.f urqmd.f $(SYSTEMFILES)

SRCGNU		=	genranf.f


# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(SYSTEMFILES) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC) )))
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
		$(FC) $(OBJECTS) $(LIBS) $(LDFLAGS) -o $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj*

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)/$(basename $(TARGET))
		strip $(INSTPATH)/$(basename $(TARGET))

# --------------- Dependencies -------------------

#++ Dependencies of a file are all those files that are *included* by the
#++ corresponding source file.
#++ To reference the object files, please prepend "$(OBJDIR)/" to the
#++ source file name and replace the ".f" by ".o".
$(OBJDIR)/1fluid.o:    options.f defs.f
$(OBJDIR)/addpart.o:   coms.f freezeout.f newpart.f
$(OBJDIR)/angdis.o:    coms.f
$(OBJDIR)/anndec.o:    comres.f comwid.f options.f newpart.f
$(OBJDIR)/blockres.o:  comres.f
$(OBJDIR)/boxprg.o:    boxinc.f coms.f comres.f options.f
$(OBJDIR)/cascinit.o:  coms.f inputs.f options.f
$(OBJDIR)/coload.o:    boxinc.f colltab.f coms.f comres.f newpart.f options.f
$(OBJDIR)/dectim.o:    coms.f colltab.f options.f
$(OBJDIR)/delpart.o:   coms.f comres.f freezeout.f options.f newpart.f
$(OBJDIR)/detbal.o:    coms.f comres.f comwid.f newpart.f options.f
$(OBJDIR)/dwidth.o:    coms.f comres.f comwid.f options.f
$(OBJDIR)/getmass.o:   comnorm.f comres.f
$(OBJDIR)/getspin.o:   comres.f
$(OBJDIR)/init.o:      boxinc.f colltab.f coms.f comres.f freezeout.f inputs.f \
                       newpart.f options.f
$(OBJDIR)/input.o:     boxinc.f coms.f comres.f inputs.f options.f
$(OBJDIR)/ityp2pdg.o:  comres.f
$(OBJDIR)/iso.o:       comres.f newpart.f
$(OBJDIR)/jdecay2.o:   newpart.f
$(OBJDIR)/make22.o:    comres.f coms.f comwid.f newpart.f options.f
$(OBJDIR)/output.o:    boxinc.f coms.f comres.f freezeout.f inputs.f newpart.f \
                       options.f
$(OBJDIR)/paulibl.o:   coms.f comres.f options.f
$(OBJDIR)/proppot.o:   boxinc.f coms.f options.f
$(OBJDIR)/saveinfo.o:  coms.f
$(OBJDIR)/scatter.o:   boxinc.f coms.f comres.f freezeout.f newpart.f options.f
$(OBJDIR)/siglookup.o: comres.f
$(OBJDIR)/string.o:    comres.f comstr.f options.f
$(OBJDIR)/tabinit.o:   comres.f comwid.f
$(OBJDIR)/uhmerge.o:   boxinc.f colltab.f coms.f comres.f inputs.f itypdata.f \
                       freezeout.f newpart.f options.f defs.f
$(OBJDIR)/upmerge.o:   newpart.f options.f
$(OBJDIR)/whichres.o:  comres.f options.f
