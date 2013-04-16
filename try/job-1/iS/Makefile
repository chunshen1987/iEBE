all	: iS.e resonance.e iInteSp.e

iS.e	:
	make -f GNUmakefile_iS
resonance.e	:
	(cd resonance8; make; make install)
iInteSp.e	:	
	(cd iInteSp; make; make install)

distclean	:
	make -f GNUmakefile_iS distclean
	(cd resonance8; make distclean)
	(cd iInteSp; make distclean)
	rm *.e
