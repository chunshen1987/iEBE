all	: iS.e resonance.e iInteSp.e

iS.e	:
	(cd src; make; make install)
resonance.e	:
	(cd resonance8; make; make install)
iInteSp.e	:	
	(cd iInteSp; make; make install)

distclean	:
	(cd src; make distclean)
	(cd resonance8; make distclean)
	(cd iInteSp; make distclean)
	rm *.e
