superMC
=====
####Monte Carlo event generator of initial density distributions for relativistic heavy-ion collisions

This code can generate fluctuating initial conditions according to Monte-Carlo Glauber model and Monte-Carlo Kharzeev-Levin-Nardi (KLN) "gluon saturation" model.

The MCKLN model part of the code was forked from a code package mckt, which was first developed by Yasushi Nara. Option for rcBK unintegrated gluon distributions was added by Javier Albacete. The current version of mckt can be found from [http://faculty.baruch.cuny.edu/naturalscience/physics/dumitru/CGC_IC.html]

Some features of this package that you should be aware of:
* All results are outputted to the "data" subfolder; make sure this folder exists before executing the code.
* All files under "data" subfolder need to be cleaned between runs since the code *append* to the existing files (which is a simpler choice to write data when more than one rapidity is used).
* To test how fast the code is, use the "Stopwatch" class (see Stopwatch.h).
* The parameter "backup_number" can be used to auto-backup intermediate results during generating averaged profiles, and these intermediate results will be deleted once the program successfully finishes. If the program is terminated in the middle and re-run, it will continue from the last "snapshot". The auto-backup requires a subfolder called "backup" to be established. Also, if the code exits abnormally but still the backup is not favored, clean all files under "backup" subfolder.
