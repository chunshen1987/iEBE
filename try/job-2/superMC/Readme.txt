README for superMC --- features you should be aware of


---------------------------------------------------------------------
09-20-2011
---------------------------------------------------------------------
.: All results are outputted to the "data" subfolder; make sure this folder exists before executing the code.
.: All files under "data" subfolder need to be cleaned between runs since the code *append* to the existing files (which is a simpler choice to write data when more than one rapidity is used).
.: To test how fast the code is, use the "Stopwatch" class (see Stopwatch.h).
.: The parameter "auto_backup_after_number_of_averaging" can be used to auto-backup intermediate results during generating averaged profiles, and these intermediate results will be deleted once the program successfully finishes. If the program is terminated in the middle and re-run, it will continue from the last "snapshot". The auto-backup requires a subfolder called "backup" to be established. Also, if the code exits abnormally but still the backup is not favored, clean all files under "backup" subfolder.
