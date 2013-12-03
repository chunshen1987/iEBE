#! /usr/bin/env bash
#  This small script is used to add additional particle information into the database after the full ebe calculations are done.
#  In order to use this script, one needs to create a /workplace folder in iebe/ directory and copy iSS and EbeCollector into 
#  this folder. In the iSS, user should modify EOS/chosen_particles.dat to choose a list of the particles that he is interested
#  in and want to be collected to the final database. Users should also modify EbeCollector accordingly to collect additional 
#  species of particles. In the end of this script, it will replace the original database and zip files in the RESULTS folder.

resultsFolder=./RESULTS
workFolder=./workplace

(
 cd ../$resultsFolder
 for ijob in `ls *.zip | grep job`
  do
     echo processing $ijob ...
     unzip $ijob > /dev/null
     rm -f $ijob
     jobfolder=`echo $ijob | cut -f 1 -d .`
     (
      cd ..
      mv $resultsFolder/$jobfolder $workFolder
      (
       cd $workFolder/$jobfolder
       for iev in `ls | grep event`
        do
          echo processing $iev ...
          mv $iev ../iSS/results
          (cd ../iSS; ./iSS.e > /dev/null)
          mv ../iSS/results $iev
        done
      )
      (cd $workFolder/EbeCollector; echo collecting data into database ...; ./EbeCollectorShell_HydroEM.py ../$jobfolder)
      (cd $workFolder; echo compressing files ...; zip -r $ijob $jobfolder > /dev/null)
     )
     mv -f ../$workFolder/$ijob ./
  done
  rm -f collected.db
)

./combineEbeDatabasesFromZippedResults.py ../$resultsFolder
