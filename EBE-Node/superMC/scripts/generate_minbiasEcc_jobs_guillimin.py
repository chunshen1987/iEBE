#!/usr/bin/env python

import sys
from os import path, mkdir
import shutil
from glob import glob

def generate_script(folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    walltime = '10:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    script.write(
"""#!/usr/bin/env bash
#PBS -N %s
#PBS -l nodes=1:ppn=1
#PBS -l walltime=%s
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ad
#PBS -q sw
#PBS -m bea
#PBS -M chunshen1987@gmail.com
#PBS -d %s

rm -fr data
mkdir data
./superMC.e

""" % (working_folder.split('/')[-1], walltime, working_folder))
    script.close()


def generate_event_folder(working_folder, event_id):
    shutil.copytree('codes/superMC', 
        path.join(path.abspath(working_folder), 'superMC_%d' % event_id))
    event_folder = path.join(working_folder, 'superMC_%d' % event_id)
    generate_script(event_folder)

if __name__ == "__main__":
    try:
        folder_name = str(sys.argv[1])
        ncore = int(sys.argv[2])
    except IOError:
        print "./generate_jobs_guillimin.py working_folder num_of_cores"
        exit(0)

    for icore in range(ncore):
        generate_event_folder(folder_name, icore)

