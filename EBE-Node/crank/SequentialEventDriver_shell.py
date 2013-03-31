#!/usr/bin/env python

from sys import argv, exit
try:
    numberOfEvents = int(argv[1])
except:
    print("Usage: SequentialEventDriver_shell.py number_of_events")
    exit()

# call the real shell
import SequentialEventDriver
SequentialEventDriver.controlParameterList["numberOfEvents"] = numberOfEvents
SequentialEventDriver.sequentialEventDriverShell()
