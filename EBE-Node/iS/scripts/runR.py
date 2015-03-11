###################    Last edited on Feb, 11 2010   #################################
#                                                              version 21 --- Zhi Qiu
# level 2
"""
	Provide functions to run external programs.
"""

from os import getcwd
import subprocess

__q_debug__ = False

def run(command, cwd=getcwd()):
	""" Invoke a command and wait for it to stop. """
	proc = subprocess.Popen(command, shell=True, cwd=cwd)
	while proc.communicate() != (None,None):
		pass

if __name__ == "__main__":
	print("Morning!")
