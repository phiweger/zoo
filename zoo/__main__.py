'''
zoo command line.
'''

from .commands import load
import sys


def main():

    cmd = load
    cmd(sys.argv[2:])  # zoo [0] cmd [1]

# $ zoo
# Wow
