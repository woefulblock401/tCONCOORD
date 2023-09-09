import sys,os

l = open(sys.argv[1],'r').read()
n = l.replace('concoord.h','tconcoord.h')
print n
