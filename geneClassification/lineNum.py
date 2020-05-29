import sys
def lineNum(Infile):
    infile = open(Infile,'r')
    currLineNum = 0
    for line in infile:
        currLineNum += 1
    sys.stdout.write(str(currLineNum-1))

lineNum(sys.argv[1])
