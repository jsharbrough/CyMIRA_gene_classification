import sys
def getLine(fofn,lineNum):
    infile = open(fofn,'r')
    currLineNum = 0
    for line in infile:
        if currLineNum == lineNum:
            realLine = line
            while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
                realLine = realLine[0:-1]
            sys.stdout.write(realLine)
            return
        currLineNum += 1


getLine(sys.argv[1],int(sys.argv[2]))
