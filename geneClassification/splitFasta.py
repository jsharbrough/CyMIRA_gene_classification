import sys
def splitFasta(fasta,numSeqs=1000):
    seqDict,seqList = buildSeqDict(fasta)
    i = 0
    fh = fasta[0:-6]
    j = 0
    fileList = []
    while i < len(seqList):
        while j < numSeqs and i < len(seqList):
            outfile = open(fh + '.' + str(i/numSeqs) + '.fasta','a')
            outfile.write(seqList[i] + seqDict[seqList[i]])
            outfile.close()
            fileName = (fh + '.' + str(i/numSeqs) + '.fasta\n')
            if fileName not in fileList:
                fileList.append(fileName)
            j += 1
            i += 1
        j = 0
    outfile = open(fh + '.targetp.fofn','w')
    for item in fileList:
        outfile.write(item)
    sys.stdout.write(fh + '.targetp.fofn')

def buildSeqDict(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line
            scaffoldList.append(seqName)
            currSeq = ''
        else:
            currSeq += line
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict, scaffoldList

if len(sys.argv)==3:
    splitFasta(sys.argv[1],sys.argv[2])
else:
    splitFasta(sys.argv[1])
