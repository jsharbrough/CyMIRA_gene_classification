import sys
def targetpOutput(targetp,fasta):
    infile = open(fasta,'r')
    geneList = []
    for line in infile:
        if line[0] == '>':
            realLine = line
            while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
                realLine = realLine[0:-1]
            lineSplit = realLine.split(' ')
            gene = lineSplit[0]
            geneList.append(gene[1:])
    infile.close()
    infile = open(targetp,'r')
    predDict = {}
    i = 0
    lines = infile.readlines()
    infile.close()
    start = False
    while i < len(lines):
        currLine = lines[i]
        if start == False:
            if currLine == '----------------------------------------------------------------------\n':
                start = True
        elif currLine == '----------------------------------------------------------------------\n':
            start = False
        else:
            lineSplit = currLine.split(' ')
            while '' in lineSplit:
                lineSplit.remove('')
            pred = 'N/A'
            mt = float(lineSplit[3])
            pt = float(lineSplit[2])
            sp = float(lineSplit[4])
            other = float(lineSplit[5])
            if mt > sp and mt > other:
                if pt > sp and pt > other:
                    pred = 'dual'
                else:
                    pred = 'mitochondrial'
            elif pt > sp and pt > other:
                pred = 'plastid'
            else:
                pred = 'non-organellar'
            predDict[lineSplit[0]] = pred
        i += 1
    sys.stdout.write('#Gene\ttargetp Prediction\n')
    for gene in geneList:
        sys.stdout.write(gene + '\t' + predDict[gene] + '\n')

targetpOutput(sys.argv[1],sys.argv[2])
