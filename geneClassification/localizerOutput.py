import sys
def localizerOutput(localizer,fasta):
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
    infile = open(localizer,'r')
    lines = infile.readlines()
    predDict = {}
    for line in lines[4:-24]:
        realLine = line
        while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
            realLine = realLine[0:-1]
        lineSplit = realLine.split('\t')
        seqNameSplit = lineSplit[0].split(' ')
        seqName = seqNameSplit[0]
        lineEnd = lineSplit[1]
        for item in lineSplit[2:]:
            lineEnd += '\t' + item
        lineEnd.replace(' ','')
        lineSplit = lineEnd.split('\t')
        if seqName in geneList:
            pred = 'non-organellar'
            if 'Y' in lineSplit[0]:
                if 'Y' in lineSplit[1]:
                    pred = 'dual'
                else:
                    pred = 'plastid'
            elif 'Y' in lineSplit[1]:
                pred = 'mitochondrial'
            predDict[seqName] = pred
    infile.close()
    sys.stdout.write('#Gene\tlocalizer Prediction\n')
    for gene in geneList:
        sys.stdout.write(gene + '\t' + predDict[gene] + '\n')

localizerOutput(sys.argv[1],sys.argv[2])
