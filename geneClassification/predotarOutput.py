import sys
def predotarOutput(predotar,fasta):
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
    infile = open(predotar,'r')
    predDict = {}
    for line in infile:
        realLine = line
        while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
            realLine = realLine[0:-1]
        lineSplit = realLine.split('\t')
        if lineSplit[0] in geneList:
            pred = 'N/A'
            if 'Discarding' not in line:
                mt = float(lineSplit[1])
                pt = float(lineSplit[2])
                er = float(lineSplit[3])
                none = float(lineSplit[4])
                if mt > er and mt > none:
                    if pt > er and pt > none:
                        pred = 'dual'
                    else:
                        pred = 'mitochondrial'
                elif pt > er and pt > none:
                    pred = 'plastid'
                else:
                    pred = 'non-organellar'
            predDict[lineSplit[0]] = pred
    infile.close()
    sys.stdout.write('#Gene\tpredotar Prediction\n')
    for gene in geneList:
        sys.stdout.write(gene + '\t' + predDict[gene] + '\n')

predotarOutput(sys.argv[1],sys.argv[2])
