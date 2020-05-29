import sys
def ipsortOutput(ipsort,fasta):
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
    infile = open(ipsort,'r')
    currLine = ''
    predDict = {}
    i = 0
    for line in infile:
        if 'Prediction:' in currLine:
            realLine = line
            while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
                realLine = realLine[0:-1]
            while realLine[0] == '\t':
                realLine = realLine[1:]
            if realLine == 'Other':
                prediction = 'non-organellar'
            elif realLine == 'Mitochondrial Transit Peptide':
                prediction = 'mitochondrial'
            elif realLine == 'Chloroplast Transit Peptide':
                prediction = 'plastid'
            elif realLine == 'Signal Peptide':
                prediction = 'non-organellar'
            predDict[geneList[i]] = prediction
            i += 1
        currLine = line
    infile.close()
    sys.stdout.write('#Gene\tiPSORT Prediction\n')
    for gene in geneList:
        sys.stdout.write(gene + '\t' + predDict[gene] + '\n')

ipsortOutput(sys.argv[1],sys.argv[2])
