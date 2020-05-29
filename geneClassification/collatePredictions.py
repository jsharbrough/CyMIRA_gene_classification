import sys
def collatePredictions(b):
    ipsortDict = predDict(b + '.processed.ipsort.out')
    targetpDict = predDict(b + '.processed.targetp.out')
    localizerDict = predDict(b + '.processed.LOCALIZER.out')
    predotarDict = predDict(b + '.processed.predotar.out')
    infile = open(b + '.fasta','r')
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
    sys.stdout.write('#Gene\tOverall Prediction\tTargetP\tiPSORT\tLOCALIZER\tPredotar\n')
    for gene in geneList:
        predictions = [targetpDict[gene],ipsortDict[gene],localizerDict[gene],predotarDict[gene]]
        if 'dual' in predictions:
            overallPrediction = 'dual'
        elif 'mitochondrial' in predictions:
            if 'plastid' in predictions:
                overallPrediction = 'dual'
            else:
                overallPrediction = 'mitochondrial'
        elif 'plastid' in predictions:
            overallPrediction = 'plastid'
        else:
            overallPrediction = 'non-organellar'
        sys.stdout.write(gene + '\t' + overallPrediction + '\t' + targetpDict[gene] + '\t' + ipsortDict[gene] + '\t' + localizerDict[gene] + '\t' + predotarDict[gene] + '\n')
    

def predDict(predictions):
    infile = open(predictions,'r')
    predictionDict = {}
    for line in infile:
        realLine = line
        while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
            realLine = realLine[0:-1]
        lineSplit = realLine.split('\t')
        predictionDict[lineSplit[0]] = lineSplit[1]
    return predictionDict

collatePredictions(sys.argv[1])
