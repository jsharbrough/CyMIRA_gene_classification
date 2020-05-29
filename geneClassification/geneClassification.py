import sys
'''

geneClassification v3.0

USAGE:
    
    python geneClassification.py CyMIRA.txt Orthogroups.txt targetingPredictions.txt > CyMIRA+targeting.txt

    python geneClassification.py help
    
'''

def geneClassification(atFile,Orthogroups,genePredictions):
    CyMIRA = buildGeneClassDict(atFile) #build CyMIRA access dictionary
    catList = ['not','ot','ot-I','ot-NI','EC','mt','pt','mt-I','mt-NI','mt-EC','pt-I','pt-NI','pt-EC','mt-DNA_RRR','mt-Mito_TAT','mt-Mitoribosome','mt-Mitoribosome;Large_Subunit','mt-Mitoribosome;Small_Subunit','mt-OXPHOS','mt-OXPHOS;Complex_I','mt-OXPHOS;Complex_III','mt-OXPHOS;Complex_IV','mt-OXPHOS;Complex_V','mt-PPR','mt-Transcription_and_Transcript_Maturation','mt-Transcription_and_Transcript_Maturation;Intron_Splicing','mt-Transcription_and_Transcript_Maturation;RNA_Polymerase','mt-Transcription_and_Transcript_Maturation;Transcript_End_Processing','mt-Transcription_and_Transcript_Maturation;mTERF','mt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification','mt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification','mt-tRNA_Aminoacylation','pt-ACCase','pt-CLP','pt-Chlororibosome','pt-Chlororibosome;Large_Subunit','pt-Chlororibosome;Small_Subunit','pt-DNA_RRR','pt-PPR','pt-Photosynthesis','pt-Photosynthesis;ATP_Synthase','pt-Photosynthesis;Cytochrome_b6f','pt-Photosynthesis;NDH','pt-Photosynthesis;PSI','pt-Photosynthesis;PSII','pt-Photosynthesis;Rubisco','pt-Transcription_and_Transcript_Maturation','pt-Transcription_and_Transcript_Maturation;Intron_Splicing','pt-Transcription_and_Transcript_Maturation;RNA_Polymerase','pt-Transcription_and_Transcript_Maturation;Sigma_Factor','pt-Transcription_and_Transcript_Maturation;Transcript_End_Processing','pt-Transcription_and_Transcript_Maturation;mTERF','pt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification','pt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification','pt-tRNA_Aminoacylation']
    orthoDict = {} #build orthogroup dictionary (OG as key, gene list as value)
    infile = open(Orthogroups,'r')
    for line in infile:
        realLine = line
        while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
            realLine = realLine[0:-1]
        lineSplit = realLine.split(' ')
        og = lineSplit[0]
        orthoDict[og[0:-1]] = lineSplit[1:]
    infile.close()
    backwardsOGDict = {} #Associate genes with orthogroup (Gene as key, OG as value)
    for og in orthoDict:
        currList = orthoDict[og]
        for gene in currList:
            backwardsOGDict[gene] = og
    predDict = {}
    geneList = []
    infile = open(genePredictions,'r')
    for line in infile:
            if line[0] != '#':
                realLine = line
                while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
                    realLine = realLine[0:-1]
                lineSplit = realLine.split('\t')
                predDict[lineSplit[0]] = lineSplit[1]
                geneList.append(lineSplit[0])
    infile.close()
    newCymira = {'pt-Chlororibosome':[],'pt-Transcription_and_Transcript_Maturation':[],'pt-Photosynthesis':[],'mt-Transcription_and_Transcript_Maturation':[],'mt-OXPHOS':[],'mt-Mitoribosome':[],'mt-OXPHOS;Complex_V': [], 'pt-Chlororibosome;Small_Subunit': [], 'pt-Photosynthesis;Cytochrome_b6f': [], 'mt-Transcription_and_Transcript_Maturation;Intron_Splicing': [], 'pt-Transcription_and_Transcript_Maturation;Intron_Splicing': [], 'mt-Transcription_and_Transcript_Maturation;RNA_Polymerase': [], 'mt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification': [], 'mt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification': [], 'pt-ACCase': [], 'pt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification': [], 'mt-OXPHOS;Complex_I': [], 'mt-Transcription_and_Transcript_Maturation;Transcript_End_Processing': [], 'pt': [], 'not': [], 'mt-Mitoribosome;Large_Subunit': [], 'pt-Photosynthesis;PSII': [], 'pt-Chlororibosome;Large_Subunit': [], 'pt-PPR': [], 'pt-Transcription_and_Transcript_Maturation;Sigma_Factor': [], 'pt-Photosynthesis;PSI': [], 'pt-DNA_RRR': [], 'pt-Transcription_and_Transcript_Maturation;mTERF': [], 'mt-Mitoribosome;Small_Subunit': [], 'mt-NI': [], 'mt-OXPHOS;Complex_III': [], 'pt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification': [], 'mt-PPR': [], 'pt-I': [], 'mt-I': [], 'mt-DNA_RRR': [], 'ot-NI': [], 'mt-Mito_TAT': [], 'ot-I': [], 'pt-Photosynthesis;ATP_Synthase': [], 'pt-Transcription_and_Transcript_Maturation;RNA_Polymerase': [], 'mt-Transcription_and_Transcript_Maturation;mTERF': [], 'mt-OXPHOS;Complex_IV': [], 'pt-tRNA_Aminoacylation': [], 'pt-NI': [], 'pt-Photosynthesis;Rubisco': [], 'pt-Transcription_and_Transcript_Maturation;Transcript_End_Processing': [], 'mt': [], 'pt-CLP': [], 'mt-tRNA_Aminoacylation': [], 'ot': [], 'pt-Photosynthesis;NDH': []}
    for gene in geneList:
        atDict = {}
        if gene in backwardsOGDict:
            og = backwardsOGDict[gene]
            atList = []
            ogList = orthoDict[og]
            for ortholog in ogList:
                if 'AT' == ortholog[0:2]:
                    ogSplit = ortholog.split('.')
                    atList.append(ogSplit[0])
            for at in atList:
                for cat in catList:
                    if at in CyMIRA[cat]:
                        if cat in atDict:
                            atDict[cat] += 1
                        else:
                            atDict[cat] = 1
            if 'ot' in atDict:
                if 'ot-I' in atDict:
                    if atDict['ot-I']/float(len(atList)) >= 0.5:
                        if 'ot-NI' in atDict:
                            del atDict['ot-NI']
                        if 'mt-I' in atDict and 'mt-NI' in atDict:
                            del atDict['mt-NI']
                        if 'pt-I' in atDict and 'pt-NI' in atDict:
                            del atDict['pt-NI']
                        if 'not' in atDict:
                            del atDict['not']
                    elif atDict['ot']/float(len(atList)) >= 0.5 and predDict[gene] != 'non-organellar':
                        if predDict[gene] == 'mitochondria':
                            if 'mt' not in atDict:
                                atDict['mt'] = 1
                        if predDict[gene] == 'plastid':
                            if 'pt' not in atDict:
                                atDict['pt'] = 1
                        if predDict[gene] == 'dual':
                            if 'mt' not in atDict:
                                atDict['mt'] = 1
                            if 'pt' not in atDict:
                                atDict['pt'] = 1
                        if 'ot-I' not in atDict and 'ot-NI' not in atDict:
                            atDict['ot-NI'] = 1
                        if 'mt' in atDict and 'mt-I' not in atDict and 'mt-NI' not in atDict:
                            atDict['mt-NI'] = 1
                        if 'pt' in atDict and 'pt-I' not in atDict and 'pt-NI' not in atDict:
                            atDict['pt-NI'] = 1
                        if 'not' in atDict:
                            del atDict['not']
                    else:
                        atDict = {'not':1}
                elif atDict['ot']/float(len(atList)) >= 0.5 and predDict[gene] != 'non-organellar':
                    if predDict[gene] == 'mitochondria':
                        if 'mt' not in atDict:
                            atDict['mt'] = 1
                    if predDict[gene] == 'plastid':
                        if 'pt' not in atDict:
                            atDict['pt'] = 1
                    if predDict[gene] == 'dual':
                        if 'mt' not in atDict:
                            atDict['mt'] = 1
                        if 'pt' not in atDict:
                            atDict['pt'] = 1
                    if 'ot-I' not in atDict and 'ot-NI' not in atDict:
                        atDict['ot-NI'] = 1
                    if 'mt' in atDict and 'mt-I' not in atDict and 'mt-NI' not in atDict:
                        atDict['mt-NI'] = 1
                    if 'pt' in atDict and 'pt-I' not in atDict and 'pt-NI' not in atDict:
                        atDict['pt-NI'] = 1
                    if 'not' in atDict:
                        del atDict['not']
                else:
                    atDict = {'not':1}
            else:
                atDict = {'not':1}
        else:
            atDict = {'not':1}
        for cat in atDict:
            if cat != 'EC' and cat != 'mt-EC' and cat != 'pt-EC':
                currList = newCymira[cat]
                currList.append(gene)
                newCymira[cat] = currList
    ECs = []
    mt_ECs = []
    pt_ECs = []
    for item in newCymira['mt-Mitoribosome']:
        if item not in ECs:
            ECs.append(item)
        if item not in mt_ECs:
            mt_ECs.append(item)
    for item in newCymira['mt-OXPHOS']:
        if item not in ECs:
            ECs.append(item)
        if item not in mt_ECs:
            mt_ECs.append(item)
    for item in newCymira['mt-Mito_TAT']:
        if item not in ECs:
            ECs.append(item)
        if item not in mt_ECs:
            mt_ECs.append(item)
    for item in newCymira['pt-ACCase']:
        if item not in ECs:
            ECs.append(item)
        if item not in pt_ECs:
            pt_ECs.append(item)
    for item in newCymira['pt-CLP']:
        if item not in ECs:
            ECs.append(item)
        if item not in pt_ECs:
            pt_ECs.append(item)
    for item in newCymira['pt-Chlororibosome']:
        if item not in ECs:
            ECs.append(item)
        if item not in pt_ECs:
            pt_ECs.append(item)
    for item in newCymira['pt-Photosynthesis']:
        if item not in ECs:
            ECs.append(item)
        if item not in pt_ECs:
            pt_ECs.append(item)
    newCymira['EC'] = ECs
    newCymira['mt-EC'] = mt_ECs
    newCymira['pt-EC'] = pt_ECs
    sys.stdout.write('Not-organelle-targeted\tOrganelle-targeted\tOrganelle-targeted_Interacting\tOrganelle-targeted_Non-interacting\tEnzyme_Complexes\tMitochondria-targeted\tPlastid-targeted\tMitochondria-targeted_Interacting\tMitochondria-targeted_Non-interacting\tMitochondria_Enzyme_Complexes\tPlastid-targeted_Interacting\tPlastid-targeted_Non-interacting\tPlastid_Enzyme_Complexes')
    for item in catList[13:]:
        sys.stdout.write('\t' + item)
    sys.stdout.write('\n')
    i = 0
    while i < len(newCymira['not']):
        for cat in catList:
            if i < len(newCymira[cat]):
                currGenes = newCymira[cat]
                sys.stdout.write(currGenes[i] + '\t')
            else:
                sys.stdout.write('\t')
        sys.stdout.write('\n')
        i += 1
            
def buildGeneClassDict(atFile):
    infile = open(atFile,'r')
    CyMIRA = {'pt-Chlororibosome':[],'pt-Transcription_and_Transcript_Maturation':[],'pt-Photosynthesis':[],'mt-Transcription_and_Transcript_Maturation':[],'mt-OXPHOS':[],'mt-Mitoribosome':[],'mt-OXPHOS;Complex_V': [], 'pt-Chlororibosome;Small_Subunit': [], 'pt-Photosynthesis;Cytochrome_b6f': [], 'mt-Transcription_and_Transcript_Maturation;Intron_Splicing': [], 'pt-Transcription_and_Transcript_Maturation;Intron_Splicing': [], 'mt-Transcription_and_Transcript_Maturation;RNA_Polymerase': [], 'mt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification': [], 'mt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification': [], 'pt-ACCase': [], 'pt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification': [], 'mt-OXPHOS;Complex_I': [], 'mt-Transcription_and_Transcript_Maturation;Transcript_End_Processing': [], 'pt': [], 'not': [], 'mt-Mitoribosome;Large_Subunit': [], 'pt-Photosynthesis;PSII': [], 'pt-Chlororibosome;Large_Subunit': [], 'pt-PPR': [], 'pt-Transcription_and_Transcript_Maturation;Sigma_Factor': [], 'pt-Photosynthesis;PSI': [], 'pt-DNA_RRR': [], 'pt-Transcription_and_Transcript_Maturation;mTERF': [], 'mt-Mitoribosome;Small_Subunit': [], 'mt-NI': [], 'mt-OXPHOS;Complex_III': [], 'pt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification': [], 'mt-PPR': [], 'pt-I': [], 'mt-I': [], 'mt-DNA_RRR': [], 'ot-NI': [], 'mt-Mito_TAT': [], 'ot-I': [], 'pt-Photosynthesis;ATP_Synthase': [], 'pt-Transcription_and_Transcript_Maturation;RNA_Polymerase': [], 'mt-Transcription_and_Transcript_Maturation;mTERF': [], 'mt-OXPHOS;Complex_IV': [], 'pt-tRNA_Aminoacylation': [], 'pt-NI': [], 'pt-Photosynthesis;Rubisco': [], 'pt-Transcription_and_Transcript_Maturation;Transcript_End_Processing': [], 'mt': [], 'pt-CLP': [], 'mt-tRNA_Aminoacylation': [], 'ot': [], 'pt-Photosynthesis;NDH': []}
    for line in infile:
        realLine = line
        while realLine[-1] == '\t' or realLine[-1] == '\n' or realLine[-1] == '\r':
            realLine = realLine[0:-1]
        lineSplit = realLine.split('\t')
        if lineSplit[1] == 'Other' or lineSplit[1] == 'No Call':
            currList = CyMIRA['not']
            currList.append(lineSplit[0])
            CyMIRA['not'] = currList
        elif lineSplit[1] == 'Mitochondria':
            currList = CyMIRA['mt']
            currList.append(lineSplit[0])
            CyMIRA['mt'] = currList
            currList = CyMIRA['ot']
            currList.append(lineSplit[0])
            CyMIRA['ot'] = currList
            if lineSplit[2] == 'No':
                currList = CyMIRA['mt-NI']
                currList.append(lineSplit[0])
                CyMIRA['mt-NI'] = currList
                currList = CyMIRA['ot-NI']
                currList.append(lineSplit[0])
                CyMIRA['ot-NI'] = currList
            else:
                currList = CyMIRA['mt-I']
                currList.append(lineSplit[0])
                CyMIRA['mt-I'] = currList
                currList = CyMIRA['ot-I']
                currList.append(lineSplit[0])
                CyMIRA['ot-I'] = currList
                if lineSplit[3] == 'DNA-RRR':
                    currList = CyMIRA['mt-DNA_RRR']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-DNA_RRR'] = currList
                elif lineSplit[3] == 'PPR':
                    currList = CyMIRA['mt-PPR']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-PPR'] = currList
                elif lineSplit[3] == 'tRNA Aminoacylation':
                    currList = CyMIRA['mt-tRNA_Aminoacylation']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-tRNA_Aminoacylation'] = currList
                elif lineSplit[3] == 'Mito TAT Complex':
                    currList = CyMIRA['mt-Mito_TAT']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-Mito_TAT'] = currList
                elif lineSplit[3] == 'Mitoribosome':
                    currList = CyMIRA['mt-Mitoribosome']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-Mitoribosome'] = currList
                    if lineSplit[4] == 'Large Subunit':
                        currList = CyMIRA['mt-Mitoribosome;Large_Subunit']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Mitoribosome;Large_Subunit'] = currList
                    elif lineSplit[4] == 'Small Subunit':
                        currList = CyMIRA['mt-Mitoribosome;Small_Subunit']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Mitoribosome;Small_Subunit'] = currList
                elif lineSplit[3] == 'OXPHOS':
                    currList = CyMIRA['mt-OXPHOS']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-OXPHOS'] = currList
                    if lineSplit[4] == 'Complex I':
                        currList = CyMIRA['mt-OXPHOS;Complex_I']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-OXPHOS;Complex_I'] = currList
                    elif lineSplit[4] == 'Complex III':
                        currList = CyMIRA['mt-OXPHOS;Complex_III']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-OXPHOS;Complex_III'] = currList
                    elif lineSplit[4] == 'Complex IV':
                        currList = CyMIRA['mt-OXPHOS;Complex_IV']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-OXPHOS;Complex_IV'] = currList
                    elif lineSplit[4] == 'Complex V':
                        currList = CyMIRA['mt-OXPHOS;Complex_V']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-OXPHOS;Complex_V'] = currList
                elif lineSplit[3] == 'Transcription and Transcript Maturation':
                    currList = CyMIRA['mt-Transcription_and_Transcript_Maturation']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-Transcription_and_Transcript_Maturation'] = currList
                    if lineSplit[4] == 'Intron Splicing':
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;Intron_Splicing']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;Intron_Splicing'] = currList
                    elif lineSplit[4] == 'RNA Polymerase':
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;RNA_Polymerase']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;RNA_Polymerase'] = currList
                    elif lineSplit[4] == 'Transcript End Processing':
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;Transcript_End_Processing']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;Transcript_End_Processing'] = currList
                    elif lineSplit[4] == 'mTERF':
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;mTERF']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;mTERF'] = currList
                    elif lineSplit[4] == 'rRNA Base Modification':
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification'] = currList
                    elif lineSplit[4] == 'tRNA Base Modification':
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification'] = currList
        elif lineSplit[1] == 'Plastid':
            currList = CyMIRA['pt']
            currList.append(lineSplit[0])
            CyMIRA['pt'] = currList
            currList = CyMIRA['ot']
            currList.append(lineSplit[0])
            CyMIRA['ot'] = currList
            if lineSplit[2] == 'No':
                currList = CyMIRA['pt-NI']
                currList.append(lineSplit[0])
                CyMIRA['pt-NI'] = currList
                currList = CyMIRA['ot-NI']
                currList.append(lineSplit[0])
                CyMIRA['ot-NI'] = currList
            else:
                currList = CyMIRA['pt-I']
                currList.append(lineSplit[0])
                CyMIRA['pt-I'] = currList
                currList = CyMIRA['ot-I']
                currList.append(lineSplit[0])
                CyMIRA['ot-I'] = currList
                if lineSplit[3] == 'DNA-RRR':
                    currList = CyMIRA['pt-DNA_RRR']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-DNA_RRR'] = currList
                elif lineSplit[3] == 'PPR':
                    currList = CyMIRA['pt-PPR']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-PPR'] = currList
                elif lineSplit[3] == 'tRNA Aminoacylation':
                    currList = CyMIRA['pt-tRNA_Aminoacylation']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-tRNA_Aminoacylation'] = currList
                elif lineSplit[3] == 'ACCase':
                    currList = CyMIRA['pt-ACCase']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-ACCase'] = currList
                elif lineSplit[3] == 'CLP':
                    currList = CyMIRA['pt-CLP']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-CLP'] = currList
                elif lineSplit[3] == 'Chlororibosome':
                    currList = CyMIRA['pt-Chlororibosome']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-Chlororibosome'] = currList
                    if lineSplit[4] == 'Large Subunit':
                        currList = CyMIRA['pt-Chlororibosome;Large_Subunit']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Chlororibosome;Large_Subunit'] = currList
                    elif lineSplit[4] == 'Small Subunit':
                        currList = CyMIRA['pt-Chlororibosome;Small_Subunit']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Chlororibosome;Small_Subunit'] = currList
                elif lineSplit[3] == 'Photosynthesis':
                    currList = CyMIRA['pt-Photosynthesis']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-Photosynthesis'] = currList
                    if lineSplit[4] == 'ATP Synthase':
                        currList = CyMIRA['pt-Photosynthesis;ATP_Synthase']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;ATP_Synthase'] = currList
                    elif lineSplit[4] == 'Cytochrome b6f':
                        currList = CyMIRA['pt-Photosynthesis;Cytochrome_b6f']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;Cytochrome_b6f'] = currList
                    elif lineSplit[4] == 'NDH':
                        currList = CyMIRA['pt-Photosynthesis;NDH']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;NDH'] = currList
                    elif lineSplit[4] == 'PSI':
                        currList = CyMIRA['pt-Photosynthesis;PSI']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;PSI'] = currList
                    elif lineSplit[4] == 'PSII':
                        currList = CyMIRA['pt-Photosynthesis;PSII']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;PSII'] = currList
                    elif lineSplit[4] == 'Rubisco':
                        currList = CyMIRA['pt-Photosynthesis;Rubisco']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;Rubisco'] = currList
                elif lineSplit[3] == 'Transcription and Transcript Maturation':
                    currList = CyMIRA['pt-Transcription_and_Transcript_Maturation']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-Transcription_and_Transcript_Maturation'] = currList
                    if lineSplit[4] == 'Intron Splicing':
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;Intron_Splicing']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;Intron_Splicing'] = currList
                    elif lineSplit[4] == 'RNA Polymerase':
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;RNA_Polymerase']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;RNA_Polymerase'] = currList
                    elif lineSplit[4] == 'Transcript End Processing':
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;Transcript_End_Processing']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;Transcript_End_Processing'] = currList
                    elif lineSplit[4] == 'mTERF':
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;mTERF']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;mTERF'] = currList
                    elif lineSplit[4] == 'rRNA Base Modification':
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification'] = currList
                    elif lineSplit[4] == 'tRNA Base Modification':
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification'] = currList
                    elif lineSplit[4] == 'Sigma Factor':
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;Sigma_Factor']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;Sigma_Factor'] = currList
        elif lineSplit[1] == 'Dual':
            currList = CyMIRA['pt']
            currList.append(lineSplit[0])
            CyMIRA['pt'] = currList
            currList = CyMIRA['mt']
            currList.append(lineSplit[0])
            CyMIRA['mt'] = currList
            currList = CyMIRA['ot']
            currList.append(lineSplit[0])
            CyMIRA['ot'] = currList
            if lineSplit[2] == 'No':
                currList = CyMIRA['pt-NI']
                currList.append(lineSplit[0])
                CyMIRA['pt-NI'] = currList
                currList = CyMIRA['mt-NI']
                currList.append(lineSplit[0])
                CyMIRA['mt-NI'] = currList
                currList = CyMIRA['ot-NI']
                currList.append(lineSplit[0])
                CyMIRA['ot-NI'] = currList
            else:
                currList = CyMIRA['ot-I']
                currList.append(lineSplit[0])
                CyMIRA['ot-I'] = currList
                if lineSplit[3] == 'DNA-RRR':
                    currList = CyMIRA['mt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-I'] = currList
                    currList = CyMIRA['pt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-I'] = currList
                    currList = CyMIRA['pt-DNA_RRR']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-DNA_RRR'] = currList
                    currList = CyMIRA['mt-DNA_RRR']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-DNA_RRR'] = currList
                elif lineSplit[3] == 'PPR':
                    currList = CyMIRA['mt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-I'] = currList
                    currList = CyMIRA['pt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-I'] = currList
                    currList = CyMIRA['pt-PPR']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-PPR'] = currList
                    currList = CyMIRA['mt-PPR']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-PPR'] = currList
                elif lineSplit[3] == 'tRNA Aminoacylation':
                    currList = CyMIRA['mt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-I'] = currList
                    currList = CyMIRA['pt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-I'] = currList
                    currList = CyMIRA['pt-tRNA_Aminoacylation']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-tRNA_Aminoacylation'] = currList
                    currList = CyMIRA['mt-tRNA_Aminoacylation']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-tRNA_Aminoacylation'] = currList
                elif lineSplit[3] == 'Mito TAT Complex':
                    currList = CyMIRA['mt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-I'] = currList
                    currList = CyMIRA['pt-NI']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-NI'] = currList
                    currList = CyMIRA['mt-Mito_TAT']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-Mito_TAT'] = currList
                elif lineSplit[3] == 'Mitoribosome':
                    currList = CyMIRA['pt-NI']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-NI'] = currList
                    currList = CyMIRA['mt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-I'] = currList
                    currList = CyMIRA['mt-Mitoribosome']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-Mitoribosome'] = currList
                    if lineSplit[4] == 'Large Subunit':
                        currList = CyMIRA['mt-Mitoribosome;Large_Subunit']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Mitoribosome;Large_Subunit'] = currList
                    elif lineSplit[4] == 'Small Subunit':
                        currList = CyMIRA['mt-Mitoribosome;Small_Subunit']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Mitoribosome;Small_Subunit'] = currList
                elif lineSplit[3] == 'OXPHOS':
                    currList = CyMIRA['pt-NI']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-NI'] = currList
                    currList = CyMIRA['mt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-I'] = currList
                    currList = CyMIRA['mt-OXPHOS']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-OXPHOS'] = currList
                    if lineSplit[4] == 'Complex I':
                        currList = CyMIRA['mt-OXPHOS;Complex_I']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-OXPHOS;Complex_I'] = currList
                    elif lineSplit[4] == 'Complex III':
                        currList = CyMIRA['mt-OXPHOS;Complex_III']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-OXPHOS;Complex_III'] = currList
                    elif lineSplit[4] == 'Complex IV':
                        currList = CyMIRA['mt-OXPHOS;Complex_IV']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-OXPHOS;Complex_IV'] = currList
                    elif lineSplit[4] == 'Complex V':
                        currList = CyMIRA['mt-OXPHOS;Complex_V']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-OXPHOS;Complex_V'] = currList
                elif lineSplit[3] == 'ACCase':
                    currList = CyMIRA['mt-NI']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-NI'] = currList
                    currList = CyMIRA['pt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-I'] = currList
                    currList = CyMIRA['pt-ACCase']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-ACCase'] = currList
                elif lineSplit[3] == 'CLP':
                    currList = CyMIRA['mt-NI']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-NI'] = currList
                    currList = CyMIRA['pt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-I'] = currList
                    currList = CyMIRA['pt-CLP']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-CLP'] = currList
                elif lineSplit[3] == 'Chlororibosome':
                    currList = CyMIRA['mt-NI']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-NI'] = currList
                    currList = CyMIRA['pt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-I'] = currList
                    currList = CyMIRA['pt-Chlororibosome']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-Chlororibosome'] = currList
                    if lineSplit[4] == 'Large Subunit':
                        currList = CyMIRA['pt-Chlororibosome;Large_Subunit']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Chlororibosome;Large_Subunit'] = currList
                    elif lineSplit[4] == 'Small Subunit':
                        currList = CyMIRA['pt-Chlororibosome;Small_Subunit']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Chlororibosome;Small_Subunit'] = currList
                elif lineSplit[3] == 'Photosynthesis':
                    currList = CyMIRA['mt-NI']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-NI'] = currList
                    currList = CyMIRA['pt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-I'] = currList
                    currList = CyMIRA['pt-Photosynthesis']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-Photosynthesis'] = currList
                    if lineSplit[4] == 'ATP Synthase':
                        currList = CyMIRA['pt-Photosynthesis;ATP_Synthase']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;ATP_Synthase'] = currList
                    elif lineSplit[4] == 'Cytochrome b6f':
                        currList = CyMIRA['pt-Photosynthesis;Cytochrome_b6f']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;Cytochrome_b6f'] = currList
                    elif lineSplit[4] == 'NDH':
                        currList = CyMIRA['pt-Photosynthesis;NDH']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;NDH'] = currList
                    elif lineSplit[4] == 'PSI':
                        currList = CyMIRA['pt-Photosynthesis;PSI']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;PSI'] = currList
                    elif lineSplit[4] == 'PSII':
                        currList = CyMIRA['pt-Photosynthesis;PSII']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;PSII'] = currList
                    elif lineSplit[4] == 'Rubisco':
                        currList = CyMIRA['pt-Photosynthesis;Rubisco']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Photosynthesis;Rubisco'] = currList
                elif lineSplit[3] == 'Transcription and Transcript Maturation':
                    currList = CyMIRA['pt-I']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-I'] = currList
                    currList = CyMIRA['pt-Transcription_and_Transcript_Maturation']
                    currList.append(lineSplit[0])
                    CyMIRA['pt-Transcription_and_Transcript_Maturation'] = currList
                    currList = CyMIRA['mt-Transcription_and_Transcript_Maturation']
                    currList.append(lineSplit[0])
                    CyMIRA['mt-Transcription_and_Transcript_Maturation'] = currList
                    if lineSplit[4] == 'Intron Splicing':
                        currList = CyMIRA['mt-I']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-I'] = currList
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;Intron_Splicing']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;Intron_Splicing'] = currList
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;Intron_Splicing']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;Intron_Splicing'] = currList
                    elif lineSplit[4] == 'RNA Polymerase':
                        currList = CyMIRA['mt-I']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-I'] = currList
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;RNA_Polymerase']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;RNA_Polymerase'] = currList
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;RNA_Polymerase']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;RNA_Polymerase'] = currList
                    elif lineSplit[4] == 'Transcript End Processing':
                        currList = CyMIRA['mt-I']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-I'] = currList
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;Transcript_End_Processing']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;Transcript_End_Processing'] = currList
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;Transcript_End_Processing']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;Transcript_End_Processing'] = currList
                    elif lineSplit[4] == 'mTERF':
                        currList = CyMIRA['mt-I']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-I'] = currList
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;mTERF']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;mTERF'] = currList
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;mTERF']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;mTERF'] = currList
                    elif lineSplit[4] == 'rRNA Base Modification':
                        currList = CyMIRA['mt-I']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-I'] = currList
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification'] = currList
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;rRNA_Base_Modification'] = currList
                    elif lineSplit[4] == 'tRNA Base Modification':
                        currList = CyMIRA['mt-I']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-I'] = currList
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification'] = currList
                        currList = CyMIRA['mt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-Transcription_and_Transcript_Maturation;tRNA_Base_Modification'] = currList
                    elif lineSplit[4] == 'Sigma Factor':
                        currList = CyMIRA['mt-NI']
                        currList.append(lineSplit[0])
                        CyMIRA['mt-NI'] = currList
                        currList = CyMIRA['pt-Transcription_and_Transcript_Maturation;Sigma_Factor']
                        currList.append(lineSplit[0])
                        CyMIRA['pt-Transcription_and_Transcript_Maturation;Sigma_Factor'] = currList
    ECs = []
    mt_ECs = []
    pt_ECs = []
    for item in CyMIRA['mt-Mitoribosome']:
        if item not in ECs:
            ECs.append(item)
        if item not in mt_ECs:
            mt_ECs.append(item)
    for item in CyMIRA['mt-OXPHOS']:
        if item not in ECs:
            ECs.append(item)
        if item not in mt_ECs:
            mt_ECs.append(item)
    for item in CyMIRA['mt-Mito_TAT']:
        if item not in ECs:
            ECs.append(item)
        if item not in mt_ECs:
            mt_ECs.append(item)
    for item in CyMIRA['pt-ACCase']:
        if item not in ECs:
            ECs.append(item)
        if item not in pt_ECs:
            pt_ECs.append(item)
    for item in CyMIRA['pt-CLP']:
        if item not in ECs:
            ECs.append(item)
        if item not in pt_ECs:
            pt_ECs.append(item)
    for item in CyMIRA['pt-Chlororibosome']:
        if item not in ECs:
            ECs.append(item)
        if item not in pt_ECs:
            pt_ECs.append(item)
    for item in CyMIRA['pt-Photosynthesis']:
        if item not in ECs:
            ECs.append(item)
        if item not in pt_ECs:
            pt_ECs.append(item)
    CyMIRA['EC'] = ECs
    CyMIRA['mt-EC'] = mt_ECs
    CyMIRA['pt-EC'] = pt_ECs
    infile.close()
    return CyMIRA
    
helpStatement='\n\ngeneClassification v3.0\n\nUSAGE:\n\n\tpython geneClassification.py CyMIRA.txt Orthogroups.txt targetingPredictions.txt > CyMIRA+targeting.txt\n\n\tpython geneClassification.py help\n\n'
if 'help' in sys.argv or len(sys.argv) != 4:
    sys.stderr.write(helpStatement)    
else:
    geneClassification(sys.argv[1],sys.argv[2],sys.argv[3])

