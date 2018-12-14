import matplotlib.pyplot as plt
import os
from collections import defaultdict
import subprocess
import numpy
import pickle
from plumbum import local
import random
from fractions import Fraction

def readByteObjectFromDefinedDirectory(directory):
    f = open(directory, 'rb')
    return pickle.load(f)

def writeByteObjectToDefinedDirectory(directory,object):
    f = open(directory , 'wb+')
    pickle.dump(object, f)

def writeListToDestination(destination, listToWrite):
    #print('Writing list to ' + destination)
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(listToWrite):
            if i != len(listToWrite) - 1:
                writer.write(listToWrite[i] + '\n')
            elif i == len(listToWrite) - 1:
                writer.write(listToWrite[i])
            i += 1

def readDefinedFileToList(filename):
    tempList = []
    with open(filename, mode='r') as reader:
        tempList = [line.rstrip() for line in reader]
    return tempList

def createDictFromFasta(fastaList):
    tempDict = {}
    i = 0
    while i < len(fastaList):
        tempDict[fastaList[i][1:]] = fastaList[i + 1]
        i += 2
    return tempDict



################ FIGURE 2 ################
def figure2RawSeqs():
    '''
    This will be the code for creating the top stacked bar graph that will show the raw sequences before
    under going QC and MED processing

    I have the choice of having this plot have already QCed data in it (not MEDed yet) or have it be completely raw.
    Given that the QC is part of SP and that in that QC we get rid of singletons I think it is a good idea to have these
    sequences be completely raw.

    The next set of plots in this sequence will have been through the QC and through MED (only in the case of SP)
    i.e. not for the OTU analyses

    I think that the colour coding of the seqs can be somewhat non-descript. I'm thinking just using repetitive
    greys for the differnt sequences. We will plot the sample by sample. We will colour specific sequences that will
    carry on through the whole figure i.e. the differnt plots.

    Probably a good idea to have a look at the SymPortal output of the Ed data and look to see which the DIVs are
    that way we can see what the actual sequences are and we can identify them and name them in the figure.

    Pseudo-code:

    get the directory that holds all of the fastas

    for each fasta in the directory:


    '''

    # # Convert all of the Ed sequences so that they are no long reverse complement
    # baseDir = '/home/humebc/projects/SymPortalMS/edData'
    # listOfFiles = []
    # for (dirpath, dirnames, filenames) in os.walk(baseDir):
    #     listOfFiles.extend(filenames)
    #     break
    # listOfFiles = sorted([file for file in listOfFiles if '.fasta' in file])
    #
    # for file in listOfFiles:
    #     localFasta = []
    #     with open('{}/{}'.format(baseDir, file)) as f:
    #         localFasta = f.read().splitlines()
    #
    #     for i in range(len(localFasta)):
    #         if localFasta[i][0] == '>':
    #             localFasta[i+1] = reverse_complement(localFasta[i+1])
    #
    #     writeListToDestination('{}/{}'.format(baseDir, file),localFasta)

    # Read in the sequences that will need to be coloured as they are the type defining seqs
    significantSeqsFasta = []

    with open('/home/humebc/projects/SymPortalMS/significantSeqs') as f:
        significantSeqsFasta = f.read().splitlines()

    # Convert the sigSeqs to a dict
    sigSeqDict = {}
    for i in range(len(significantSeqsFasta)):
        if significantSeqsFasta[i][0] == '>':
            sigSeqDict[significantSeqsFasta[i][1:]] = significantSeqsFasta[i+1]

    baseDir = '/home/humebc/projects/SymPortalMS/edData'

    listOfFiles = []
    for (dirpath, dirnames, filenames) in os.walk(baseDir):
        listOfFiles.extend(filenames)
        break
    listOfFiles = sorted([file for file in listOfFiles if '.fasta' in file])

    # We esentially need to make something like a name and fasta combo.
    # Lets get all of the sequences found in all of the samples and hold them in a default dict
    # of sequence: abundance
    # Then we can sort this
    # Then we can work in this order.
    # we can have a second dict that is the sample, followed by a list that will hold the abundance of each of the
    # sequences in order of the sorted list.
    # to get this we will go for every sequence in sorted list
    # for every fasta in the basedir
    # create a temp default dict that is the seqs as keys and 0 as values.
    # then for each sequence add this to the temp dict
    # then write pickle out this dict

    # a dict that is sequence key and then a list for value which has 0 for each fasta
    # then for each fasta,
    # for seq value in the fasta's dict
    # dict[seq][filelist.index(fasta)] = value

    # then we will have the dict populated in a way that will let us plot the bars

    # for key value pair in dict plot bars


    # First go through all fastas and find all of the seqs and their abundances over all of the samples
    # we will hold this in a dictionary
    #TODO This is really slow and should do unique seqs
    seqTotalsDict = defaultdict(float)
    for file in listOfFiles:
        print(file)
        localFasta = []
        with open('{}/{}'.format(baseDir, file)) as f:
            localFasta = f.read().splitlines()

        # This returns a dict of the clade C seqs and their abundance
        CAbundDict = blastForC(localFasta)

        # It will take quite some time creating the C only fastas so lets print them out


        totalSeqs = sum(CAbundDict.values())
        for key in CAbundDict.keys():
                seqTotalsDict[key] += ((CAbundDict[key]/totalSeqs) * 100)

    # Here we can blast the dict so that we are only working with the C sequences
    # then further down when we get the abundances of the sequences in each of the fasta we will only
    # count on sequences that are in this master list that is only the clade C seqs


    # Here we have the dict populated with seq as key and the abundance across all seqs
    # Get a sorted list of the most abundant seqs
    sortedSeqs = [k for k in sorted(seqTotalsDict, key=seqTotalsDict.get, reverse=True)]

    # now per file
    plot_info_dict = {seq:[0 for file in listOfFiles] for seq in sortedSeqs}
    print('Populating plot info Dict')
    for file in listOfFiles:
        print('Sample: {}'.format(file))
        tempDict = defaultdict(float)
        localFasta = []
        with open('{}/{}'.format(baseDir, file)) as f:
            localFasta = f.read().splitlines()

        localCAbundDict = blastForC(localFasta)
        #Debug#
        firstTotal = sum(localCAbundDict.values())
        #####
        # When we wer doing the consolidating we were losing the C3 seqs,
        # Also for this raw figure seqs figure we don't really want to be consolidating as we want
        # to show that this is what MED is good for.
        # consolidateSeqs(localCAbundDict, sigSeqDict, file)
        totalSeqs = sum(localCAbundDict.values())
        for key in localCAbundDict.keys():
                tempDict[key] += ((localCAbundDict[key] / totalSeqs) * 100)
        # here we have the local tempDict populated
        # Now transfer that data into the plot_info_dict
        for seq, abund in tempDict.items():
            plot_info_dict[seq][listOfFiles.index(file)] = abund

    # Here we should have all of the plotting info done.



    colourPalette = ['#E8AB52', '#5FA96E', '#9680D3', '#D8C542', '#D47070', '#59A6D4', '#D76FBC', '#F6724C',
                     '#5FC2B6', '#80D1E4', '#6BD78A', '#B7D456']
    greyPalette = ['#D0CFD4', 	'#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

    colourDict = {a: b for a, b in zip(sigSeqDict.keys(), colourPalette)}
    # Set up plot
    plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1,1,1)
    ax.set_ylim(bottom=0, top=100)
    # 8 = plot, plot, space, plot, plot, space, plot, plot
    ind = range(len(listOfFiles))
    width = 1

    # Counters to work out which colour or grey we are on
    colourCounter = 0
    greyCounter = 0

    # The bottom for each sequence plot will be updated each time we plot
    bottom = [0 for file in listOfFiles]

    # Do the plotting
    plotColour = None
    # Need to plot in the order of the sorted list
    # Only plot the first 250 seqs as this is the majority of the volume
    # After that plot the rest of the volume as fixed height boxes
    # NB the greyFileHeight should be the same as the max value (or not smaller) of the max random generated
    greyFileHeight = 0.1
    for seq in sortedSeqs:
        seqIndex = sortedSeqs.index(seq)
        if seqIndex < 250:
            print('Plotting seq {} out of {}'.format(seqIndex, len(sortedSeqs)))
            # Work out whether the seq should be colour or grey
            matchTup = checkIfSeqSig(seq, sigSeqDict)
            if matchTup[0]:
                # then this is a sig seq and it should be a colour
                # Only plot the colours once, even if the lesser sequences match a sigseq again
                if colourCounter in range(len(colourPalette)):
                    plotColour = colourDict[matchTup[1]]
                    colourCounter += 1
                else:
                    plotColour = greyPalette[greyCounter % 4]
                    greyCounter += 1
            else:
                plotColour = greyPalette[greyCounter % 4]
                greyCounter += 1

            ax.bar(ind, plot_info_dict[seq], width, bottom=bottom, color=plotColour)
            bottom = [L + M for L, M in zip(bottom, plot_info_dict[seq])]
        else:

            #Make a custeom values list
            print('Plotting seq {} out of {}'.format(seqIndex, len(sortedSeqs)))
            plotValuesList = []
            for bValue in bottom:
                if bValue < 100:
                    # If this is not the last bar then plot a standard bar
                    if bValue < 100 - greyFileHeight:
                        plotValuesList.append(numpy.random.choice(numpy.arange(0.01, 0.1, 0.01), p=[0.1, 0.3, 0.3, 0.1, 0.1, 0.025, 0.025, 0.025, 0.025]))
                    # If we're onto the last bar for the sample, take it up to 100
                    elif bValue > 100 - greyFileHeight:
                        plotValuesList.append(100-bValue)
                else:
                    plotValuesList.append(0)
            if sum(plotValuesList) == 0:
                break
            # now plot
            if checkIfSeqSig(seq, sigSeqDict):
                # then this is a sig seq and it should be a colour
                try:
                    plotColour = colourPalette[colourCounter]
                    colourCounter += 1
                except:
                    plotColour = greyPalette[greyCounter % 4]
                    greyCounter += 1
            else:
                plotColour = greyPalette[greyCounter % 4]
                greyCounter += 1

            ax.bar(ind, plotValuesList, width, bottom=bottom, color=plotColour)
            bottom = [L + M for L, M in zip(bottom, plotValuesList)]

    # Here we should have it plotted and we'll just need to put the finishing touches to it.
    # we should debug to here
    debug = 'to here'
    return

def figure2_afterMED():
    apples = 'aser'
    baseDir = '/home/humebc/projects/SymPortalMS/edData'
    listOfDirs = []
    for (dirpath, dirnames, filenames) in os.walk(baseDir):
        listOfDirs.extend(dirnames)
        break
    listOfDirs = sorted([dir for dir in listOfDirs if dir not in  ['transFiles', 'arifOTU']])

    # read in the significant sequences again
    # Read in the sequences that will need to be coloured as they are the type defining seqs
    significantSeqsFasta = []

    with open('/home/humebc/projects/SymPortalMS/significantSeqs') as f:
        significantSeqsFasta = f.read().splitlines()

    # Convert the sigSeqs to a dict
    sigSeqDict = {}
    for i in range(len(significantSeqsFasta)):
        if significantSeqsFasta[i][0] == '>':
            sigSeqDict[significantSeqsFasta[i][1:]] = significantSeqsFasta[i + 1]



    totalSeqsDict = defaultdict(float)
    for dir in listOfDirs:
        # fasta of node to node sequences
        nodeFasta = readDefinedFileToList('{}/{}/NODE-REPRESENTATIVES.fasta'.format(baseDir, dir))
        nodeDict = {}
        for i in range(len(nodeFasta)):
            if nodeFasta[i][0] == '>':
                nodeDict[nodeFasta[i][1:].split('|')[0]] = nodeFasta[i + 1].replace('-','')

        # abundances of each node in the sample as dict
        countTab = readDefinedFileToList('{}/{}/MATRIX-COUNT.txt'.format(baseDir, dir))
        nodes = countTab[0].split('\t')[1:]
        counts = [int(x) for x in countTab[1].split('\t')[1:]]
        total = sum(counts)
        nodeAbundDict = {}
        for i in range(len(nodes)):
            nodeAbundDict[nodes[i]] = counts[i]

        # convert the nodes to abund dict to a seq to abund dict
        seqToAbundDict = {}
        for key, value in nodeAbundDict.items():
            seqToAbundDict[nodeDict[key]] = (value / total) * 100

        # for the sequences that are matches to the sigSeqs, make sure that any seqs that are subsets of each other
        # are collapsed
        # for each of the significatn sequences
        # get a list of the localseqs that are a subset or vice versa
        # then collapse these seqs keeping the longest sequence as the representative and the abundance of the totals
        for sigName, sigSeq in sigSeqDict.items():
            seqsToCollapse = []
            for localSeq, localAbund in seqToAbundDict.items():
                if sigSeq in localSeq or localSeq in sigSeq:
                    seqsToCollapse.append(localSeq)
            # Here we have a list of seqs that need to be collapsed
            if len(seqsToCollapse) > 1:
                newTot = sum([seqToAbundDict[seqName] for seqName in seqsToCollapse])
                # delete all seqs and replace with the C3 sig seq
                for localSeq in seqsToCollapse:
                        del seqToAbundDict[localSeq]
                # add the sigSeq in question and the value
                seqToAbundDict[sigSeq] = newTot




        #Pickle this dict out so that we don't have to remake it
        writeByteObjectToDefinedDirectory('{}/{}/dir.seqToAbundDict'.format(baseDir, dir), seqToAbundDict)

        #Add this samples sequences to the total seqs dict
        for key, value in seqToAbundDict.items():
            totalSeqsDict[key] += value

    # Here we have the total seqsDict
    # Now sort the list and populate the plot info list according to this order
    sortedSeqs = [k for k in sorted(totalSeqsDict, key=totalSeqsDict.get, reverse=True)]

    # now per file
    plot_info_dict = {seq: [0 for dir in listOfDirs] for seq in sortedSeqs}

    # now populate the plot_info_dict from the sample by sample dicts that we pickled earlier
    for dir in listOfDirs:
        seqToAbundDict = readByteObjectFromDefinedDirectory('{}/{}/dir.seqToAbundDict'.format(baseDir, dir))

        for seq, abund in seqToAbundDict.items():
            plot_info_dict[seq][listOfDirs.index(dir)] = abund

    # and now we are ready to plot using that info


    # Now plot

    # Colour palettes
    colourPalette = ['#E8AB52', '#5FA96E', '#9680D3', '#D8C542', '#D47070', '#59A6D4', '#D76FBC', '#F6724C',
                     '#5FC2B6', '#80D1E4', '#6BD78A', '#B7D456']
    greyPalette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

    # Colour palette dict
    colourDict = {a:b for a, b in zip(sigSeqDict.keys(), colourPalette)}
    # Set up plot
    plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1, 1, 1)
    ax.set_ylim(bottom=0, top=100)
    # 8 = plot, plot, space, plot, plot, space, plot, plot
    ind = range(len(listOfDirs))
    width = 1

    # Counters to work out which colour or grey we are on
    colourCounter = 0
    greyCounter = 0

    # The bottom for each sequence plot will be updated each time we plot
    bottom = [0 for file in listOfDirs]

    # Do the plotting
    plotColour = None
    # Need to plot in the order of the sorted list
    # Only plot the first 250 seqs as this is the majority of the volume
    # After that plot the rest of the volume as fixed height boxes
    greyFileHeight = 0.5
    plottedSigSeqs = []
    for seq in sortedSeqs:
        seqIndex = sortedSeqs.index(seq)
        if seqIndex < 250:
            # print('Plotting seq {} out of {}'.format(seqIndex, len(sortedSeqs)))
            # Work out whether the seq should be colour or grey
            matchTup = checkIfSeqSig(seq, sigSeqDict)
            if matchTup[0]:
                # then this is a sig seq and it should be a colour
                # Only plot the colours once, even if the lesser sequences match a sigseq again
                if matchTup[1] not in plottedSigSeqs:
                    plotColour = colourDict[matchTup[1]]
                    plottedSigSeqs.append(matchTup[1])
                    print('{}: {}'.format(matchTup[1], colourDict[matchTup[1]]))
                else:
                    plotColour = greyPalette[greyCounter % 4]
                    greyCounter += 1
            else:
                plotColour = greyPalette[greyCounter % 4]
                greyCounter += 1

            ax.bar(ind, plot_info_dict[seq], width, bottom=bottom, color=plotColour)
            bottom = [L + M for L, M in zip(bottom, plot_info_dict[seq])]
        else:
            # Make a custeom values list
            print('Plotting seq {} out of {}'.format(seqIndex, len(sortedSeqs)))
            plotValuesList = []
            for bValue in bottom:
                if bValue < 100:
                    # If this is not the last bar then plot a standard bar
                    if bValue < 100 - greyFileHeight:
                        plotValuesList.append(numpy.random.choice(numpy.arange(0.05, 0.5, 0.05),
                                                                  p=[0.1, 0.3, 0.3, 0.1, 0.1, 0.025, 0.025, 0.025,
                                                                     0.025]))
                    # If we're onto the last bar for the sample, take it up to 100
                    elif bValue > 100 - greyFileHeight:
                        plotValuesList.append(100 - bValue)
                else:
                    plotValuesList.append(0)
            if sum(plotValuesList) == 0:
                break
            # now plot
            if checkIfSeqSig(seq, sigSeqDict):
                # then this is a sig seq and it should be a colour
                try:
                    plotColour = colourPalette[colourCounter]
                    colourCounter += 1
                except:
                    plotColour = greyPalette[greyCounter % 4]
                    greyCounter += 1
            else:
                plotColour = greyPalette[greyCounter % 4]
                greyCounter += 1

            ax.bar(ind, plotValuesList, width, bottom=bottom, color=plotColour)
            bottom = [L + M for L, M in zip(bottom, plotValuesList)]
    return

def arifClustering():
    '''This is the code we were writing to do the cluster analysis. I wrote code that works that gets a collection of
    all of the difference sequences that were returned from all of the samples. I then align this in mafft using
    the auto function. I then use mothur to create a distance matrix and finally I cluster and get representatives
    for the clustered sequences.'''
    # we need to get a master list of sequences with one representative of each sequence
    # then align in mafft
    # then dist
    # then cluster

    baseDir = '/home/humebc/projects/SymPortalMS/edData'
    listOfDirs = []
    for (dirpath, dirnames, filenames) in os.walk(baseDir):
        listOfDirs.extend(dirnames)
        break
    listOfDirs = sorted([dir for dir in listOfDirs if dir not in  ['transFiles', 'arifOTU']])


    ####### I HAVE COMMENTED THIS OUT AS IT ENDED UP AS ALL OF THE SEQS BEING CLUSTERED TOGETHER ########
    ### The code does work though ###
    # globalSeqSet = set()
    # for dir in listOfDirs:
    #     localFasta = readDefinedFileToList('{}/{}/{}.CladeC.good.unique.abund.pcr.good.unique.fasta'.format(baseDir, dir, dir))
    #
    #     for i in range(len(localFasta)):
    #         if localFasta[i][0] == '>':
    #             globalSeqSet.add(localFasta[i + 1])
    #
    # # Here we have a globalSetOfSequences
    # # convert this into a fasta that we can cluster
    # arifFasta = []
    # sequenceCounter = 0
    # for seq in globalSeqSet:
    #     arifFasta.extend(['>seq{}'.format(sequenceCounter), seq])
    #     sequenceCounter += 1
    #
    # # here we have a global fasta
    # # Now make a directory where we can do all of this work so we don't mess up the baseDir
    # os.makedirs('{}/arifOTU'.format(baseDir), exist_ok=True)
    #
    # # write out the fasta
    # writeListToDestination('{}/arifOTU/master.fasta'.format(baseDir), arifFasta)
    #
    # # now align
    # # http://plumbum.readthedocs.io/en/latest/local_commands.html#pipelining
    # mafft = local["mafft"]
    # inFile = '{}/arifOTU/master.fasta'.format(baseDir)
    # outFile = '{}/arifOTU/master.aligned.fasta'.format(baseDir)
    #
    # # now run mafft including the redirect
    # (mafft['--auto', '--thread', '16', inFile] > outFile)()
    #
    # mothurBatch = ['dist.seqs(fasta={}/arifOTU/master.aligned.fasta, countends=F, output=lt)'.format(baseDir),
    #                'cluster(phylip={}/arifOTU/master.aligned.phylip.dist, method=average, cutoff=0.03)'.format(baseDir)]
    #
    # writeListToDestination('{}/arifOTU/mothurBatch'.format(baseDir), mothurBatch)
    #
    # completedProcess = subprocess.run(['mothur', '{}/arifOTU/mothurBatch'.format(baseDir)])
    ##################################################################
    # So I have looked at the data that this output and all of the sequences clustered together at < 0.02 and the

    # cutoff was changed to 0.0125.
    # this means that the plot can be very simple indeed. Just a single bar needs to be plotted for each sample
    # C3.
    plotInfo = [100 for dir in listOfDirs]
    plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1, 1, 1)
    ax.set_ylim(bottom=0, top=100)
    # 8 = plot, plot, space, plot, plot, space, plot, plot
    ind = range(len(listOfDirs))
    width = 1
    bottom = [0 for dir in listOfDirs]
    ax.bar(ind, plotInfo, width, bottom=bottom, color='#E8AB52')

def QCOfEdSeqs():
    ''' This function performs the Mothur QC on the sequences on a sample by sample basis, producing a redundant
    i.e. deuniqued, .fasta file at the end which is written out. The MED analysis and clustering can then be run
    on these files'''

    baseDir = '/home/humebc/projects/SymPortalMS/edData/'

    listOfFiles = []
    for (dirpath, dirnames, filenames) in os.walk(baseDir):
        listOfFiles.extend(filenames)
        break
    listOfFiles = sorted([file for file in listOfFiles if '.fasta' in file])

    # # For some reason the Ed files have the rev primer removed but not the Fwd primer
    # # I will add the reverse primer on so that it can be removed by the seq.pcr in mothur
    # # the degenerate nucleotides was causing a problem in mothur so I have replced them.
    # for file in listOfFiles:
    #     print(file)
    #     localFasta = []
    #     with open('{}/{}'.format(baseDir, file)) as f:
    #         localFasta = f.read().splitlines()
    #     revPrimer = 'CGGGTTCWCTTGTYTGACTTCATGC'
    #     for i in range(len(localFasta)):
    #         if localFasta[i][0] == '>':
    #             localFasta[i+1] = localFasta[i+1].replace('CGGGTTCACTTGTCTGACTTCATGC', 'GCATGAAGTCAGACAAGTGAACCCG')
    #     writeListToDestination('{}/{}'.format(baseDir, file), localFasta)

    pear = 'apple'
    for file in listOfFiles:

        rootName = file.replace('.fasta', '') + '.CladeC'

        os.makedirs('{}/{}'.format(baseDir, file.replace('.fasta', '')), exist_ok=True)
        outputDir = '{}{}/'.format(baseDir, file.replace('.fasta', ''))
        os.chdir('{}'.format(outputDir))

        originalFasta = readDefinedFileToList('{}{}'.format(baseDir, file))
        # TODO, check to see if the  C fasta already exisits.
        # Extract a new fasta from the originalFasta that is only the C seqs and write this out for MOTHUR to use
        extractCSeqs(originalFasta, outputDir)
        primerFwdSeq = 'GAATTGCAGAACTCCGTGAACC'  # Written 5'-->3'
        primerRevSeq = 'CGGGTTCWCTTGTYTGACTTCATGC'  # Written 5'-->3'

        oligoFile = [
            r'#SYM_VAR_5.8S2',
            'forward\t{0}'.format(primerFwdSeq),
            r'#SYM_VAR_REV',
            'reverse\t{0}'.format(primerRevSeq)
        ]

        writeListToDestination('{}primers.oligos'.format(outputDir), oligoFile)

        mBatchFile = [
            r'set.dir(input={0})'.format(outputDir),
            r'set.dir(output={0})'.format(outputDir),
            r'summary.seqs(fasta={}{}.fasta)'.format(outputDir, rootName),
            # NB that the maxabmbig has to be set to 2 as I have 2 degenerate nucleotides in the revprimer seqs
            r'screen.seqs(fasta={}{}.fasta, maxambig=0, maxhomop=5)'.format(
                outputDir, rootName),
            r'summary.seqs(fasta={0}{1}.good.fasta)'.format(outputDir, rootName),
            r'unique.seqs(fasta={0}{1}.good.fasta)'.format(outputDir, rootName),
            r'summary.seqs(fasta={0}{1}.good.unique.fasta, name={0}{1}.good.names)'.format(
                outputDir, rootName),
            r'split.abund(cutoff=2, fasta={0}{1}.good.unique.fasta, name={0}{1}.good.names)'.format(
                outputDir, rootName),
            r'summary.seqs(fasta={0}{1}.good.unique.abund.fasta, name={0}{1}.good.abund.names)'.format(
                outputDir, rootName),
            r'pcr.seqs(fasta={0}{1}.good.unique.abund.fasta, name={0}{1}.good.abund.names, oligos={0}primers.oligos, pdiffs=2, rdiffs=2)'.format(
                outputDir, rootName),
            r'summary.seqs(fasta={0}{1}.good.unique.abund.pcr.fasta, name={0}{1}.good.abund.pcr.names)'.format(
                outputDir, rootName),
            r'unique.seqs(fasta={0}{1}.good.unique.abund.pcr.fasta, name={0}{1}.good.abund.pcr.names)'.format(
                outputDir, rootName),
            r'summary.seqs(fasta={0}{1}.good.unique.abund.pcr.unique.fasta, name={0}{1}.good.unique.abund.pcr.names)'.format(
                outputDir, rootName),
            r'screen.seqs(fasta={0}{1}.good.unique.abund.pcr.fasta, name={0}{1}.good.abund.pcr.names,  minlength=184, maxlength=310)'.format(
                outputDir, rootName),
            r'summary.seqs(fasta={0}{1}.good.unique.abund.pcr.good.fasta, name={0}{1}.good.abund.pcr.good.names)'.format(
                outputDir, rootName),
            r'unique.seqs(fasta={0}{1}.good.unique.abund.pcr.good.fasta, name={0}{1}.good.abund.pcr.good.names)'.format(
                outputDir, rootName),
            r'summary.seqs(fasta={0}{1}.good.unique.abund.pcr.good.unique.fasta, name={0}{1}.good.unique.abund.pcr.good.names)'.format(
                outputDir, rootName),
            r'deunique.seqs(fasta={0}{1}.good.unique.abund.pcr.good.unique.fasta, name={0}{1}.good.unique.abund.pcr.good.names)'.format(outputDir, rootName)
        ]
        pathToBatchFile = '{}mBatchFile{}'.format(outputDir, rootName)
        writeListToDestination(pathToBatchFile, mBatchFile)
        completedProcess = subprocess.run(['mothur', r'{0}'.format(pathToBatchFile)])


    # here we have completed all of the inital QC for the analysisprocessed all of the .fastas and we have written out the fasta.
    # I think we should now do the 3% clustering and the MED analyses in seperate functions as this is where the
    # three figures start to separate


    return

def MEDQCedSeqs():
    baseDir = '/home/humebc/projects/SymPortalMS/edData'
    listOfDirs = []
    for (dirpath, dirnames, filenames) in os.walk(baseDir):
        listOfDirs.extend(dirnames)
        break
    listOfDirs = sorted([dir for dir in listOfDirs if dir != 'transFiles'])

    # we want to perform MED on each of the samples
    for dir in listOfDirs:
        pathToFile = '{}/{}/{}.CladeC.good.unique.abund.pcr.good.unique.redundant.fasta'.format(baseDir, dir, dir)
        # Sanity check at debug to how sample information is extracted
        # completedProcess = subprocess.run([r'o-get-sample-info-from-fasta', r'{0}'.format(pathToFile)])
        completedProcess = subprocess.run([r'o-pad-with-gaps', r'{0}'.format(pathToFile)])
        # Decompose
        pathToFile = '{}/{}/{}.CladeC.good.unique.abund.pcr.good.unique.redundant.fasta-PADDED-WITH-GAPS'.format(baseDir, dir, dir)
        outputDir = '{}/{}/'.format(baseDir, dir)
        completedProcess = subprocess.run(
            [r'decompose', '--skip-gexf-files', '--skip-gen-figures', '--skip-gen-html', '--skip-check-input', '-o',
             outputDir, pathToFile])

    # Now work out what datastructure we want to hold the plot info in
    # I think it would be good to be able to use the figure2rawseq code for plotting.
    # for this we only need to create MEDed fastas from the MED outputs
    # then we can feed these straight into plotting
    debaug = 'asdf'

    # Do the same as we did for the last plot
    # first go through each of the files and add the rel abund of each seq to a total seq dict
    # then sort this dict to get an order of the sequences
    # then go back through each of the files again but this time add the info to a plot info dict
    # that has the seqs as keys in order and 0 for each of the samples


    # for each file
    # read in the nodes fasta and make dictionary

def consolidateSeqs(abundDict, sigSeqsDict, file):
    # Because the Ed sequences are longer than the ITS2 seqs that are significant, a single sig seq
    # can match multiple Ed seqs. We will need to control for this. Preferably by clustering all of the
    # sequences that match an Ed sequence and assigning a representative sequence with a total abundance of the
    # sequences it represents
    print('Consolidating Seqs for {}'.format(file))
    for sigSeq in sigSeqsDict.keys():
        tempList = []
        total = 0
        for edSeq, abund in abundDict.items():
            # Get a list of the edSeqs that match the sigSeq
            if sigSeqsDict[sigSeq] in edSeq:
                # then this is a match
                tempList.append((edSeq, abund))
                total += abund
        if tempList:
            # Here we have a list of each edSeq that represents the
            # This should give us a sorted list according to the abundance
            sorted(tempList, key=lambda x: (-x[1], x[0]))

            # # Debug #
            # maxAbund = abundDict[tempList[0][0]]
            # ########
            representativeSeq = tempList[0][0]
            representativeAbund = total
            # Delete all the other seqs from the abundDict
            for seqTup in tempList:
                if seqTup[0] != representativeSeq:
                    del abundDict[seqTup[0]]
            # Change the abundance of the now representativeSeq
            abundDict[representativeSeq] = representativeAbund

    return abundDict

def checkIfSeqSig(seqtocheck, sigseqdict):
    '''The purpose of this is simply to take one of Eds sequences, reverse complement it and then check if
    any of the sig seqences fit into it. If they do then this is one of the significant sequences and
    this should return true
    Eddo! why always with the reverse sequencing ha ha ha ha
    Just to keep me on my toes right :)
    '''

    # revCompSeq = reverse_complement(seqtocheck)

    for seqName, sigSeq in sigseqdict.items():
        # if sigSeq in revCompSeq:
        if sigSeq in seqtocheck or seqtocheck in sigSeq:
            return (True, seqName)
    return (False, seqName)

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def blastForC(fasta):
    ''' This method will take a fasta file in the form of a list
    blast it against the ITS2 clade database
    and create a new fasta that only contains clade C sequences
    '''

    # Create a unique names and fasta set to speed up
    # create a defult dict to count the seqs and get abundances
    fastaDD = defaultdict(int)
    for i in range(len(fasta)):
        if fasta[i][0] == '>':
            fastaDD[fasta[i + 1]] += 1

    # now write the unique dict as a fasta
    # blast it and keep only the seqs that are C in a new dict
    # output this dict, we can work with it
    counter = 0
    newFasta = []
    converterDict = {}
    for key in fastaDD.keys():
        newFasta.extend(['>seq{}'.format(counter), key])
        # Here we need a dict to convert the counter seq to the sequence
        converterDict[key] = 'seq{}'.format(counter)
        counter += 1

    convertDict = createDictFromFasta(newFasta)


    writeListToDestination('{}/blastIn.fasta'.format('/home/humebc/projects/SymPortalMS'), newFasta)
    # Print out the .ncbirc file in the cwd that shows blast were to find the DBs in question
    ncbirc = ['[BLAST]', 'BLASTDB={}'.format('/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload')]
    writeListToDestination('{}/.ncbirc'.format('/home/humebc/projects/SymPortalMS'), ncbirc)

    outputpath = '{}/symClade.out'.format('/home/humebc/projects/SymPortalMS')
    outputFmt = '6 qseqid sseqid'
    inputPath = '{}/blastIn.fasta'.format('/home/humebc/projects/SymPortalMS')
    completedProcess = subprocess.run(
        ['blastn', '-out', outputpath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'ITS2Clade.fas',
         '-max_target_seqs', '1', '-num_threads', '10'])

    done = []
    # Read in blast output
    blastOutput = readDefinedFileToList('{}/symClade.out'.format('/home/humebc/projects/SymPortalMS'))

    # create dict from the blastOutput
    blastDict = {line.split('\t')[0]: line.split('\t')[1][5] for line in blastOutput}

    outputDict = {}
    for key in fastaDD.keys():
        if converterDict[key] in blastDict.keys():
            if blastDict[converterDict[key]] == 'C':
                outputDict[key] = fastaDD[key]
        else:
            continue




    return outputDict

def extractCSeqs(fasta, outputDir):
    ''' This method will take a fasta file in the form of a list
        blast it against the ITS2 clade database
        and create a new fasta that only contains clade C sequences
        '''
    rootName = outputDir.split('/')[-2]

    fastaOutName = rootName + '.CladeC.fasta'

    # Check to see if the C fasta for this sample already exists
    try:
        # readDefinedFileToList('bob')
        readDefinedFileToList('{}{}'.format(outputDir, fastaOutName))
        return '{}{}'.format(outputDir, fastaOutName)
    except:
        # Create a unique names and fasta set to speed up
        # create a defult dict to count the seqs and get abundances
        fastaDD = defaultdict(int)
        for i in range(len(fasta)):
            if fasta[i][0] == '>':
                fastaDD[fasta[i + 1]] += 1

        # now write the unique dict as a fasta
        # blast it and keep only the seqs that are C in a new dict
        # output this dict, we can work with it
        counter = 0
        newFasta = []
        converterDict = {}
        for key in fastaDD.keys():
            newFasta.extend(['>seq{}'.format(counter), key])
            # Here we need a dict to convert the counter seq to the sequence
            converterDict[key] = 'seq{}'.format(counter)
            counter += 1

        convertDict = createDictFromFasta(newFasta)

        writeListToDestination('{}blastIn.fasta'.format(outputDir), newFasta)
        # Print out the .ncbirc file in the cwd that shows blast were to find the DBs in question
        ncbirc = ['[BLAST]', 'BLASTDB={}'.format('/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload')]
        writeListToDestination('{}.ncbirc'.format(outputDir), ncbirc)

        outputpath = '{}symClade.out'.format(outputDir)
        outputFmt = '6 qseqid sseqid'
        inputPath = '{}blastIn.fasta'.format(outputDir)
        completedProcess = subprocess.run(
            ['blastn', '-out', outputpath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'ITS2Clade.fas',
             '-max_target_seqs', '1', '-num_threads', '10'])

        done = []
        # Read in blast output
        blastOutput = readDefinedFileToList('{}symClade.out'.format(outputDir))

        # create dict from the blastOutput
        blastDict = {line.split('\t')[0]: line.split('\t')[1][5] for line in blastOutput}
        outputDict = {}
        for key in fastaDD.keys():
            if converterDict[key] in blastDict.keys():
                if blastDict[converterDict[key]] == 'C':
                    outputDict[key] = fastaDD[key]
            else:
                continue

        # Here write out a redundant cladeC fasta
        # name should be e.g. A01.CladeC.fasta
        fastaOut = []
        seqCount = 0
        for key, abund in outputDict.items():
            for i in range(abund):
                fastaOut.extend(['>seq{}'.format(seqCount), key])
                seqCount += 1

        writeListToDestination('{}{}'.format(outputDir, fastaOutName), fastaOut)
        return '{}{}'.format(outputDir, fastaOutName)
##########################################

################# FIGURE 3 ##############

def fig3_fixed_types():
    ''' The aim of this function will be to generate figure three. It will have a type class and a sample class.
    I will aim to butcher the dynamic_types function as much as possible to save time. Each sample will be generated
    with the same amount of the same types each time. Samples will be a mix of clade A, C and D types but mainly, C.
    We will actually need to make about 8 differnt plots for this figure.

    1) the raw dataplot
    2, 3, 4) the cladaly separated raw dataplots
    5, 6, 7 the post QC and MED plots of reach of the clades
    8, 9, 10) the the predicted ITS2 type profiles for each clade
    11) the reconstructed output which is the type profiles for each clade combined in the single sample plot'''

    sample_list = generate_samples()

    ###### GET ORDER OF SEQUENCES FOR PLOTTING ######
    # Go through each of the samples to get a master abund dict of the sequences to get an order
    seqAbundDict = defaultdict(float)
    for sample in sample_list:
        for v_type_profile in sample.list_of_ITS2_type_objects:
            for seq_name, seq_abs_raw_abund in v_type_profile.raw_seq_abs_abundances.items():
                seqAbundDict[seq_name] += seq_abs_raw_abund

    # Now sort the dict to get list of sequences in order of abundance
    sortedSeqs = [k for k in sorted(seqAbundDict, key=seqAbundDict.get, reverse=True)]
    ###################################################

    ######### GENERATE AND POPULATE plot_info_dict ########
    # Now populate a plot_info_dict for doing the plotting
    plot_info_dict = {seq: [0 for sample in sample_list] for seq in sortedSeqs}
    for sample in sample_list:
        for seq in plot_info_dict.keys():
            for v_type_profile in sample.list_of_ITS2_type_objects:
            # Here we are populating a given samples worth of the plot_info_dict
            # we only want to populate according to what is in the sample
             # for each of the seqs we need to populate info for
                if seq in v_type_profile.raw_seq_abs_abundances:
                    plot_info_dict[seq][sample_list.index(sample)] = v_type_profile.raw_seq_abs_abundances[seq]


    ######################################################

    ######### GENERATE COLOUR PALETTES ########
    # Now we have the info ready to be plotted
    # Colour palettes
    clade_palette = {'C': '#F9ED32', 'A': '#4BAADF', 'D': '#6D6E71'}
    c_palette = ['#F9ED32', '#BC802B', '#836126', '#F89B58', '#FBAC1F', '#F58E82', '#F15F48', '#EE2F24', '#A52022',
                 '#F09B20',
                 '#F9F292', '#F9EE49', '#E7E621', '#AAB036']
    a_palette = ['#4BAADF', '#A8D598', '#9AD5C5', '#75CBC8', '#7CC243', '#21B24B', '#0A7B3E', '#5CC7D4',
                 '#00A6C5', '#0066AB', '#236FB7', '#2A3588', '#645EA9', '#895DA6']

    colourPalette = ['#E8AB52', '#5FA96E', '#9680D3', '#D8C542', '#D47070', '#59A6D4', '#D76FBC', '#F6724C',
                     '#5FC2B6', '#80D1E4', '#6BD78A', '#B7D456']
    greyPalette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']
    grey_palette_two = ['#4C4C4C', '#666666', '#7F7F7F', '#999999', '#B2B2B2', '#CCCCCC', '#E5E5E5']
    ###########################################

    ######### GENERATE CLADAL PLOT ########
    # First plot the clades with the simple colouring
    # Set up plot
    bob = plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1, 1, 1)
    ax.set_ylim(bottom=0, top=1000)
    ind = range(len(sample_list))
    width = 1

    # Now plot
    # We have to plot in order of clade A, C then D first
    bottom = [0 for sample in sample_list]
    colourCounter = 0
    totSeqs = len(sortedSeqs)
    for clade in ['A', 'C', 'D']:
        for seq in sortedSeqs:
            if clade in seq:
                print('Printing {} out of {} bars'.format(colourCounter, totSeqs))
                ax.bar(ind, plot_info_dict[seq], width, bottom=bottom, color=grey_palette_two[colourCounter % 7])
                bottom = [L + M for L, M in zip(bottom, plot_info_dict[seq])]
                colourCounter += 1
                #######################################

    # Here we should have the first set of bars plotted
    print(str(os.getcwd()))
    # plt.show()
    plt.savefig('fig3_raw.svg', format="svg")
    plt.close(bob)
    debugToHere = 'asdf'



    ######### CLADAL SEPARATED RAW PLOTS ##############
    # Now do the next three plots which will be the cladaly separated raw
    for clade in ['A', 'C', 'D']:

        # Now sort the dict to get list of sequences in order of abundance
        sortedSeqs_clade = [seq_name for seq_name in sortedSeqs if clade in seq_name]
        ###################################################

        ######### GENERATE AND POPULATE plot_info_dict ########
        # Now populate a plot_info_dict for doing the plotting
        plot_info_dict = {seq: [0 for sample in sample_list] for seq in sortedSeqs_clade}
        for sample in sample_list:
            # Here we need to work out the total rel abundances of types of the clade_in_question
            total_raw_seqs_in_samples_type = sum([type_abundance for type_name, type_abundance in sample.type_abundance_dict.items() if clade in type_name])
            for seq in plot_info_dict.keys():
                for v_type_profile in sample.list_of_ITS2_type_objects:
                    # Here we are populating a given samples worth of the plot_info_dict
                    # we want to be sure to still normalise to 1 or 1000 in this case
                    if seq in v_type_profile.raw_seq_abs_abundances:
                        plot_info_dict[seq][sample_list.index(sample)] = sample.sequencing_depth * ((v_type_profile.raw_seq_abs_abundances[seq] / sample.sequencing_depth) / total_raw_seqs_in_samples_type)

        # Here we should be ready to plot the info for the second set of plots
                    ######## GENERATE CLADAL PLOT ########
        # First plot the clades with the simple colouring
        # Set up plot
        bob = plt.figure(num=1, figsize=(8, 10), dpi=80)
        ax = plt.subplot(1, 1, 1)
        ax.set_ylim(bottom=0, top=1000)
        ind = range(len(sample_list))
        width = 1

        # Now plot
        # We have to plot in order of clade A, C then D first
        bottom = [0 for sample in sample_list]
        colourCounter = 0
        totSeqs = len(sortedSeqs_clade)
        for seq in sortedSeqs_clade:
            print('Printing {} out of {} bars'.format(colourCounter, totSeqs))
            ax.bar(ind, plot_info_dict[seq], width, bottom=bottom,
                   color=grey_palette_two[colourCounter % 7])
            bottom = [L + M for L, M in zip(bottom, plot_info_dict[seq])]
            colourCounter += 1

        # Here we should have the first set of bars plotted
        print(str(os.getcwd()))
        # plt.show()
        plt.savefig('fig3_raw_{}.svg'.format(clade), format="svg")
        debugToHere = 'asdf'
        plt.close(bob)
    #######################################

    ###### DIV PLOTS ############
    # These plots ae in theory exactly the same plots as the cladally separated plots
    # except that they have gone through MED and they now have the DIVs coloured
    # while the non-DIV sequences remain grey
    # We can re-use a lot of the code above.
    # we should first get a list of sequences by order of abundance from the
    # post-qc and MEDed sequence dictionaries

    DIV_palette = {'A1': '#332288', 'A1aw':'#6699CC', 'A1av':'#88CCEE', 'A1ax':'#44AA99', 'C15':'#64C204', 'C15x':'#117733', 'C3':'#999933', 'C3a':'#FEC44F', 'C3ab':'#EC7014', 'C3c':'#993404', 'C3cc':'#DDCC77', 'C1':'#661100', 'C1ab':'#FF6666', 'C1e':'#CC6677', 'D1':'#AA4466', 'D4':'#882255', 'D6':'#AA4499'}

    ###### GET ORDER OF SEQUENCES FOR PLOTTING ######
    # Go through each of the samples to get a master abund dict of the sequences to get an order
    seqAbundDict = defaultdict(float)
    for sample in sample_list:
        for v_type_profile in sample.list_of_ITS2_type_objects:
            for seq_name, seq_abs_qc_MED_abund in v_type_profile.post_qc_abs_aundances.items():
                seqAbundDict[seq_name] += seq_abs_qc_MED_abund

    # Now sort the dict to get list of sequences in order of abundance
    sortedSeqs = [k for k in sorted(seqAbundDict, key=seqAbundDict.get, reverse=True)]

    # Now do the next three plots which will be the cladaly separated raw
    for clade in ['A', 'C', 'D']:

        # Now sort the dict to get list of sequences in order of abundance
        sortedSeqs_clade = [seq_name for seq_name in sortedSeqs if clade in seq_name]
        ###################################################

        ######### GENERATE AND POPULATE plot_info_dict ########
        # Now populate a plot_info_dict for doing the plotting
        plot_info_dict = {seq: [0 for sample in sample_list] for seq in sortedSeqs_clade}
        for sample in sample_list:
            total_post_qc_MED_seqs_in_samples_type = sum([type_abund for type_name, type_abund in sample.type_abundance_dict.items() if clade in type_name])
            for seq in plot_info_dict.keys():
                for v_type_profile in sample.list_of_ITS2_type_objects:
                # Here we are populating a given samples worth of the plot_info_dict
                # we want to be sure to still normalise to 1 or 1000 in this case
                    if seq in v_type_profile.post_qc_abs_aundances:
                        plot_info_dict[seq][sample_list.index(sample)] = sample.sequencing_depth * ((v_type_profile.post_qc_abs_aundances[seq] / sample.sequencing_depth) / total_post_qc_MED_seqs_in_samples_type)

                    # Here we should be ready to plot the info for the second set of plots
                    ######### GENERATE CLADAL PLOT ########
        # First plot the clades with the simple colouring
        # Set up plot
        bob = plt.figure(num=1, figsize=(8, 10), dpi=80)
        ax = plt.subplot(1, 1, 1)
        ax.set_ylim(bottom=0, top=1000)
        ind = range(len(sample_list))
        width = 1

        # Now plot
        # We have to plot in order of clade A, C then D first
        bottom = [0 for sample in sample_list]
        colourCounter = 0
        totSeqs = len(sortedSeqs_clade)
        for seq in sortedSeqs_clade:
            print('Printing {} out of {} bars'.format(colourCounter, totSeqs))
            # If non-Div use grey palette else use the DIV palette
            if 'Unk' in seq:
                ax.bar(ind, plot_info_dict[seq], width, bottom=bottom,
                       color=grey_palette_two[colourCounter % 7])

                colourCounter += 1
            else:
                ax.bar(ind, plot_info_dict[seq], width, bottom=bottom,
                       color=DIV_palette[seq])
            bottom = [L + M for L, M in zip(bottom, plot_info_dict[seq])]
        # Here we should have the first set of bars plotted
        print(str(os.getcwd()))

        # plt.show()
        plt.savefig('fig3_MED_post_qc_{}.svg'.format(clade), format="svg")
        debugToHere = 'asdf'
        plt.close(bob)

    ###################################################

    ######## Clade separated type profile plots ########
    # Here we will want to go clade by clade again
    # then sample by sample

    type_palette = {'A1-A1aw': '#332288', 'A1-A1av-A1ax': '#6699CC', 'C15-C15x': '#88CCEE',
                    'C3-C3a-C3ab': '#44AA99', 'C3-C3c-C3cc': '#64C204',
                   'C1-C1ab-C1e': '#117733', 'D1-D4-D6': '#999933'}

    for clade in ['A', 'C', 'D']:
        # Get a list of all of the type names of the given clade
        # so that we can make a plot_info_dict
        types_rel_abundances_dict = defaultdict(float)
        for sample in sample_list:
            for type_name, type_rel_abund in sample.type_abundance_dict.items():
                if clade in type_name:
                    types_rel_abundances_dict[type_name] += type_rel_abund

        sorted_types = [k for k in sorted(types_rel_abundances_dict, key=types_rel_abundances_dict.get, reverse=True)]
        soreted_seqs_types = [type_name for type_name in sorted_types if clade in type_name]
        # Here we have the names of the types that are found in all of the samples
        # now go back through and populate a plot_info_dict
        plot_info_dict_type = {type_name: [0 for samples in sample_list] for type_name in soreted_seqs_types}
        for sample in sample_list:
            # Here need to get the total rel abundances of type of the given clade found in the sample
            rel_abund_total_for_types_of_the_given_clade_within_this_sample = 0
            for type_profile in sample.list_of_ITS2_type_objects:
                if type_profile.clade == clade:
                    rel_abund_total_for_types_of_the_given_clade_within_this_sample += sample.type_abundance_dict[type_profile.name]

            for type_profile_name in plot_info_dict_type.keys():
                if type_profile_name in sample.type_abundance_dict:
                    # we must remember to normalise the type abundances here.We know the total of the rel abundances
                    # of the types of the clade in quesiton in the sequecnes in question so we will use this
                    plot_info_dict_type[type_profile_name][sample_list.index(sample)] = sample.sequencing_depth * (sample.type_abundance_dict[type_profile_name] / rel_abund_total_for_types_of_the_given_clade_within_this_sample)

        # Here we have the plot_info_dict populated with some very simple data which is the type abundance
        bob = plt.figure(num=1, figsize=(8, 10), dpi=80)
        ax = plt.subplot(1, 1, 1)
        ax.set_ylim(bottom=0, top=1000)
        ind = range(len(sample_list))
        width = 1

        # Now plot
        # We have to plot in order of clade A, C then D first
        bottom = [0 for sample in sample_list]


        for type_name in soreted_seqs_types:

            # If non-Div use grey palette else use the DIV palette


            ax.bar(ind, plot_info_dict_type[type_name], width, bottom=bottom, color=type_palette[type_name])
            bottom = [L + M for L, M in zip(bottom, plot_info_dict_type[type_name])]
        # Here we should have the first set of bars plotted
        print(str(os.getcwd()))

        # plt.show()
        plt.savefig('fig3_ITS2_type_profile_{}.svg'.format(clade), format="svg")
        debugToHere = 'asdf'
        plt.close(bob)

    ######## Final plot #######
    # this plot should be the re-combined type plots.
    # we should be able to borrow heavily from above and simply not separate by clade

        ######## Clade separated type profile plots ########
        # Here we will want to go clade by clade again
        # then sample by sample

    type_palette = {'A1-A1aw': '#332288', 'A1-A1av-A1ax': '#6699CC', 'C15-C15x': '#88CCEE',
                    'C3-C3a-C3ab': '#44AA99', 'C3-C3c-C3cc': '#64C204',
                    'C1-C1ab-C1e': '#117733', 'D1-D4-D6': '#999933'}


        # Get a list of all of the type names of the given clade
        # so that we can make a plot_info_dict
    types_rel_abundances_dict = defaultdict(float)
    for sample in sample_list:
        for type_name, type_rel_abund in sample.type_abundance_dict.items():
            types_rel_abundances_dict[type_name] += type_rel_abund

    sorted_types = [k for k in
                    sorted(types_rel_abundances_dict, key=types_rel_abundances_dict.get, reverse=True)]

    # Here we have the names of the types that are found in all of the samples
    # now go back through and populate a plot_info_dict
    plot_info_dict_type = {type_name: [0 for samples in sample_list] for type_name in sorted_types}
    for sample in sample_list:
        for type_profile_name in plot_info_dict_type.keys():
            if type_profile_name in sample.type_abundance_dict:
                # we must remember to normalise the type abundances here.We know the total of the rel abundances
                # of the types of the clade in quesiton in the sequecnes in question so we will use this
                plot_info_dict_type[type_profile_name][sample_list.index(sample)] = \
                    sample.sequencing_depth * sample.type_abundance_dict[type_profile_name]

    # Here we have the plot_info_dict populated with some very simple data which is the type abundance
    bob = plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1, 1, 1)
    ax.set_ylim(bottom=0, top=1000)
    ind = range(len(sample_list))
    width = 1

    # We should plot this in the same order as the original raw data was plotted
    # in terms of clade so that it has exactly the same shape.
    bottom = [0 for sample in sample_list]
    for clade in ['A', 'C', 'D']:

        # Now plot
        # We have to plot in order of clade A, C then D first


        for type_name in sorted_types:
            # If non-Div use grey palette else use the DIV palette
            if clade in type_name:

                ax.bar(ind, plot_info_dict_type[type_name], width, bottom=bottom, color=type_palette[type_name])
                bottom = [L + M for L, M in zip(bottom, plot_info_dict_type[type_name])]
        # Here we should have the first set of bars plotted
    print(str(os.getcwd()))

    # plt.show()
    plt.savefig('fig3_ITS2_type_profile_final.svg', format="svg")
    debugToHere = 'asdf'
    plt.close(bob)



def fig3_dynamic_types():
    ####################### GENERATE TYPES ######################
    ''' 110118 I got a litle carried away with this code and was making a truely virtual and somewhat randomly
    generated set of sequences which was, in all frankness, over the top. I am now going to work on the functin
    above, fig3_fixed_types. Please see its comments for what it is about'''
    # First generate the types that we'll be using:
    # 2 A types and 9 C types



    sample_list = generateSamples()

    ###### GET ORDER OF SEQUENCES FOR PLOTTING ######
    # Go through each of the samples to get a master abund dict of the sequences to get an order
    seqAbundDict = defaultdict(float)
    for sample in sample_list:
        for seq, abund in sample.seqCounts.items():
            seqAbundDict[seq] += abund


    # Now sort the dict to get list of sequences in order of abundance
    sortedSeqs = [k for k in sorted(seqAbundDict, key=seqAbundDict.get, reverse=True)]
    ###################################################

    ######### GENERATE AND POPULATE plot_info_dict ########
    # Now populate a plot_info_dict for doing the plotting
    plot_info_dict = {seq: [0 for sample in sample_list] for seq in sortedSeqs}
    for sample in sample_list:
        for seq, abund in sample.seqCounts.items():
            plot_info_dict[seq][sample_list.index(sample)] = abund

    ######################################################

    ######### GENERATE COLOUR PALETTES ########
    # Now we have the info ready to be plotted
    # Colour palettes
    clade_palette = {'C': '#F9ED32', 'A': '#4BAADF', 'D': '#6D6E71'}
    c_palette = ['#F9ED32', '#BC802B', '#836126', '#F89B58', '#FBAC1F', '#F58E82', '#F15F48', '#EE2F24', '#A52022',   '#F09B20',
                  '#F9F292', '#F9EE49', '#E7E621', '#AAB036']
    a_palette = ['#4BAADF','#A8D598', '#9AD5C5', '#75CBC8',   '#7CC243', '#21B24B', '#0A7B3E',   '#5CC7D4',
                 '#00A6C5', '#0066AB', '#236FB7', '#2A3588', '#645EA9', '#895DA6']


    colourPalette = ['#E8AB52', '#5FA96E', '#9680D3', '#D8C542', '#D47070', '#59A6D4', '#D76FBC', '#F6724C',
                     '#5FC2B6', '#80D1E4', '#6BD78A', '#B7D456']
    greyPalette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

    ###########################################

    ######### GENERATE CLADAL PLOT ########
    # First plot the clades with the simple colouring
    # Set up plot
    plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1, 1, 1)
    ax.set_ylim(bottom=0, top=1000)
    ind = range(len(sample_list))
    width = 1


    # Now plot
    # We have to plot in order of clade A, C then D first
    bottom = [0 for sample in sample_list]
    colourCounter = 0
    totSeqs = len(sortedSeqs)
    for clade in ['A', 'C', 'D']:
        for seq in sortedSeqs:
            if clade in seq:
                print('Printing {} out of {} bars'.format(colourCounter, totSeqs))
                ax.bar(ind, plot_info_dict[seq], width, bottom=bottom, color=greyPalette[colourCounter % 6])
                bottom = [L + M for L, M in zip(bottom, plot_info_dict[seq])]
                colourCounter += 1
        #######################################

    # Here we should have the first set of bars plotted
    print(str(os.getcwd()))
    plt.savefig('fig3_raw.svg', format="svg")
    debugToHere = 'asdf'


    #TODO at this point we need to simulate putting the samples through QC and


    ############# GENERATE SUBCLADAL PLOTS ##############
    # now we need to plot a bar for each sample that is only for clade A and clade C
    # clade D did not appear at sufficient abundance, per our design, to be considered.
    # The clade D sequences will be added back in at the end
    # we will plot two separate plots, one for each clade
    # we will need to produce two separate plot_info_dicts
    # For each of the samples, if there are sequences, they should be normalised to 1
    # i.e. we should be plotting proportions
    # DO C
    # First get list of sequences that are clade C and sort according to abundance
    seqAbundDict_C = defaultdict(float)
    for sample in sample_list:
        for seq, abund in sample.seqCounts.items():
            if 'C' in seq:
                seqAbundDict_C[seq] += abund

    # Now sort the dict to get list of sequences in order of abundance
    sortedSeqs_C = [k for k in sorted(seqAbundDict_C, key=seqAbundDict_C.get, reverse=True)]

    # Now generte plot_info_dict
    plot_info_dict_C = {seq: [0 for sample in sample_list] for seq in sortedSeqs_C}
    for sample in sample_list:
        tot_C = sum([sample.seqCounts[seq] for seq in sample.seqCounts.keys() if 'C' in seq])
        for seq, abund in sample.seqCounts.items():
            if 'C' in seq:
                plot_info_dict_C[seq][sample_list.index(sample)] = (abund/tot_C) * 100

    # DO A
    # First get list of sequences that are clade C and sort according to abundance
    seqAbundDict_A = defaultdict(float)
    for sample in sample_list:
        for seq, abund in sample.seqCounts.items():
            if 'A' in seq:
                seqAbundDict_A[seq] += abund

    # Now sort the dict to get list of sequences in order of abundance
    sortedSeqs_A = [k for k in sorted(seqAbundDict_A, key=seqAbundDict_A.get, reverse=True)]

    # Now generte plot_info_dict
    plot_info_dict_A = {seq: [0 for sample in sample_list] for seq in sortedSeqs_A}
    for sample in sample_list:
        tot_A = sum([sample.seqCounts[seq] for seq in sample.seqCounts.keys() if 'A' in seq])
        for seq, abund in sample.seqCounts.items():
            if 'A' in seq:
                plot_info_dict_A[seq][sample_list.index(sample)] = (abund/tot_A) * 100

    # Now we can start plotting up the two subcladal plotsinfodicts
    # Plot C
    plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1, 3, 2)
    ax.set_ylim(bottom=0, top=100)
    ind = range(len(sample_list))
    width = 1

    # Now plot
    # We have to plot in order of clade A, C then D first
    bottom = [0 for sample in sample_list]
    colourCounter = 0
    for seq in sortedSeqs_C:
        ax.bar(ind, plot_info_dict_C[seq], width, bottom=bottom, color=c_palette[colourCounter % len(c_palette)])
        bottom = [L + M for L, M in zip(bottom, plot_info_dict_C[seq])]
        colourCounter += 1

    # Plot A
    plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1, 3, 3)
    ax.set_ylim(bottom=0, top=100)
    ind = range(len(sample_list))
    width = 1

    # Now plot
    # We have to plot in order of clade A, C then D first
    bottom = [0 for sample in sample_list]
    colourCounter = 0
    for seq in sortedSeqs_A:
        ax.bar(ind, plot_info_dict_A[seq], width, bottom=bottom, color=a_palette[colourCounter % len(a_palette)])
        bottom = [L + M for L, M in zip(bottom, plot_info_dict_A[seq])]
        colourCounter += 1

    adsf = 'asdf'
    #####################################################




def generate_samples():



    #TODO make it so that if types are added, then they will appear in the sample. I think we need to minimis the variation a bit
    # To do this we can insist that the generated abundance is within a certain range, i.e. above 0.
    # Generate Samples and add to sample_list
    sample_list = []
    # One A1-A1aw, One C15-C15x
    vSamp1 = virtual_sample(
        sample_name='vSamp1',
        type_abundance_dict={'A1-A1aw' : 0.25, 'C3-C3a-C3ab': 0.5, 'D1-D4-D6' : 0.25},
        sequencing_depth=1000)
    # One A1-A1aw, One C15-C15x
    vSamp2 = virtual_sample(
        sample_name='vSamp2',
        type_abundance_dict={'C3-C3a-C3ab': 0.85, 'D1-D4-D6' : 0.15},
        sequencing_depth=1000)
    # One A1-A1aw, One C15-C15x
    vSamp3 = virtual_sample(
        sample_name='vSamp3',
        type_abundance_dict={'A1-A1aw' : 0.4, 'C3-C3a-C3ab': 0.6},
        sequencing_depth=1000)
    # One A1-A1aw, Two C15-C15x, C3-C3a-C3ab
    vSamp4 = virtual_sample(
        sample_name='vSamp4',
        type_abundance_dict={'C3-C3a-C3ab': 1},
        sequencing_depth=1000)
    # One A1-A1aw, Two C15-C15x, C3-C3a-C3ab
    vSamp5 = virtual_sample(
        sample_name='vSamp5',
        type_abundance_dict={'C3-C3c-C3cc': 1},
        sequencing_depth=1000)
    # One A1-A1aw, One C15-C15e
    vSamp6 = virtual_sample(
        sample_name='vSamp6',
        type_abundance_dict={'A1-A1aw' : 0.2, 'C3-C3c-C3cc': 0.8},
        sequencing_depth=1000)
    # One A1-A1aw, One C15-C15e
    vSamp7 = virtual_sample(
        sample_name='vSamp7',
        type_abundance_dict={'C3-C3c-C3cc': 1},
        sequencing_depth=1000)
    # One A1-A1aw, One C3-C3a-C3ab
    vSamp8 = virtual_sample(
        sample_name='vSamp8',
        type_abundance_dict={'C3-C3a-C3ab': 1},
        sequencing_depth=1000)
    # One A1-A1aw, One C3-C3a-C3ab
    vSamp9 = virtual_sample(
        sample_name='vSamp9',
        type_abundance_dict={'C3-C3a-C3ab': 0.8, 'D1-D4-D6' : 0.2},
        sequencing_depth=1000)
    # One A1-A1av-A1ax, One C3-C3a-C3ab
    vSamp10 = virtual_sample(
        sample_name='vSamp10',
        type_abundance_dict={'A1-A1aw' : 0.4, 'C3-C3a-C3ab': 0.6},
        sequencing_depth=1000)
    # One A1-A1av-A1ax, One C3-C3a-C3ab
    vSamp11 = virtual_sample(
        sample_name='vSamp11',
        type_abundance_dict={'A1-A1av-A1ax' : 1},
        sequencing_depth=1000)
    # One A1-A1av-A1ax, One C3-C3a-C3ab
    vSamp12 = virtual_sample(
        sample_name='vSamp12',
        type_abundance_dict={'A1-A1av-A1ax' : 1},
        sequencing_depth=1000)
    # One A1-A1aw, One C15-C115c
    vSamp13 = virtual_sample(
        sample_name='vSamp13',
        type_abundance_dict={'C15-C15x' : 0.4, 'C3-C3a-C3ab': 0.6},
        sequencing_depth=1000)
    # One A1-A1av-A1ax, One C15-C15x
    vSamp14 = virtual_sample(
        sample_name='vSamp14',
        type_abundance_dict={'C15-C15x': 0.9, 'C3-C3a-C3ab': 0.1},
        sequencing_depth=1000)
    # One A1-A1aw, Two C1-C2-C4, C3-C3c-C3cc
    vSamp15 = virtual_sample(
        sample_name='vSamp15',
        type_abundance_dict={'C15-C15x' : 1},
        sequencing_depth=1000)
    # One A1-A1av-A1ax, two C1-C2-C4, C3-C21,
    vSamp16 = virtual_sample(
        sample_name='vSamp16',
        type_abundance_dict={'C1-C1ab-C1e': 1},
        sequencing_depth=1000)
    # One A1-A1av-A1ax, One C3-C3a-C3ab
    vSamp17 = virtual_sample(
        sample_name='vSamp17',
        type_abundance_dict={'C1-C1ab-C1e': 1},
        sequencing_depth=1000)
    # One A1-A1av-A1ax, One C1-C1cc
    vSamp18 = virtual_sample(
        sample_name='vSamp18',
        type_abundance_dict={'A1-A1aw': 0.05, 'C1-C1ab-C1e': 0.9, 'D1-D4-D6': 0.05},
        sequencing_depth=1000)

    vSamp19 = virtual_sample(
        sample_name='vSamp19',
        type_abundance_dict={'C3-C3a-C3ab': 1},
        sequencing_depth=1000)

    vSamp20 = virtual_sample(
        sample_name='vSamp20',
        type_abundance_dict={'A1-A1aw': 0.15, 'C3-C3a-C3ab': 0.85},
        sequencing_depth=1000)

    vSamp21 = virtual_sample(
        sample_name='vSamp21',
        type_abundance_dict={'A1-A1aw': 0.2, 'C3-C3a-C3ab': 0.80},
        sequencing_depth=1000)

    vSamp22 = virtual_sample(
        sample_name='vSamp22',
        type_abundance_dict={'C15-C15x' : 0.2, 'C3-C3a-C3ab': 0.80},
        sequencing_depth=1000)

    vSamp23 = virtual_sample(
        sample_name='vSamp23',
        type_abundance_dict={'C3-C3a-C3ab': 1},
        sequencing_depth=1000)

    vSamp24 = virtual_sample(
        sample_name='vSamp24',
        type_abundance_dict={'C3-C3a-C3ab': 1},
        sequencing_depth=1000)

    vSamp25 = virtual_sample(
        sample_name='vSamp25',
        type_abundance_dict={'A1-A1aw': 0.3, 'C1-C1ab-C1e': 0.7},
        sequencing_depth=1000)

    vSamp26 = virtual_sample(
        sample_name='vSamp26',
        type_abundance_dict={'C3-C3a-C3ab': 1},
        sequencing_depth=1000)

    vSamp27 = virtual_sample(
        sample_name='vSamp27',
        type_abundance_dict={'D1-D4-D6': 1},
        sequencing_depth=1000)

    vSamp28 = virtual_sample(
        sample_name='vSamp28',
        type_abundance_dict={'A1-A1aw': 0.2, 'C3-C3a-C3ab': 0.8},
        sequencing_depth=1000)

    vSamp29 = virtual_sample(
        sample_name='vSamp29',
        type_abundance_dict={'A1-A1av-A1ax': 1},
        sequencing_depth=1000)

    vSamp30 = virtual_sample(
        sample_name='vSamp30',
        type_abundance_dict={'C3-C3a-C3ab': 1},
        sequencing_depth=1000)

    sample_list.extend(
        [vSamp1, vSamp2, vSamp3, vSamp4, vSamp5, vSamp6, vSamp7, vSamp8, vSamp9, vSamp10,
         vSamp11, vSamp12, vSamp13, vSamp14, vSamp15, vSamp16, vSamp17, vSamp18, vSamp19,
         vSamp20, vSamp21, vSamp22, vSamp23, vSamp24, vSamp25, vSamp26, vSamp27, vSamp28,
         vSamp29, vSamp30])

    return sample_list


class virtual_sample():
    ''' The virtual sample will be generated simply by inputing a dict which will be the name of a type
    and the percentage that we want this sample to contain of that type.'''


    def __init__(self, sample_name, type_abundance_dict, sequencing_depth):
        self.name = sample_name
        self.sequencing_depth = sequencing_depth
        self.type_abundance_dict = type_abundance_dict
        self.list_of_ITS2_type_objects = [virtual_ITS2_type_profile(type_name) for type_name, type_abund in self.type_abundance_dict.items()]
        # Generate the absolute sequence abundances both raw and post-qc and MED
        # These get saved to the ITS2_type_profile instances
        self.generate_and_assign_raw_seqs_to_samples_types()
        self.generate_and_assign_post_qc_seqs_to_samples_types()


    def __str__(self):
        return self.name    

    # returns the counts of sequences in a dict of sequence to counts
    # This will work on the assumption that there is a read depth of 1000 seqs
    def generate_and_assign_raw_seqs_to_samples_types(self):
        ''' For each of the types, simply multiply the relative abundances by the sequencing_depth'''
        for i in range(len(self.list_of_ITS2_type_objects)):
            raw_seq_absolute_abundances = []
            for seq_tup in self.list_of_ITS2_type_objects[i].seq_abund_profile:
                seq_name, rel_abund = seq_tup[0], seq_tup[1]
                absolute_abundance_of_seq = rel_abund * self.sequencing_depth * self.type_abundance_dict[self.list_of_ITS2_type_objects[i].name]
                raw_seq_absolute_abundances.append((seq_name, absolute_abundance_of_seq))
            self.list_of_ITS2_type_objects[i].raw_seq_abs_abundances = {seq_name : seq_abund for seq_name, seq_abund in raw_seq_absolute_abundances}

    def generate_and_assign_post_qc_seqs_to_samples_types(self):
        ''' The purpose of this is to make the DIV raw abundances grow slightly to simulate the sequences going
        through the MED and QC process. In reality we would also lose sequences but we are not going to take this
        into account here. we are only bothered about giving the impression visually that MED and qc has occured.
        As such we will increase each of the DIV raw abundaces by a fixed amount and then regenerate the
        non-DIV abundances to make up the remainer of the sequences. To make this realistic, when we recalculate the
        abundances for the non DIV sequences we should make sure that we use much larger percentages so that we end
        up with fewer, better represented non_DIV seqs that will represent the MED nodes.
        We basically need to 1/4 the proportion of sequences that are non-DIV. We should work out the proportion
        that we should be gaining and add this onto the abundances of the DIVs'''

        # For each of the types in the samples

        for i in range(len(self.list_of_ITS2_type_objects)):
            # Firstly go through and identify how much of the types rel abundances are taken up by non-DIVs
            current_total_DIV_proportion = 0
            for j in range(len(self.list_of_ITS2_type_objects[i].seq_abund_profile)): # for each seq, abund pair
                seq_name, seq_abund = self.list_of_ITS2_type_objects[i].seq_abund_profile[j]
                if 'Unk' not in seq_name: # then this is a DIV seq
                    current_total_DIV_proportion += seq_abund
            current_total_non_DIV_proportion = 1 - current_total_DIV_proportion

            # Divide this by 4. Multiply by 3. This is the relproportions total that we can add to the DIV abundances
            # We should add this abundance to the current abundances in proportion to the DIVs current abundances
            proportion_of_non_DIV_seqs_to_collapse_into_DIVs = .75 * current_total_non_DIV_proportion

            # calculate the abundance currently taken up by the DIVs.
            # Then for each DIV in the type, add the total relproportion to be redistributed * the relabundance of the
            # DIV / total abundances taken up by the div). THis will then be the new absolute abundance of the sequence

            MED_adjusted_DIV_abundances = []
            DIV_adjusted_rel_abund_total = 0
            for j in range(len(self.list_of_ITS2_type_objects[i].seq_abund_profile)): # for each seq, abund pair
                seq_name, seq_abund = self.list_of_ITS2_type_objects[i].seq_abund_profile[j]
                # Only do this for the DIVs
                if 'Unk' not in seq_name:
                    adjusted_rel_abund = seq_abund + (proportion_of_non_DIV_seqs_to_collapse_into_DIVs * (seq_abund/current_total_DIV_proportion))
                    DIV_adjusted_rel_abund_total += adjusted_rel_abund
                    MED_adjusted_DIV_abundances.append((seq_name, adjusted_rel_abund))

            # Then we add up these new abundances and carry them into a new generate MED non-DIV samples
            # use the same distributions as used previously to genreate non-DIV samples but this time increase them slightly
            # once these have been calcualted add them to the types post_qc_abs_aundances parameter.
            total = DIV_adjusted_rel_abund_total
            seqNameCounter = 0
            while total < 1:
                increment = numpy.random.choice(numpy.arange(0.005, 0.025, 0.005), p=[0.25, 0.25, 0.25, 0.25])
                # Make sure that the gap isn't smaller than the smallest increment
                if (1 - total) < 0.005:
                    MED_adjusted_DIV_abundances.append(('{}_Unk_{}'.format(self.list_of_ITS2_type_objects[i].name, seqNameCounter), 1 - total))
                    seqNameCounter += 1
                    total += 1 - total
                
                elif total + increment <= 1:  # Then the increment is a good size, add
                    MED_adjusted_DIV_abundances.append(('{}_Unk_{}'.format(self.list_of_ITS2_type_objects[i].name, seqNameCounter), increment))
                    seqNameCounter += 1
                    total += increment
            
            # Here we have a full tup list of seq to adjusted_for_MED_rel_abundance
            # We can now convert these to absolute abundances
            raw_seq_absolute_abundances = []
            for seq_name, rel_abund in MED_adjusted_DIV_abundances:
                absolute_abundance_of_seq = rel_abund * self.sequencing_depth * self.type_abundance_dict[self.list_of_ITS2_type_objects[i].name]
                raw_seq_absolute_abundances.append((seq_name, absolute_abundance_of_seq))
            self.list_of_ITS2_type_objects[i].post_qc_abs_aundances = {seq_name: seq_abund for seq_name, seq_abund in raw_seq_absolute_abundances}

        
class virtual_ITS2_type_profile():
    # This will be a virutal type profile
    # when we are generating samples we will define what percentage of which type we wish to have in the sample
    # The actual sequence abundances will be generated for each of the types in this class. The sequences
    # will be generated on a mean and standard deviation basis.
    def __init__(self, typeName):
        self.name = typeName
        self.seqs, self.means = self.getSeqList()
        self.seq_abund_profile = self.generate_relative_abundances_of_type_sequences()
        self.clade = self.seqs[0][0]
        # lists of tuples where each tuple is the name of the sequence with its absolute abundance
        self.raw_seq_abs_abundances  = {}
        self.post_qc_abs_aundances = {}

    def __str__(self):
        return self.name

    def getSeqList(self):
        # CLADE A types
        # For each type have the list of DIV names and the average DIV abundances
        if self.name =='A1-A1aw':
            return (['A1', 'A1aw'], [0.45, 0.25])
        elif self.name == 'A1-A1av-A1ax':
            return (['A1', 'A1av', 'A1ax'], [0.4, 0.25, 0.1])
        # Clade C types
        elif self.name == 'C15-C15x':
            return (['C15', 'C15x'], [0.5, 0.2])
        elif self.name == 'C15-C15e':
            return (['C15', 'C15e'], [0.4, 0.25])
        elif self.name == 'C15-C115c':
            return (['C15', 'C115c'], [0.45, 0.2])
        elif self.name == 'C1-C1c-C1d':
            return (['C1', 'C1c', 'C1d'], [0.3, 0.25, 0.2])
        elif self.name == 'C1-C1cc':
            return (['C1', 'C1cc', 'C1cd', 'C1ce'], [0.6, 0.07, 0.05, 0.05])
        elif self.name == 'C1-C1ab-C1e':
            return (['C1', 'C1ab', 'C1e'], [0.4, 0.2, 0.18])
        elif self.name == 'C3-C3a-C3ab':
            return (['C3', 'C3a', 'C3ab'], [0.3, 0.28, 0.18])
        elif self.name == 'C3-C21':
            return (['C3', 'C21'], [0.35, 0.28])
        elif self.name == 'C3-C3c-C3cc':
            return (['C3', 'C3c', 'C3cc'], [0.3, 0.25, 0.15])
        # Clade D types
        elif self.name == 'D1-D4-D6':
            return (['D1', 'D4', 'D6'], [0.25, 0.25, 0.15])

    def generate_relative_abundances_of_type_sequences(self):
        '''Generate a tuple list that is sequence to relative proportion
        where sequenes are defined by those present in self. seqs which is a list
        the relative proportions will be generated according to a normal distribution around
        the means given in self.means with a fixed standard deviation of the mean (0.1)'''

        # make sure that the named sequences don't make up more than 90% of the sequences

        list_of_random_proportions = []
        while 1:
            list_of_random_proportions = []
            for seq, abund in zip(self.seqs, self.means):
                # Generate the named seqs first and then generate the unknown seqs
                # Consider that the named seqs must be present, so if a negative value is generated, then run again
                while 1:
                    temp = random.gauss(abund, 0.01)
                    if temp > 0:
                        list_of_random_proportions.append(temp)
                        break
            if sum(list_of_random_proportions) < 0.9:
                break

        # I want to have three different gradients for creating the unknown sequences
        # I will have one third of the unknown proportion generated with bigger, one same, one less
        # Now generate the unknown seqs
        total = sum(list_of_random_proportions)
        seqNameCounter = 0
        
        total_to_fill_with_unkn = 1 - total
        unknown_total = 0
        while total < 1:
            increment = int()
            if unknown_total < total_to_fill_with_unkn * 0.333: # Then do the bigger fills
                increment = numpy.random.choice(numpy.arange(0.01, 0.1, 0.01),
                                                p=[0.025, 0.025, 0.025, 0.025, 0.1, 0.1, 0.1, 0.2, 0.4])
            elif unknown_total > total_to_fill_with_unkn * 0.666: # Then do the smaller fills
                increment = numpy.random.choice(numpy.arange(0.0005, 0.005, 0.0005),
                                                p=[0.2, 0.2, 0.1, 0.1, 0.1, 0.125, 0.125, 0.025, 0.025])
            else: # do the middle fills
                increment = numpy.random.choice(numpy.arange(0.001, 0.01, 0.001),
                                                p=[0.4, 0.2, 0.1, 0.1, 0.1, 0.025, 0.025, 0.025, 0.025])
            # Make sure that the gap isn't smaller than the smallest increment
            if (1 - total) < 0.001:
                list_of_random_proportions.append((1 - total))
                self.seqs.append('{}_Unk_{}'.format(self.name, seqNameCounter))
                seqNameCounter += 1
                unknown_total += 1 - total
            elif total + increment <= 1: # Then the increment is a good size, add
                list_of_random_proportions.append(increment)
                self.seqs.append('{}_Unk_{}'.format(self.name, seqNameCounter))
                seqNameCounter += 1
                unknown_total += increment
            total = sum(list_of_random_proportions)


        relAbundTups = []
        for i in range(len(self.seqs)):
            relAbundTups.append((self.seqs[i], list_of_random_proportions[i]))

        return relAbundTups




fig3_fixed_types()