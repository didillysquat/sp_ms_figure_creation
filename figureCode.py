import matplotlib.pyplot as plt
import os
from collections import defaultdict
import subprocess
import numpy
import pickle
from plumbum import local
import random
from Bio import SeqIO
import pandas as pd
from multiprocessing import Queue, Process, Manager
from Bio.Align import AlignInfo
import Bio.AlignIO

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
#### These two are the two used for the SP figures in the end
def figure2RawSeqs_remake_with_fastq():
    with open('/home/humebc/projects/SymPortalMS/making_example_data_set/SP_DEMO_DATA/unzipped/stability.trim.contigs.pcr.names') as f:
        name_file = f.read().splitlines()


    # The first thing we need to do is to work out which of the representative sequences are clade C sequences
    # we will use blast to do this.
    # We will need the fasta that is associated with the .names file to use as input.

    # write out the .ncbirc file pointing towards the ITS2clade.fas
    ncbirc = ['[BLAST]', 'BLASTDB={}'.format('/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload')]
    ncbirc_path = '/home/humebc/projects/SymPortalMS/.ncbirc'
    with open(ncbirc_path, 'w') as f:
        for item in ncbirc:
            f.write("{}\n".format(item))


    # setup the blast
    outputpath = '/home/humebc/projects/SymPortalMS/sym_clade.out'
    outputFmt = '6 qseqid sseqid'
    inputPath = '/home/humebc/projects/SymPortalMS/making_example_data_set' \
                '/SP_DEMO_DATA/unzipped/stability.trim.contigs.unique.rc.pcr.fasta'
    completedProcess = subprocess.run(
        ['blastn', '-out', outputpath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'ITS2Clade.fas',
         '-max_target_seqs', '1', '-num_threads', '10'])


    # Read in blast output
    with open(outputpath, 'r') as f:
        blast_output = f.read().splitlines()

    # get list of representative sequences that are clade C
    clade_C_rep_seqs = [line.split('\t')[0] for line in blast_output if line.split('\t')[1][5] == 'C']





    # create a dict that asociates each sequence to its representative sequecne in the .names file
    # at the same time use a default dict to get an ordered list of the representative sequences
    # in the order of most sequences representated first.
    # This will be the order in which we plot the sequences

    # Before generating from scratch, try to see if there is a pickled seq_to_rep_dict and sorted_rep_seq_list already
    # If generating from scratch, pickle once created to save time in future
    seq_to_rep_dict_path = '/home/humebc/projects/SymPortalMS/seq_to_rep_dict_path'
    sorted_rep_seq_list_path = '/home/humebc/projects/SymPortalMS/sorted_rep_seq_list'
    try:
        seq_to_rep_dict = pickle.load(open(seq_to_rep_dict_path, 'rb'))
        sorted_rep_seq_list = pickle.load(open(sorted_rep_seq_list_path, 'rb'))

    except:
        seq_to_rep_dict = {}
        rep_seq_default_dict = defaultdict(int)
        for line in name_file:
            if line.split('\t')[0] in clade_C_rep_seqs:
                rep_seq_default_dict[line.split('\t')[0]] += len(line.split('\t')[1].split(','))
                print('Populating {} sequences for the representative sequence: {}'.format(len(line.split('\t')[1].split(',')), line.split('\t')[0]))
                temp_dict = {seq : line.split('\t')[0] for seq in line.split('\t')[1].split(',')}
                seq_to_rep_dict.update(temp_dict)

        sorted_rep_seq_list = sorted(rep_seq_default_dict, key=rep_seq_default_dict.__getitem__, reverse=True)


        pickle.dump(seq_to_rep_dict, open(seq_to_rep_dict_path, 'wb+'))
        pickle.dump(sorted_rep_seq_list, open(sorted_rep_seq_list_path, 'wb+'))


    # get a list of the samples in order that will look good graphically
    # we will order the samples according first to their species
    # secondly by the ITS2 type profile they contained
    # and thirdly by their sample name

    # read in the sample list genreated previousy
    with open('/home/humebc/projects/SymPortalMS/sample_order.txt') as f:
        sample_order_list = f.read().splitlines()

    # create the dict that will hold all of the info to be plotted
    # Initially this will be a 2D list. In each list we be a relative abund of each representative
    # sequence in the order of the sorted_rep_seq_list. These lists will be added in the order of the samples in
    # sample_order_list
    # this will then be transformed so that we have a list per representative sequence in the order of
    # sorted_rep_seq_list. Each list will then have an element which is the relabund of the seq in each sample
    # in the order of sample_order_list. This is the datastructure required for the matplot lib stacked bar
    # We will first try this non-multiprocessed but it may be very slow.

    # do the pickle check again
    plotting_info_holder_list_path = '/home/humebc/projects/SymPortalMS/plotting_info_holder_list'
    try:
        plotting_info_holder_list = pickle.load(open(plotting_info_holder_list_path, 'rb'))
    except:
        plotting_info_holder_list_to_transpose = []
        for sample in sample_order_list:
            print('Populating plotting_info_holder_list for sample: {}'.format(sample))
            # create a counter dict to keep track of abundance of sequences
            counter_dict = defaultdict(int)
            counter_dict.update((rep_seq, 0) for rep_seq in sorted_rep_seq_list)
            sample_fastq_path = '/home/humebc/projects/SymPortalMS/making_example_data_set/SP_DEMO_DATA/unzipped/{}.1_its2_only.fastq'.format(sample)
            # use the SeqIO to easily parse through the fastq format rather than reinventing wheel
            sample_fastq_seqs = list(SeqIO.parse(sample_fastq_path, "fastq"))

            # for each sequence in the fastq add 1 to its represenatives value in the counter_dict
            # get the representative sequence from the seq_to_rep_dict created earlier
            for record in sample_fastq_seqs:
                # only for clade C seqs
                if record.description[:-2] in seq_to_rep_dict:
                    counter_dict[seq_to_rep_dict[record.description[:-2]]] += 1


            # now create a list of the normalised abundances in order of sorted_rep_seq_list
            tot_seqs = sum(counter_dict.values())
            norm_abunds_for_sample = [counter_dict[rep_seq]/tot_seqs for rep_seq in sorted_rep_seq_list]

            # now add this list to the the plotting_info_holder_list
            plotting_info_holder_list_to_transpose.append(norm_abunds_for_sample)

        # here we have the plotting_info_holder_list populated
        # we should now transpose it to be compatible with the matplotlib
        plotting_info_holder_list = list(map(list, zip(*plotting_info_holder_list_to_transpose)))
        pickle.dump(plotting_info_holder_list, open(plotting_info_holder_list_path, 'wb+'))


    ###### PLOTTING
    # the greys for plotting
    greyPalette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']


    # Set up plot
    plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1, 1, 1)
    ax.set_ylim(bottom=0, top=1)

    ind = range(len(sample_order_list))
    width = 1

    # Counters to work out which grey we are on
    greyCounter = 0

    # The bottom for each sequence plot will be updated each time we plot
    bottom = [0 for smpl in sample_order_list]


    # Only plot the first 250 seqs as this is the majority of the volume
    # After that plot the rest of the volume as fixed height boxes
    # NB the greyFileHeight should be the same as the max value (or not smaller) of the max random generated
    greyFileHeight = 0.01
    for i in range(len(sorted_rep_seq_list)):
        if i < 250:

            print('Plotting seq {} out of {}'.format(i, len(sorted_rep_seq_list)))

            plotColour = greyPalette[greyCounter % 4]
            greyCounter += 1

            ax.bar(ind, plotting_info_holder_list[i], width, bottom=bottom, color=plotColour)
            bottom = [L + M for L, M in zip(bottom, plotting_info_holder_list[i])]
        else:


            print('Plotting seq {} out of {}'.format(i, len(sorted_rep_seq_list)))
            plotValuesList = []
            for bValue in bottom:
                if bValue < 100:
                    # If this is not the last bar then plot a standard bar
                    if bValue < 100 - greyFileHeight:
                        plotValuesList.append(numpy.random.choice(numpy.arange(0.0001, 0.001, 0.0001),
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

            plotColour = greyPalette[greyCounter % 4]
            greyCounter += 1

            ax.bar(ind, plotValuesList, width, bottom=bottom, color=plotColour)
            bottom = [L + M for L, M in zip(bottom, plotValuesList)]


    # Here we should have it plotted and we'll just need to put the finishing touches to it.
    plt.savefig("/home/humebc/projects/SymPortalMS/fig_3_raw_100718.svg")
    plt.show()
    plt.savefig("/home/humebc/projects/SymPortalMS/fig_3_raw_100718.svg")

    return

def figure2_afterMED_with_SP_output():
    ''' This will produce the coloured stacked bar chat that represents the sequences per sample after they
    have been through the SP QC including the MED and the blasting.'''

    # read in the list that contains the DIVs that were used to create the ITS2 types.
    # These will be coloured in the plot.
    with open('/home/humebc/projects/SymPortalMS/significantSeqs', 'r') as f:
        significant_seqs_fasta_list = f.read().splitlines()

    significant_seqs_list = []
    for line in significant_seqs_fasta_list:
        if line.startswith('>'):
            significant_seqs_list.append(line[1:])

    # we will use a SymPortal output to make the graph
    # this will save us a lot of time.
    # we will import this as a pandas dataframe
    plotting_data_df = pd.read_csv('/home/humebc/projects/SymPortalMS/SP_output_for_figures/30_141.DIVs.absolute.txt',
                                    delimiter='\t', header=0, index_col=0)

    # now change the order of the index of the df according to the predefined order (same as raw plot)
    # read in the sample list genreated previousy
    with open('/home/humebc/projects/SymPortalMS/sample_order.txt') as f:
        sample_order_list = f.read().splitlines()

    plotting_data_df = plotting_data_df.reindex(sample_order_list)


    # now drop the columns that contain the noName summaries and the QC info
    # the first of the columns we want to drop
    first_column_index = plotting_data_df.columns.get_loc("noName Clade A")
    # last of the columns we want to drop
    last_columnn_index = plotting_data_df.columns.get_loc("post_med_unique")
    plotting_data_df.drop(plotting_data_df.columns[first_column_index:last_columnn_index + 1], axis=1, inplace=True)

    # we also need to drop the colums that are not clade C seqs
    index_to_drop = []
    for seq_name in list(plotting_data_df):
        if 'C' not in seq_name:
            index_to_drop.append(plotting_data_df.columns.get_loc(seq_name))

    plotting_data_df.drop(plotting_data_df.columns[index_to_drop], axis=1, inplace=True)

    # finally we need to normalise each of the seq abundances to the total seq abundances found in each sample
    # i.e. divide each cell by the sum of its row
    plotting_data_df = plotting_data_df.div(plotting_data_df.sum(axis=1), axis=0)


    # we can directly use the plotting_data_df for plotting
    # we will go column by column when plotting
    # therefore at this stage we have all of the info we need.
    # we just need to work out which of the sub bars should be coloured and what colour they should be


    # Colour palettes
    colourPalette = ['#E8AB52', '#5FA96E', '#9680D3', '#D8C542', '#D47070', '#59A6D4', '#D76FBC', '#F6724C',
                     '#5FC2B6', '#80D1E4', '#6BD78A', '#B7D456']
    greyPalette = ['#D0CFD4', '#89888D', '#4A4A4C', '#8A8C82', '#D4D5D0', '#53544F']

    # create a seq to colour dictionary that we can use when plotting
    seq_to_colour_dict = {}
    colour_counter = 0
    grey_counter = 0

    # for each sequence
    for sequence_to_plot in list(plotting_data_df):
        if sequence_to_plot in significant_seqs_list:
            seq_to_colour_dict[sequence_to_plot] = colourPalette[colour_counter]
            colour_counter += 1
        else:
            seq_to_colour_dict[sequence_to_plot] = greyPalette[grey_counter%6]
            grey_counter += 1

    # we are now ready to plot
    plt.figure(num=1, figsize=(8, 10), dpi=80)
    ax = plt.subplot(1, 1, 1)
    ax.set_ylim(bottom=0, top=1)

    ind = range(len(sample_order_list))
    width = 1

    # The bottom for each sequence plot will be updated each time we plot
    bottom = [0 for smpl in sample_order_list]


    ordered_seq_list = list(plotting_data_df)
    for i in range(len(ordered_seq_list)):
        seq_name_in_Q = ordered_seq_list[i]
        print('Plotting seq {} out of {}'.format(i, len(ordered_seq_list)))
        ax.bar(ind, plotting_data_df[seq_name_in_Q], width, bottom=bottom, color=seq_to_colour_dict[seq_name_in_Q])
        bottom = [L + M for L, M in zip(bottom, plotting_data_df[seq_name_in_Q])]


    # Here we should have it plotted and we'll just need to put the finishing touches to it.
    plt.savefig("/home/humebc/projects/SymPortalMS/fig_3_colour_090718.svg")
    plt.show()
    return
##########################################

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
    ''' The aim of this function will be to generate figure two. It will have a type class and a sample class.
    Each sample will be generated with the same amount of the same types each time. Samples will be a mix of
    clade A, C and D types but mainly, C.
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
    plt.show(bob)
    # plt.savefig('fig3_raw.svg', format="svg")
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
        plt.show(bob)
        # plt.savefig('fig3_raw_{}.svg'.format(clade), format="svg")

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

        plt.show(bob)
        # plt.savefig('fig3_MED_post_qc_{}.svg'.format(clade), format="svg")

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
        plt.show(bob)
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


    # plt.savefig('fig3_ITS2_type_profile_final.svg', format="svg")
    plt.show(bob)
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

def psba_phylogeny():
    ''' The aim of this method is to make a phylgoentic tree of the samples from Ed's study using the psba seqs.
    We will work with the psba fastas that he has in the dryad repo from his project
    https://datadryad.org/resource/doi:10.5061/dryad.h6s54.
    We should go file by file and read in the sequences using Bio.AlignIO.read()
    Then:
    from Bio.Align import AlignInfo
    summary_align = AlignInfo.SummaryInfo(alignment)
    we can then calculate a rough consensus sequence from this summary object
    consensus = summary_align.dumb_consensus()
    '''


    # working directory with the fasta files of the psba sequences in them
    wkd = '/home/humebc/projects/SymPortalMS/making_example_data_set/SP_DEMO_DATA/psbA_processing'

    # get a list of the fasta files
    list_of_filenames = []
    for (dirpath, dirnames, filenames) in os.walk(wkd):
        list_of_filenames.extend([fasta_file for fasta_file in filenames if '.fasta' in fasta_file])
        break

    # setup MP
    # Create the queues that will hold the cc information
    taskQueue = Queue()

    # output queue to put the dictionaries that will make the fasta and name files
    output_manager = Manager()
    output_manager_list = output_manager.list()

    # populate the input queue
    for file_name in list_of_filenames:
        taskQueue.put(file_name)

    num_proc = 20

    # place the stop cues
    for n in range(num_proc):
        taskQueue.put('STOP')

    allProcesses = []

    for n in range(num_proc):
        p = Process(target=consensus_worker,
                    args=(taskQueue, output_manager_list, wkd))
        allProcesses.append(p)
        p.start()

    # wait for the processes to finish before processing the output list
    for p in allProcesses:
        p.join()

    # process the output files from the MP which will be tuples of the sequence name and sequence
    consensus_sequence_fasta = []
    for tup_fasta_info in list(output_manager_list):
        consensus_sequence_fasta.extend([tup_fasta_info[0], tup_fasta_info[1]])


    # Here we should ahve the consensus_seqeunce_fasta popuated and we can now write it out
    with open('{}/{}'.format(wkd, 'consensus_psbs.fasta'), 'w') as f:
        for line in consensus_sequence_fasta:
            f.write('{}\n'.format(line))


def consensus_worker(input, output_list, wkd):
    for fasta_filename in iter(input.get, 'STOP'):
        print('Getting consensus for {}'.format(fasta_filename))
        alignment_object = Bio.AlignIO.read('{}/{}'.format(wkd, fasta_filename), 'fasta')
        summary_object = AlignInfo.SummaryInfo(alignment_object)
        consensus = summary_object.dumb_consensus()
        output_list.append(('>{}'.format(fasta_filename), str(consensus)))

def psba_phylogeny_from_scratch():
    ''' When we tried making the phylogeny from the seuences ed had in dryad it was a nice split between the
    three species. However, the seqs were only 220 bp long and so I think we can do better.
    Here we will work on individual forward and reverse reads.
    We will put them through mothur screening and then blast to get rid of the As.
    Then we can have a look at the name files and see if there are sequences that really stand out,
    Then perhaps we work with those. FOr making the blast dictionaries we should look at the most abundant seqs in
    one sample from each species that is known to have low A clade. So, C10, A02 and Y12.'''

    # working directory with the fasta files of the psba sequences in them
    wkd = '/home/humebc/projects/SymPortalMS/making_example_data_set/SP_DEMO_DATA/psbA_processing_raw_fastq'

    # get a list of the fasta files
    list_of_filenames = []
    for (dirpath, dirnames, filenames) in os.walk(wkd):
        list_of_filenames.extend([fastq_file for fastq_file in filenames if '.fastq' in fastq_file])
        break

    if not list_of_filenames:
        # this may be caused due to the fact that we are runnin this for the second time and the fastq files have
        # already been moved into their respective directories
        # in which case get names from the dir names
        sample_names_list = []
        for (dirpath, dirnames, filenames) in os.walk(wkd):
            sample_names_list.extend(dirnames)
            break
    else:
        # get a list of the sample names
        sample_names_set = set()
        for file_name in list_of_filenames:
            sample_names_set.add(file_name.split('.')[0])
        sample_names_list = list(sample_names_set)

        # make a directory for each of the samples
        for smp in sample_names_list:
            os.makedirs('{}/{}'.format(wkd, smp), exist_ok=True)


    # set up the MP queues
    input_queue = Queue()

    # now move the files into their respective directories
    # and add two oligo files into each directory
    # also add the mBatchFile_1 and mBatchFile_2
    oligo_1 = ['forward AAGCAYCCAATRTAGAGACGATTTGYTGTGG']
    oligo_2 = ['forward GGWATGGAAGTVATGCATGAAAGAAAYGC']
    for smp in sample_names_list:
        smp_wkd = '{}/{}'.format(wkd, smp)
        # move both the fastq files for the sample into the directory
        if os.path.exists('{}/{}.1.fastq'.format(wkd, smp)):
            subprocess.run([r'mv', r'{}/{}.1.fastq'.format(wkd, smp),
                            r'{}/{}/{}.1.fastq'.format(wkd, smp, smp)],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if os.path.exists('{}/{}.2.fastq'.format(wkd, smp)):
            subprocess.run([r'mv', r'{}/{}.2.fastq'.format(wkd, smp),
                            r'{}/{}/{}.2.fastq'.format(wkd, smp, smp)],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # write out the oligo files that will be used for the pcr step
        # write out oligo_1
        # this will be the oligo file for the .1 files (even though its the reverse primer)
        with open('{}/{}/primers.oligos_1'.format(wkd, smp), 'w') as f:
            for line in oligo_1:
                f.write('{}\n'.format(line))
        # write out oligo_2
        # this will be the oligo file for the .2 files (even though its the reverse primer)
        with open('{}/{}/primers.oligos_2'.format(wkd, smp), 'w') as f:
            for line in oligo_2:
                f.write('{}\n'.format(line))

        # put in the mothur batch files that will be run for each of the fastqs
        mBatch_contents_1 = [
            'set.dir(input={})'.format(smp_wkd),
            'set.dir(output={})'.format(smp_wkd),
            'fastq.info(fastq={}/{}.1.fastq)'.format(smp_wkd, smp),
            'trim.seqs(fasta={0}/{1}.1.fasta, qwindowaverage=25, qwindowsize=4, '
            'qfile={0}/{1}.1.qual)'.format(smp_wkd, smp),
            'pcr.seqs(fasta={0}/{1}.1.trim.fasta, oligos={0}/primers.oligos_1, pdiffs=3)'.format(smp_wkd, smp),
            'unique.seqs(fasta={}/{}.1.trim.pcr.fasta)'.format(smp_wkd, smp)
        ]

        mBatch_contents_2 = [
            'set.dir(input={})'.format(smp_wkd),
            'set.dir(output={})'.format(smp_wkd),
            'fastq.info(fastq={}/{}.2.fastq)'.format(smp_wkd, smp),
            'trim.seqs(fasta={0}/{1}.2.fasta, qwindowaverage=25, qwindowsize=4, '
            'qfile={0}/{1}.2.qual)'.format(smp_wkd, smp),
            'pcr.seqs(fasta={0}/{1}.2.trim.fasta, oligos={0}/primers.oligos_2, pdiffs=3)'.format(smp_wkd, smp),
            'unique.seqs(fasta={}/{}.2.trim.pcr.fasta)'.format(smp_wkd, smp)
        ]

        with open('{}/mothur_batch_1'.format(smp_wkd), 'w') as f:
            for line in mBatch_contents_1:
                f.write('{}\n'.format(line))

        with open('{}/mothur_batch_2'.format(smp_wkd), 'w') as f:
            for line in mBatch_contents_2:
                f.write('{}\n'.format(line))

        # At this point we should be ready to run in this directory
        # We will add the sample name to the input list
        # once we have all of the samples added we will start the MPing
        input_queue.put(smp)

    num_proc = 24
    for n in range(num_proc):
        input_queue.put('STOP')

    allProcesses = []

    for n in range(num_proc):
        p = Process(target=mothur_worker,
                    args=(input_queue, wkd))
        allProcesses.append(p)
        p.start()

    # wait for the processes to finish before processing the output list
    for p in allProcesses:
        p.join()

    # at this point we have the mothur files run and so now we have the .names and .fasta files for each sample
    # now we need to get rid of sequecnes that aren't clade C
    # to do this we should run all sequencs against a blast database
    # to create the blast database we need to identify example C psbA sequences
    # best way to do that is to look at the most common sequences in samples that contain very little A.
    # according to the SP ITS2 sequence output these samples are:
    samples_to_check = ['C10', 'A02', 'Y12']

    # for each of the these samples, read in the fasta and .names files
    # find the two most abundant sequeneces (two for the .1 files and two for the .2 file)
    # and add these to the fasta_list that we will turn into
    # a blastn database using makeblastdb
    # then we can go sample by sample and run the sequences in each samples fasta file against this database
    # and only keep the matches.
    # Given the high variability of the psbA region, it should be that clade A sequences simply don't return matches
    # once we have gotten rid of the other clade sequences it will be time to work out how to get a reliable
    # consensus for each of the samples' fwd and rev reads.

    # add the abundant sequences from each of the 'low A' samples
    abundant_seqs_list = []
    for sample_to_check in samples_to_check:
        smp_wkd = '{}/{}'.format(wkd, sample_to_check)

        # for both the .1 and the .2 files
        for file_num in [1,2]:
            # read in fasta
            with open('{}/{}.{}.trim.pcr.unique.fasta'.format(smp_wkd, sample_to_check, file_num)) as f:
                fasta_file = [line.rstrip() for line in f]
            # make fasta dict
            fasta_dict = {fasta_file[i].split('/')[0][1:] : fasta_file[i + 1] for i in range(len(fasta_file)) if fasta_file[i].startswith('>')}

            # read in .names file
            with open('{}/{}.{}.trim.pcr.names'.format(smp_wkd, sample_to_check, file_num)) as f:
                names_file = [line.rstrip() for line in f]
            # make .names dict
            names_dict = {line.split('\t')[0][:-2] : len(line.split('\t')[1].split(',')) for line in names_file}


            # get list of most abundant representative sequences
            sorted_rep_seq_list = sorted(names_dict, key=names_dict.__getitem__, reverse=True)
            # add first and second most abundant
            abundant_seqs_list.extend(['>{}'.format(sorted_rep_seq_list[0]), fasta_dict[sorted_rep_seq_list[0]]])
            abundant_seqs_list.extend(['>{}'.format(sorted_rep_seq_list[1]), fasta_dict[sorted_rep_seq_list[1]]])

    # now write out fasta
    with open('{}/psbA_C.fasta'.format(wkd), 'w') as f:
        for line in abundant_seqs_list:
            f.write('{}\n'.format(line))

    # create blast db
    subprocess.run(
        ['makeblastdb', '-in', '{}/psbA_C.fasta'.format(wkd), '-dbtype', 'nucl', '-title', 'psbA_C'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    apples = 'asdf'
    # now we can cycle back through the samples performing a blast to get a list of samples that are
    # most likely clade C.
    # when in the sample directory, get the most two most abundant of the seqs that are Cs for .1 and .2 file
    # then confirm that each of these is a subset of each other (the difference should just be due to length)
    # then get rid of the longer sequence.
    # add the shorter sequence for each file to a file master alignment (perhaps have an alginment for each file).
    # Once we have the sequences in the master alignement, for each set of the .1 and .2, crop to the shortest seq
    # then we can concatenate and work with this concatenated sequences for the phylogenetic analysis

    # this will hold the fastas for each of the files e.g. .1 and .2 files
    fasta_holder_list = [[],[]]

    for smp in sample_names_list:
        smp_wkd = '{}/{}'.format(wkd, smp)
        for file_num in [1,2]:
            # print('Identifying most abundant C seqs for {}:{}'.format(smp, file_num))
            input_fasta_path = clean_fasta_of_primer_match_info(path_to_dirty_fasta='{}/{}.{}.trim.pcr.unique.fasta'.format(smp_wkd, smp, file_num))
            blastOutputPath = '{}/blast_{}.out'.format(smp_wkd, file_num)
            outputFmt = '6 qseqid sseqid evalue pident qcovs sstrand'

            # create and write out .ncbirc file
            ncbircFile = []
            ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(wkd)])
            with open('{}/.ncbirc'.format(smp_wkd), 'w') as f:
                for line in ncbircFile:
                    f.write('{}\n'.format(line))

            os.chdir(smp_wkd)

            # Run local blast
            subprocess.run([
                'blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', input_fasta_path, '-db', 'psbA_C.fasta',
                 '-max_target_seqs', '1', '-num_threads', '24'])

            # Read in ouput file
            with open(blastOutputPath, 'r') as f:
                blast_output_file = [line.rstrip() for line in f]

            # Get list of the C sequences
            c_seqs = [line.split('\t')[0][:-2] for line in blast_output_file]

            # Now grab the two most abundant sequences for the sample/file that are c_seqs
            # read in fasta
            with open('{}/{}.{}.trim.pcr.unique.fasta'.format(smp_wkd, smp, file_num), 'r') as f:
                fasta_file = [line.rstrip() for line in f]
            # make fasta dict
            fasta_dict = {fasta_file[i].split('/')[0][1:]: fasta_file[i + 1] for i in range(len(fasta_file)) if
                          fasta_file[i].startswith('>')}

            # read in .names file
            with open('{}/{}.{}.trim.pcr.names'.format(smp_wkd, smp, file_num), 'r') as f:
                names_file = [line.rstrip() for line in f]
            # make .names dict
            names_dict = {line.split('\t')[0][:-2]: len(line.split('\t')[1].split(',')) for line in names_file}

            print('Sample {}:{} has {} unique seqs'.format(smp, file_num, len(fasta_dict)))
            # get list of most abundant representative sequences
            sorted_rep_seq_list = sorted(names_dict, key=names_dict.__getitem__, reverse=True)

            # verify that the most abundant and second most abundant seqs are subset of each other
            count = 0
            temp_abundant_seq = []
            for i in range(len(sorted_rep_seq_list)):
                if count < 2:
                    if sorted_rep_seq_list[i] in c_seqs:
                        temp_abundant_seq.append(fasta_dict[sorted_rep_seq_list[i]])
                        count += 1
                else:
                    break
            # Here we have the two most abundant sequences that are C sequences in the temp_abundant_seq list
            # check to see subset and add shortest to the fasta for the sample and file
            if temp_abundant_seq[0] in temp_abundant_seq[1] or temp_abundant_seq[1] in temp_abundant_seq[0]:
                fasta_holder_list[file_num - 1].extend(['>{}'.format(smp), min(temp_abundant_seq, key=len)])
            else:
                print('{} sample did not agree'.format(smp))

    # here we have the fasta_holder_list which has two fastas in it one for each file

    paths_of_cropped_alignments = []
    for file_num in [1,2]:

        # write out the fastas
        with open('{}/psbAncr_fasta_{}.fasta'.format(wkd, file_num), 'w') as f:
            for line in fasta_holder_list[file_num -1]:
                f.write('{}\n'.format(line))

        # align with MAFFT # easier to do this calling process.run rather than the biopython wrapper
        mafft = local["mafft"]
        in_file = '{}/psbAncr_fasta_{}.fasta'.format(wkd, file_num)
        out_file = '{}/psbAncr_aligned_{}.fasta'.format(wkd, file_num)
        # now run mafft including the redirect
        (mafft['--auto', '--thread', '16', in_file] > out_file)()

        # crop to the shortest seq
        cropped_alignment_path = crop_to_shortest_seq(mafft_aligned_fasta_path=out_file)
        paths_of_cropped_alignments.append(cropped_alignment_path)


    # concatenate the two alignments
    # best way to do this is likely is just to read in the two alignments and iterate creating a new alignment
    # read in _1
    with open(paths_of_cropped_alignments[0], 'r') as f:
        fasta_1 = [line.rstrip() for line in f]
    # read in _2
    with open(paths_of_cropped_alignments[0], 'r') as f:
        fasta_2 = [line.rstrip() for line in f]
    # create concatenated alignment
    concat_alignment_fasta = []
    for i in range(len(fasta_1)):
        if fasta_1[i].startswith('>'):
            concat_alignment_fasta.append(fasta_1[i])
        else:
            concat_alignment_fasta.append(fasta_1[i] + fasta_2[i])
    # write out the concat_alignment_fasta
    with open('{}/concat_alignment.fasta', 'w') as f:
        for line in concat_alignment_fasta:
            f.write('{}\n'.format(line))

    apples  = 'asdf'
    # run nucleotide model
    # make ML tree

def convert_interleaved_to_sequencial_fasta_two(fasta_in):
    fasta_out = []

    for i in range(len(fasta_in)):

        if fasta_in[i].startswith('>'):
            if fasta_out:
                # if the fasta is not empty then this is not the first
                fasta_out.append(temp_seq_str)
            #else then this is the first sequence and there is no need to add the seq.
            temp_seq_str = ''
            fasta_out.append(fasta_in[i])
        else:
            temp_seq_str = temp_seq_str + fasta_in[i]
    #finally we need to add in the last sequence
    fasta_out.append(temp_seq_str)
    return fasta_out

def crop_to_shortest_seq(mafft_aligned_fasta_path):
    # at this point we have the aligned .fasta written to the output directory
    # at this point we need to trim the fasta.
    # I was going to use trimAl but this doesn't actually have an option to clean up the ends of alignments
    # instead, read in the alignment as a TwoD list to a pandas dataframe
    # then delete the begining and end columns that contain gap sites
    aligned_fasta_interleaved = readDefinedFileToList(mafft_aligned_fasta_path)
    aligned_fasta = convert_interleaved_to_sequencial_fasta_two(aligned_fasta_interleaved)
    array_list = []
    for i in range(1, len(aligned_fasta), 2):
        array_list.append(list(aligned_fasta[i]))

    # make into pandas dataframe
    alignment_df = pd.DataFrame(array_list)

    # go from either end deleting any columns that have a gap
    columns_to_drop = []
    for i in list(alignment_df):
        # if there is a gap in the column at the beginning
        if '-' in list(alignment_df[i]) or '*' in list(alignment_df[i]):
            columns_to_drop.append(i)
        else:
            break
    for i in reversed(list(alignment_df)):
        # if there is a gap in the column at the end
        if '-' in list(alignment_df[i]) or '*' in list(alignment_df[i]):
            columns_to_drop.append(i)
        else:
            break

    # get a list that is the columns indices that we want to keep
    col_to_keep = [col_index for col_index in list(alignment_df) if col_index not in columns_to_drop]

    # drop the gap columns
    alignment_df = alignment_df[col_to_keep]

    # here we have the pandas dataframe with the gap columns removed
    # convert back to a fasta and write out
    cropped_fasta = []
    alignment_index_labels = list(alignment_df.index)
    for i in range(len(alignment_index_labels)):
        seq_name = '>{}'.format(alignment_index_labels[i])
        aa_seq = ''.join(alignment_df.loc[alignment_index_labels[i]])
        cropped_fasta.extend([seq_name, aa_seq])

    # here we have the cropped and aligned fasta
    # write it out
    aligned_cropped_fasta_path = mafft_aligned_fasta_path.replace('.fasta', '_cropped.fasta')
    with open(aligned_cropped_fasta_path, 'w') as f:
        for line in cropped_fasta:
            f.write('{}\n'.format(line))

    return aligned_cropped_fasta_path


def mothur_worker(input_queue, wkd):
    for smp in iter(input_queue.get, 'STOP'):
        # if the .names and .fasta file that we're trying to produce already exists
        # then we have already undertaken the mothur processing and we can move on
        one_name = os.path.exists('{}/{}/{}.1.trim.pcr.names'.format(wkd, smp, smp))
        two_name = os.path.exists('{}/{}/{}.2.trim.pcr.names'.format(wkd, smp, smp))
        if  one_name and two_name:
            continue
        else:
            # run the mothur
            subprocess.run(['mothur', '{}/{}/mothur_batch_1'.format(wkd, smp)])
            subprocess.run(['mothur', '{}/{}/mothur_batch_2'.format(wkd, smp)])

def clean_fasta_of_primer_match_info(path_to_dirty_fasta):
    #read dirty fasta to list
    with open(path_to_dirty_fasta, 'r') as f:
        dirty_fasta = [line.rstrip() for line in f]

    clean_fasta = [line.split('\t')[0] if line.startswith('>') else line for line in dirty_fasta]

    #write out clean fasta
    path_to_clean_fasta = path_to_dirty_fasta.replace('.fasta', '_clean.fasta')
    with open(path_to_clean_fasta, 'w') as f:
        for line in clean_fasta:
            f.write('{}\n'.format(line))
    return path_to_clean_fasta


def trouble_shooting():
    wkd = '/home/humebc/projects/SymPortalMS/making_example_data_set/SP_DEMO_DATA/psbA_processing_raw_fastq/Y04'

    # read in to list the pcr output fasta for _1
    with open('{}/Y04.1.trim.pcr.fasta'.format(wkd), 'r') as f:
        fasta_pcr_1 = [line.rstrip() for line in f]
    # use this to grab the sequences from the .trim fasta of _2
    with open('{}/Y04.2.trim.fasta'.format(wkd), 'r') as f:
        fasta_trim_2 = [line.rstrip() for line in f]
    psba_seq_names = [line.split('/')[0][1:] for line in fasta_pcr_1 if line.startswith('>')]
    # pull in the original _2 .fasta
    with open('{}/Y04.2.fasta'.format(wkd), 'r') as f:
        fasta_2 = [line.rstrip() for line in f]
    # fasta_2_dict = {fasta_2.split('\t')[0][1:] : fasta_2[i + 1] for i in range(len(fasta_2)) if fasta_2[i].startswith('>')}
    # pull in the original _2 qual as well
    with open('{}/Y04.2.qual'.format(wkd), 'r') as f:
        qual_2 = [line.rstrip() for line in f]

    _2_psba_fasta = []
    _2_psba_qual = []
    for i in range(len(fasta_2)):
        if fasta_2[i].startswith('>'):
            if fasta_2[i].split('/')[0][1:] in psba_seq_names:
                # then this is one of the seqs that worked for the _2
                _2_psba_fasta.extend([fasta_2[i], fasta_2[i + 1]])
    # populate new qual
    for i in range(len(qual_2)) -1:
        if qual_2[i].startswith('>'):
            if qual_2[i].split('/')[0][1:] in psba_seq_names:
                # then this is one of the seqs that worked for the _2
                _2_psba_qual.extend([qual_2[i], qual_2[i + 1]])



    # # write these out as their own fasta and look for the primer sequence
    # _2_psba_seqs_fasta = []
    # for i in range(len(fasta_trim_2)):
    #     if fasta_trim_2[i].startswith('>'):
    #         if fasta_trim_2[i].split('/')[0][1:] in psba_seq_names:
    #             # then this is one of the seqs that worked for the _2
    #             _2_psba_seqs_fasta.extend([fasta_trim_2[i], fasta_trim_2[i+1]])

    # now write out
    with open('{}/good_2.fasta'.format(wkd), 'w') as f:
        for line in _2_psba_seqs_fasta:
            f.write('{}\n'.format(line))

    apples = 'asdf'

def SP_MS_scatter_plot_ordination():
    # read in the .csv as pandas array
    plot_info_df = pd.read_csv('/home/humebc/phylogeneticSoftware/SymPortal_interim_250318/outputs/ordination/141/between_samples/mothur/C/PCoA_coords.csv', header=0, index_col=0)

    # colours
    colour_list = []
    for smp in list(plot_info_df.index)[:-1]:
        if str(smp).startswith('A'):
            colour_list.append('black')
        elif str(smp).startswith('C'):
            colour_list.append('green')
        else:
            colour_list.append('blue')

    # plot the first and second columns
    plt.scatter(plot_info_df['PC1'][:-1], plot_info_df['PC2'][:-1], c=colour_list)
    plt.xlabel('PC1 46%')
    plt.ylabel('PC2 39%')
    plt.text(0.1, -0.55, 'Acropora spp.', style='italic')
    plt.text(0.2, 0.1, 'Cyphastrea spp.', style='italic')
    plt.text(-0.38, 0.2, 'Platygyra spp.', style='italic')
    plt.text(plot_info_df['PC1']['A08'] + 0.025, plot_info_df['PC2']['A08'], 'A08')
    plt.savefig('/home/humebc/projects/SymPortalMS/ordination_110718.svg', format='svg')
    plt.show()

    return


def generate_supp_fig_1_for_sp_ms():
    """ This will generate the supp_fig_1 for the ms"""
    # This code will be used to start to investigate some database stats.
    # Firstly I will look for some light, high end stats for the PS MS i.e. number of samples, studies types etc.
    # I think it would also be nice to have a look at the frequency distribution of the types
    # I think that there is probably a very citable paper to be had from SP with regards to looking more deeply
    # at the data

    # get the lastest dataAnalysis
    latest_data_analysis_object = DataAnalysis.objects.get(id=44)
    how_many_data_sets = len(latest_data_analysis_object.list_of_data_set_uids.split(','))
    clade_collection_objects = list(latest_data_analysis_object.get_clade_collections())
    samples_of_clade_collections = list(DataSetSample.objects.filter(cladecollection__in=clade_collection_objects).distinct())
    # list of analysisTypes
    analysis_types_from_data_analysis_object = AnalysisType.objects.filter(data_analysis_from=latest_data_analysis_object)

    num_of_clade_col_types_associated_to_analysis_types = len(CladeCollectionType.objects.filter(analysis_type_of__in=analysis_types_from_data_analysis_object))

    analysis_type_abundance_list = []
    total_clade_col_types_counter = 0
    for analysis_type in analysis_types_from_data_analysis_object:
        analysis_type_abundance_list.append((analysis_type.name, len(CladeCollectionType.objects.filter(analysis_type_of=analysis_type))))

    sorted_analysis_type_abundance_list = sorted(analysis_type_abundance_list, key=lambda x: x[1], reverse=True)

    sorted_for_plotting = []
    count = 0
    for sat in sorted_analysis_type_abundance_list:
        total_clade_col_types_counter += sat[1]
        cum_value = total_clade_col_types_counter/num_of_clade_col_types_associated_to_analysis_types
        sorted_for_plotting.append((sat[0],cum_value))
        count += 1
        if cum_value > .5:
            print('cumulative_value=' + str(count))

    x_cumulative = range(len(sorted_for_plotting))
    y_cumulative = [a[1] for a in sorted_for_plotting]

    fig, ax1 = plt.subplots()
    ax1.scatter(x_cumulative, y_cumulative, c='black', marker='.')
    ax1.set_ylabel("Cumulative relative abundance as proportion\nof total ITS2 type profile instances")
    ax1.set_xlabel("Rank abundance of ITS2 type profiles")

    ax2 = ax1.twinx()

    # plt.scatter(x_cumulative, y_cumulative, c='black', marker='.')

    x_abs = range(len(sorted_analysis_type_abundance_list))
    y_abs = [x[1] for x in sorted_analysis_type_abundance_list]

    ax2.scatter(x_abs, y_abs, c='blue', marker='.')
    ax2.set_ylabel("Number of samples containing the ITS2 type profile")
    ax1.spines['right'].set_color('blue')
    # plt.scatter(x_abs, y_abs, c='blue', marker='.')
    fig.show()
    plt.savefig('supp_fig1.svg')
    plt.savefig('supp_fig1.png')
psba_phylogeny_from_scratch()