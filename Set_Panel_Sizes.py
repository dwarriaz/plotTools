#! usr/bin/python3

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import argparse
import re

############## COMMAND LINE ARGS ##############

parser = argparse.ArgumentParser()
parser.add_argument('--outputFile','-o',default='output/UBE2C.pdf',type=str,action='store',help='output file goes here')
parser.add_argument('--pslFiles','-i',default='data/DMSO_H358.hg38.psl,data/AMG_H358.hg38.psl',type=str,action='store',help='path/`to/pslFile.psl')
parser.add_argument('--gtfFile','-g',default='data/gencode.v45.chr_patch_hapl_scaff.annotation.gtf',type=str,action='store',help='/path/to/genome.gtf')
parser.add_argument('--genoCoord','-c',default='chr20:45812644-45816952',type=str,action='store',help='Genomic Coordinates')
parser.add_argument('--rmskTable','-r',default='data/rmsk.table',type=str,action='store',help='rmsk Table')
parser.add_argument('--gene','-G',default='',type=str,action='store',help='Gene of interest')


args = parser.parse_args()

outFile=args.outputFile

############ PSL PLOTTING ################

class psl:
    def pslData(rawFile):
        with open(rawFile,'r') as pslFile:

            rawdata = []

            for line in pslFile:
                if line.split('\t')[13] == chromosome and int(line.split('\t')[15]) > start and int(line.split('\t')[16]) <  stop:
                    data = line.strip('\n').split('\t')
                    rawdata.append([int(data[15]),int(data[16]),int(data[17]),data[20],data[18]])
                    rawdata.sort(key = lambda x:x[0])

        middleLayerValues = []
        while len(rawdata):
            currentLayer = []
            for alignment in rawdata:
                if len(currentLayer) == 0:
                    currentLayer.append(alignment)
                else:
                    if alignment[0] > currentLayer[-1][1]:
                        currentLayer.append(alignment)
            for alignment in currentLayer:
                rawdata.remove(alignment)
            middleLayerValues.append(currentLayer)
        
        return middleLayerValues
        

    def pslPlot(middleLayerValues,panelObject,color):    
        for y in range(0,len(middleLayerValues)):
            if y <= 12:
                for data in middleLayerValues[y]:
                    first = mplpatches.Rectangle((data[0],y+0.25),(data[1]-data[0]),0.1,
                                                facecolor = color,
                                                edgecolor = 'black',
                                                linewidth = 0.25)
                    panelObject.add_patch(first)
                    blockStart = data[3].split(',')
                    blockwidth = data[4].split(',')
                    for x in range(0,data[2]):
                        second = mplpatches.Rectangle((int(blockStart[x]),y),int(blockwidth[x]),0.5,
                                                    facecolor = color,
                                                    edgecolor = 'black',
                                                    linewidth = 0.25)
                        panelObject.add_patch(second)
            panelObject.set_xlim(float(start),float(stop))
            minimum_Y = min([12,len(middleLayerValues)])
            panelObject.set_ylim(0,minimum_Y+1.8)
            #panelObject.set_xlabel(title,fontsize=6)

########### GTF PLOTTING #################

class gtf:
    def gtfDataProcessing(rawFile,gene=''):
        gtfData = {}
        with open(rawFile,'r') as gtfFile:
            for line in gtfFile:
                if line.split('\t')[0] == chromosome and int(line.split('\t')[4]) > start and int(line.split('\t')[3]) < stop:
                    transcript = line.split('\t')[8].split(';')[1]
                    if len(gene) != 0 and line.split('\t')[1] == gene:
                        if transcript in gtfData.keys():
                            gtfData[transcript].append(line.split('\t'))
                        else:
                            gtfData.update({transcript:[line.split('\t')]})
                    else:
                        if transcript in gtfData.keys():
                            gtfData[transcript].append(line.split('\t'))
                        else:
                            gtfData.update({transcript:[line.split('\t')]})

        widthDict =  {'transcript':0.05, 'exon':0.25, 'CDS':0.5}

        plotValues = []

        for key,value in gtfData.items():
            transcript = []
            for line in value:
                if line[2] in widthDict.keys():
                    transcript.append((line[3],line[4],line[2]))
            if len(transcript) != 0:
                plotValues.append(transcript)
        plotValues.sort(key = lambda x:x[0][0])

        return plotValues

    def gtfPlot(plotValues,panelObject):
        layerValues = []
        widthDict =  {'transcript':0.05, 'exon':0.25, 'CDS':0.5}
        while len(plotValues) and len(layerValues) <= 8:
            currentLayer = []
            for transcript in plotValues:
                if len(currentLayer) == 0:
                    currentLayer.append(transcript)
                else:
                    if transcript[0][0] > currentLayer[-1][0][1]:
                        currentLayer.append(transcript)
            for tx in currentLayer:
                plotValues.remove(tx)
            layerValues.append(currentLayer)
        panelObject.set_xlim(float(start),float(stop))
        panelObject.set_ylim(0,len(layerValues)+1.8)
        #panelObject.set_xlabel(title, fontsize = 6, weight = 'bold')
        panelObject.tick_params(labelleft=False,size=1)
        xticks = [xticks for xticks in range(start,stop,1000)]
        xticks.append(stop)
        panelObject.set_xticks(xticks)
        panelObject.set_xticklabels(xticks,fontsize = 2,weight='bold')
        panelObject.set_yticks([])
        for y in range(0,len(layerValues)):
            for transcript in layerValues[y]:
                for exon in transcript:
                    if exon[2] in widthDict.keys():
                        width = int(exon[1])-int(exon[0])
                        temp = mplpatches.Rectangle((float(exon[0]),(1+y-float(widthDict[exon[2]]/2))),float(width),float(widthDict[exon[2]]),
                                                    facecolor = 'gray',
                                                    edgecolor = 'black',
                                                    linewidth = 0.25)
                        panelObject.add_patch(temp)

########### RMSK PLOTTING ################

class genomeBrowser:
    def rmskData(rmskTable):
        rmskData = {}
        with open(rmskTable,'r') as rmskFile:
            next(rmskFile)
            for line in rmskFile:
                if line.split('\t')[5] == chromosome and int(line.split('\t')[6]) >= start and int(line.split('\t')[7]) <= stop:
                    repeat = line.split('\t')
                    rmskData.update({line.split('\t')[10]:repeat})

        repClassDict =  {'SINE':[Golden_brown, 0], 'LINE':[Wine, 1], 'LTR':[Jasmine, 2]} ### change colors later
        
        layer_limit = set()
        for repName, repeat in rmskData.items():
            if repeat[11] in repClassDict.keys():
                layer_limit.add(repClassDict[repeat[11]][1])
                y_axis = repClassDict[repeat[11]][1]

        return (max(list(layer_limit)),rmskData)

    def rmskPlot(rmskData,panelObject,y_max): 

        repClassDict =  {'SINE':[Golden_brown, 0], 'LINE':[Wine, 1], 'LTR':[Jasmine, 2]}
        repClassLabels = ['SINE','LINE','LTR']

        panelObject.set_xlim(float(start),float(stop))
        panelObject.set_yticks([x+.5 for x in range(0,y_max)])
        panelObject.set_yticklabels([repClassLabels[x] for x in range(0,y_max)], fontsize = 2, weight='bold')
        panelObject.set_xticks([])

        
        layer_limit = set()
        for repName, repeat in rmskData.items():
            print(repName)
            if repeat[11] in repClassDict.keys():
                layer_limit.add(repClassDict[repeat[11]][1])
                y_axis = repClassDict[repeat[11]][1]       
                width = int(repeat[7])-int(repeat[6])
                temp = mplpatches.Rectangle((float(repeat[6]),y_axis+0.02),float(width),.75,
                        facecolor = repClassDict[repeat[11]][0] if repeat[11] in repClassDict.keys() else '#ceb1be',
                        edgecolor = 'black',
                        linewidth = 0.25)
                panelObject.add_patch(temp)
        panelObject.set_ylim(0,max(list(layer_limit))+1)


############# DEFINE FILES ###############


pslFile1 = args.pslFiles.split(',')[0]
pslFile2 = args.pslFiles.split(',')[1]


########## PARSE INPUT ##############

chromosome, start, stop = re.findall(r'(\w+):(\d+)-(\d+)', args.genoCoord)[0]
start = int(start) - 100
stop = int(stop) + 100

############ COLORS ###############
Golden_brown = (153/255, 98/255, 30/255)
Jasmine = (255/255, 208/255, 123/255)
Wine = (109/255, 59/255, 71/255)
Lapis_lazuli = (35/255, 87/255, 137/255)
Cal_poly_green = (44/255, 85/255, 48/255)


######### DETERMINE RATIOS ##########

gtfPlottingData = gtf.gtfDataProcessing(args.gtfFile,args.gene)
psl1PlottingData = psl.pslData(pslFile1)
psl2PlottingData = psl.pslData(pslFile2)
rmskPlottingData = genomeBrowser.rmskData(args.rmskTable)
print(rmskPlottingData)

total = len(gtfPlottingData)+len(psl1PlottingData)+len(psl2PlottingData)+rmskPlottingData[0]+1


############## FIGURE FORMATTING ##############

figureHeight = 1
figureWidth = 6.5
panelWidth = 6
plt.figure(figsize = (figureWidth,figureHeight))

########### PANEL SIZING AND RATIO #############

psl2StackHeight = 0.02
psl1StackHeight = psl2StackHeight + (len(psl2PlottingData)/total)
rmskStackHeight = psl1StackHeight + (len(psl1PlottingData)/total) - .07
gtfStackHeight = rmskStackHeight + ((rmskPlottingData[0]+1)/total) + .12

psl2 = plt.axes([0.04615384615384615, 0.02, 0.9230769230769231, 0.29411764705882354])
psl1 = plt.axes([0.04615384615384615, 0.31411764705882356, 0.9230769230769231, 0.17647058823529413])
rmask = plt.axes([0.04615384615384615, 0.44058823529411765, 0.9230769230769231, 0.058823529411764705])
Gencode = plt.axes([0.04615384615384615, 0.5994117647058823, 0.9230769230769231, 0.32454361054766734])


#Gencode = plt.axes([.2/figureWidth,3.5/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
#psl1 = plt.axes([.2/figureWidth,1.9/figureHeight,panelWidth/figureWidth,1/figureHeight])
#psl2 = plt.axes([.2/figureWidth,0.781/figureHeight,panelWidth/figureWidth,1/figureHeight])
#rmask = plt.axes([.2/figureWidth,2.7/figureHeight,panelWidth/figureWidth,0.5/figureHeight])


########## FILE PLOTTING ############

gtfPlottingData = gtf.gtfDataProcessing(args.gtfFile,args.gene)

gtf.gtfPlot(gtfPlottingData,Gencode)
psl.pslPlot(psl1PlottingData,psl1,Lapis_lazuli)
psl.pslPlot(psl2PlottingData,psl2,Cal_poly_green)
genomeBrowser.rmskPlot(rmskPlottingData[1],rmask,rmskPlottingData[0]+1)


for panel in psl1,psl2:
    panel.tick_params(bottom=False, labelbottom=False,left=False, labelleft=False,
                   right=False, labelright=False,top=False, labeltop=False)
    panel.margins(x=0)
    for orientation in ['top','bottom','right','left']:
        panel.spines[orientation].set_visible(False)
for panel in Gencode,rmask:
    for orientation in ['top','bottom','right','left']:
            panel.spines[orientation].set_visible(False)
    panel.margins(x=0)
rmask.tick_params(bottom=False, labelbottom=False,left=False, labelleft=True,
                   right=False, labelright=False,top=False, labeltop=False)


plt.savefig(outFile,dpi = 2400)