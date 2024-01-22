#! usr/bin/python3

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import argparse
import re

############## COMMAND LINE ARGS ##############

parser = argparse.ArgumentParser()
parser.add_argument('--outputFile','-o',default='UBE2C.pdf',type=str,action='store',help='output file goes here')
parser.add_argument('--pslFiles','-i',default='DMSO_H358.hg38.psl,AMG_H358.hg38.psl',type=str,action='store',help='path/`to/pslFile.psl')
parser.add_argument('--gtfFile','-g',default='gencode.v45.chr_patch_hapl_scaff.annotation.gtf',type=str,action='store',help='/path/to/genome.gtf')
parser.add_argument('--genoCoord','-c',default='chr20:45812644-45816952',type=str,action='store',help='Genomic Coordinates')
parser.add_argument('--rmskTable','-r',default='temp.table',type=str,action='store',help='rmsk Table')

args = parser.parse_args()

outFile=args.outputFile

############ PSL PLOTTING ################

class psl:
    def pslPlot(rawFile,panelObject,title,color):
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
        for y in range(0,len(middleLayerValues)):
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
        panelObject.set_ylim(0,len(middleLayerValues)+1.8)
        panelObject.set_xlabel(title,fontsize=6)

########### GTF PLOTTING #################

class gtf:
    def gtfPlot(rawFile,panelObject,title):
        gtfData = {}
        with open(rawFile,'r') as gtfFile:
            for line in gtfFile:
                if line.split('\t')[0] == chromosome and int(line.split('\t')[4]) > start and int(line.split('\t')[3]) < stop:
                    transcript = line.split('\t')[8].split(';')[1]
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
        layerValues = []
        while len(plotValues):
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
        panelObject.set_xlabel(title, fontsize = 6, weight = 'bold')
        panelObject.tick_params(labelleft=False)
        panelObject.set_xticks([start,stop])
        panelObject.set_xticklabels([start,stop],fontsize = 4,weight='bold')
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
    def rmskPlot(rmskTable,panelObject):
        rmskData = {}
        with open(rmskTable,'r') as rmskFile:
            next(rmskFile)
            for line in rmskFile:
                if line.split('\t')[5] == chromosome and int(line.split('\t')[6]) > start and int(line.split('\t')[7]) < stop:
                    repeat = line.split('\t')
                    rmskData.update({line.split('\t')[10]:repeat})

        repClassDict =  {'SINE':(Golden_brown, 0), 'LINE':(Wine, 1), 'LTR':(Jasmine, 2)} ### change colors later
        
        panelObject.set_xlim(float(start),float(stop))
        panelObject.set_ylim(0,1.8)
        panelObject.set_xlabel('Repeat Masker', fontsize = 6, weight = 'bold')
        panelObject.set_yticks([1,2,3])
        panelObject.set_yticklabels(['SINE','LINE','LTR'], fontsize = 6)
        panelObject.set_xticks([])
        #panelObject.set_title('Repeat Masker',fontsize = 6,weight='bold')
        for repName, repeat in rmskData.items():
            print(repeat[11])
            if repeat[11] in repClassDict.keys():
                y_axis  = repClassDict[repeat[11]][1]
                width = int(repeat[7])-int(repeat[6])
                temp = mplpatches.Rectangle((float(repeat[6]),(0.5+y_axis)),float(width),.75,
                        facecolor = repClassDict[repeat[11]][0] if repeat[11] in repClassDict.keys() else '#ceb1be',
                        edgecolor = 'black',
                        linewidth = 0.25)
                panelObject.add_patch(temp)


############# DEFINE FILES ###############

gtfFile = args.gtfFile
pslFile1 = args.pslFiles.split(',')[0]
pslFile2 = args.pslFiles.split(',')[1]
args.rmskTable


############## FIGURE FORMATTING ##############

figureHeight = 4.5
figureWidth = 8.5
panelHeight = 1.3
panelWidth = 7.5
plt.figure(figsize = (figureWidth,figureHeight))



############ COLORS ###############
Golden_brown = (153/255, 98/255, 30/255)
Jasmine = (255/255, 208/255, 123/255)
Wine = (109/255, 59/255, 71/255)
Lapis_lazuli = (35/255, 87/255, 137/255)
Cal_poly_green = (44/255, 85/255, 48/255)


########### PANEL SIZING AND RATIO #############

Gencode = plt.axes([.3/figureWidth,3/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
psl1 = plt.axes([.3/figureWidth,1.4/figureHeight,panelWidth/figureWidth,1/figureHeight])
psl2 = plt.axes([.3/figureWidth,0.281/figureHeight,panelWidth/figureWidth,1/figureHeight])
rmask = plt.axes([.3/figureWidth,2.2/figureHeight,panelWidth/figureWidth,0.5/figureHeight])

#Gencode = plt.axes([.2/figureWidth,3.5/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
#psl1 = plt.axes([.2/figureWidth,1.9/figureHeight,panelWidth/figureWidth,1/figureHeight])
#psl2 = plt.axes([.2/figureWidth,0.781/figureHeight,panelWidth/figureWidth,1/figureHeight])
#rmask = plt.axes([.2/figureWidth,2.7/figureHeight,panelWidth/figureWidth,0.5/figureHeight])

########## PARSE INPUT ##############

chromosome, start, stop = re.findall(r'(\w+):(\d+)-(\d+)', args.genoCoord)[0]
start = int(start) - 100
stop = int(stop) + 100

########## FILE PLOTTING ############

gtf.gtfPlot(args.gtfFile,Gencode,'GENCODE V45')
psl.pslPlot(args.pslFiles.split(',')[0],psl1,'DMSO-Treated',Lapis_lazuli)
psl.pslPlot(args.pslFiles.split(',')[1],psl2,'AMG-Treated',Cal_poly_green)
genomeBrowser.rmskPlot(args.rmskTable,rmask)


for panel in psl1,psl2:
    panel.tick_params(bottom=False, labelbottom=False,left=False, labelleft=False,
                   right=False, labelright=False,top=False, labeltop=False)
    for orientation in ['top','bottom','right','left']:
        panel.spines[orientation].set_visible(False)
for panel in Gencode,rmask:
    for orientation in ['top','bottom','right','left']:
            panel.spines[orientation].set_visible(False)
rmask.tick_params(bottom=False, labelbottom=False,left=False, labelleft=True,
                   right=False, labelright=False,top=False, labeltop=False)


plt.savefig(outFile,dpi = 2400)