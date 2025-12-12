#!/usr/bin/python3
import matplotlib.pyplot as plt
from .myclasses import CLASS_FOR_RUNRESULT
import numpy
import pandas
import random

def Findtargetedsample(SampleName, sample):
    if isinstance(sample, int):
        if 0 <= sample < len(SampleName):
            print("Show '%s' Sampling Result" % SampleName[sample])
            return sample
        else:
            raise IndexError("The sample index %d is out of range. Valid range is 0 to %d." % (sample, len(SampleName)-1))
    elif isinstance(sample, str):
        if sample in SampleName:
            print("Show '%s' Sampling Result" % sample)
            return SampleName.index(sample)
        else:
            raise ValueError("The sample name '%s' is not found in the sample list." % sample)
    else:
        raise TypeError("The 'sample' argument should be an integer index or a string name.")


def DrawPlotFunction1(ResultObject, sample = 0, PlotSavePath = ''):
   sample = Findtargetedsample(ResultObject.SampleName, sample)
   CellTypeLabels = ResultObject.CellType
   #mcmc_samples = (ResultObject.FlatSamplingResult)[sample]
   mcmc_samples = (ResultObject.McmcSamplingResult)[sample, :, :, :]
   if PlotSavePath == "":
      DrawPlotFunction1Function1(mcmc_samples = mcmc_samples, CellTypeLabels = CellTypeLabels )
   else :
      DrawPlotFunction1Function1(
         mcmc_samples = mcmc_samples, 
         CellTypeLabels = CellTypeLabels,
         PlotSavePath = PlotSavePath )

def DrawPlotFunction1Function1(
   mcmc_samples, 
   CellTypeLabels,
   PlotSavePath = "/public/home/yangjingmin/immucellAI2.0/backup0106/learning/myresult/SampingByEmceeShowingIterationProcess.png"
   ):
   print("Set parameters 'PlotSavePath='', Or default path/filename will used...")
   fig, axes = plt.subplots( len(CellTypeLabels), sharex=True, figsize=(15, 3 * len(CellTypeLabels)))
   for ii in range(len(CellTypeLabels)):
      ax = axes[ii]
      ax.plot(mcmc_samples[:, :, ii], "k", alpha =0.3, rasterized=True)
      xlength = mcmc_samples.shape[0]
      ax.set_xlim( 0, xlength )
      ax.set_ylabel(CellTypeLabels[ii])
   axes[-1].set_xlabel("step number")
   plt.savefig(PlotSavePath) 



def MpltDrawCircle(
   axes,
   center = (0, 0), 
   radius = 2, 
   thin = 0.2,
   CircleRange = (0, 1),
   color = "lightblue",
   label = "xxx" ):
   xcenter, ycenter = center[0], center[1]
   theta = numpy.arange( CircleRange[0] * 2 * numpy.pi, (CircleRange[0] + CircleRange[1]) * 2 * numpy.pi, 0.01) 
   radiuslist =  numpy.arange(radius, radius + thin, 0.01)
   
   #xxx = numpy.concatenate( xcenter + radiuslist[:, None] * numpy.cos(theta) )
   #yyy = numpy.concatenate( ycenter + radiuslist[:, None] * numpy.sin(theta) )
   xxx, yyy =[], []
   for iii in range(radiuslist.shape[0]):
     xxxiii = (xcenter + radiuslist[iii] * numpy.cos(theta)).tolist()
     yyyiii = (ycenter + radiuslist[iii] * numpy.sin(theta)).tolist()
     if iii%2 == 1:
        xxx += xxxiii[::-1]
        yyy += yyyiii[::-1]
     else: 
        xxx += xxxiii
        yyy += yyyiii
   axes.plot(xxx, yyy, color = color, label = label )

def randomcolor():
   colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
   color ="#"+''.join([random.choice(colorArr) for i in range(6)])
   return color

def DrawPlotFunction2(
    ResultObject,
    sample = 0,
    PlotSavePath = "myresult/CellTypeClassificationRatioCircle.png"):

    sampleint = Findtargetedsample(ResultObject.SampleName, sample)
    SumCellTypeRatio = ResultObject.SumCellTypeRatio
    SelectSample = SumCellTypeRatio.iloc[sampleint, :]
    CellTypeCateogryContent = ResultObject.CellTypeCateogryContent
    #plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.rcParams['axes.unicode_minus'] = False

    figure = plt.figure()
    axes = figure.add_subplot(111)
    radius, thin = 0.5, 1
    designedlayer = [ 0 for ii in range(8) ]
    for keyword, oneCellTypeNode in CellTypeCateogryContent.items():

       celllayer = int( oneCellTypeNode["RankedLayer"] )
       if pandas.isnull(SelectSample[keyword]):continue
       if celllayer <1 : continue
 
       #color = [ tuple(numpy.random.choice(range(0, 255), size=3)) ]
       MpltDrawCircle(axes, 
          radius = radius + ( celllayer - 1 ) *0.1 + sum([(thin - 0.2 *ii) for ii in range(celllayer-1) ]), 
          thin = thin - 0.2 *(celllayer -1),
          CircleRange = (designedlayer[celllayer -1], SelectSample[keyword]),
          #color = color[0].
          color = randomcolor(),
          label = keyword)
       designedlayer[celllayer -1] += SelectSample[keyword]        
   
    axes.yaxis.set_major_locator(plt.NullLocator())
    axes.xaxis.set_major_locator(plt.NullLocator())
    axes.spines['top'].set_color('none')
    axes.spines['bottom'].set_color('none')
    axes.spines['left'].set_color('none')
    axes.spines['right'].set_color('none')

    axes.axis("equal")
    
    plt.legend(loc = 'upper right', ncol=1, bbox_to_anchor=(1, 1), prop={ 'weight': 'bold', 'size': 6} )
    plt.title(" Celltype ratio")
    #plt.show()
    plt.savefig(PlotSavePath)


def DrawPlotFunction3():
   print("Under achievement...., wait for upate")

