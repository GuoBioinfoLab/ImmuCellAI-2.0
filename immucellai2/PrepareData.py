#!/usr/bin/python3
from .myclasses import CLASS_FOR_RUN
from .ObtainCategory import ObtainCellTypeCategory
import pandas
import os
import importlib.util
import re
import sys
import multiprocessing as mp

def Obtainmyconfigpath():
   script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
   return os.path.join(script_dir, "immucellai2", "myconfig", "")

def get_package_dir(package_name):
    spec = importlib.util.find_spec(package_name)
    if spec is None:
        raise ModuleNotFoundError(f"Package '{package_name}' not found")
    path = spec.origin or spec.submodule_search_locations[0]
    return os.path.dirname(path) 

def SelectGeneForDeconvolution(DFReferenceProfile, FileCoveredGenes="", Method="UsedMarker"):
    print("Select the gene for the following deconvolution...")
    GeneUsedForDeconvolution = []
    DFReferenceProfileGenes = DFReferenceProfile.index.values
    if Method == "UsedMarker":
        if FileCoveredGenes == "":
            pkg_marker_path = os.path.join(get_package_dir("immucellai2"), "myconfig/MarkerUsedDeconvolution.txt")
            if os.path.exists(pkg_marker_path):
                FileCoveredGenes = pkg_marker_path
                print(f"[INFO] Found marker file in immucellai2 package: {FileCoveredGenes}")
            GeneUsedForDeconvolution0 = (pandas.read_table(FileCoveredGenes, sep= "\t")).iloc[0].to_list()
            GeneUsedForDeconvolution = list(set(GeneUsedForDeconvolution0).intersection(set(DFReferenceProfileGenes)))
    return GeneUsedForDeconvolution

def CelltypeCategoryCheck(FileCellTypeCategory = "", celltypelist = [] ):
   print("Check the Celltype covered by configfile")
   if FileCellTypeCategory == "":
      try:
         with resources.path("immucellai2.myconfig", "Celltype.category") as config_path:
            FileCellTypeCategory = str(config_path)
      except Exception as e:
         FileCellTypeCategory = Obtainmyconfigpath() + "Celltype.category"
         print(f"[WARNING] Using fallback path for config file: {FileCellTypeCategory}")
   try:
      obtaincontent = ObtainCellTypeCateogry(FileCellTypeCategory)
   except Exception as e:
      print(f"Error reading cell type config: {str(e)}")
      raise
   Allcelltype = []
   for keyword, oneCellTypeNode in obtaincontent.items():
      Allcelltype += [ keyword ] + oneCellTypeNode["AlsoKnownAs"] + oneCellTypeNode["RelatedNode"]["HisChidNode"]
   for onecelltype in celltypelist:
      if onecelltype not in Allcelltype:
         raise ValueError( "EEROR: reference matrix celltpe'{0}' NOT IN configfile, please CHECK...".format(onecelltype))
   return FileCellTypeCategory
   
def InitialCellTypeRatioCheck(InitialCellTypeRatio, FileInitialCellTypeRatio = "", ncelltype = 0):
   print("Check the celltype ratio initialization method...")
   if InitialCellTypeRatio[1] != "prior":
      return
   if FileInitialCellTypeRatio == "":
      FileInitialCellTypeRatio = Obtainmyconfigpath() + "myCellTypeRatio.initial"
   Expactedcelltypenum = (pandas.read(FileInitialCellTypeRatio, sep = "\t", header = 0, index_col = 0)).shape[1]
   if Expactedcelltypenum <1:
      raise ValueError("FAILED")
   elif Expactedcelltypenum in [ ncelltype, ncelltype -1 ]:
      return FileInitialCellTypeRatio
   else:
      InitialCellTypeRatio = 'randn'     

def PrepareData(FileReferenceProfile , 
   FileSampleExpressionProfile , 
   EnvironmentConfig = ("", "") ,
   FileCoveredGenes = "" ,
   FileCellTypeCategory = "" ,
   FileInitialCellTypeRatio = "" ,
   InitialCellTypeRatio = ('Normone', 'randn')):
   print("prepare for RunObject...")
   if FileReferenceProfile.shape[1] < 2:
      print("warning: When open Reference File, might sep = ' ' not '\t'")
   print("celltype reference raw matrix:\n", FileReferenceProfile.iloc[0:4, 0:4])
   ReferenceCelltype = {}
   for oneCellType in FileReferenceProfile.columns.values.tolist():
      numbertail = re.findall("\.[0-9]*$", oneCellType)
      oneCellType0 = oneCellType
      if numbertail != []: oneCellType = oneCellType[:-len(numbertail)]
      if oneCellType in ReferenceCelltype.keys():
         ReferenceCelltype[oneCellType].append(ReferenceCelltype[oneCellType])
      else: ReferenceCelltype[oneCellType] = [oneCellType0]
   DFReferenceProfile = pandas.DataFrame(columns = list(ReferenceCelltype.keys()),
       index = FileReferenceProfile.index.values)
   for celltype in  DFReferenceProfile.columns.values:
        DFReferenceProfile[celltype] = (
           FileReferenceProfile.loc[:, ReferenceCelltype[celltype] ]).mean(axis = 1)
   print("celltype reference matrix:\n", DFReferenceProfile.iloc[0:4, 0:4]) 
   DFSampleExpressionProfile = pandas.read_table(FileSampleExpressionProfile, sep = "\t", header = 0, index_col = 0)
   print(" initialize a Object For running...")  
   print("environment config(cpus, threads): ", EnvironmentConfig)
   GeneUsedForDeconvolution = SelectGeneForDeconvolution(DFReferenceProfile)
   #FileCellTypeCategory = CelltypeCategoryCheck(FileCellTypeCategory, celltypelist = list(ReferenceCelltype.keys()))
   FileInitialCellTypeRatio = InitialCellTypeRatioCheck(InitialCellTypeRatio, 
      FileInitialCellTypeRatio, ncelltype = DFReferenceProfile.shape[1]) 
   DFReferenceProfile0 = DFReferenceProfile.loc[GeneUsedForDeconvolution, ]
   DFReferenceProfile0 = DFReferenceProfile0[DFReferenceProfile0.index.isin(DFSampleExpressionProfile.index)]   
   selected_DFSampleExpressionProfile = DFSampleExpressionProfile.loc[DFReferenceProfile0.index]
   selected_DFSampleExpressionProfile = selected_DFSampleExpressionProfile.transpose() 
   SampleList = list(selected_DFSampleExpressionProfile.index) 
   return CLASS_FOR_RUN(
      DFReferenceProfile0, 
      selected_DFSampleExpressionProfile, 
      SampleList,
      EnvironmentConfig,
      InitialCellTypeRatio = InitialCellTypeRatio,
      FileCellTypeCategory = FileCellTypeCategory,
      FileInitialCellTypeRatio = FileInitialCellTypeRatio,) 
