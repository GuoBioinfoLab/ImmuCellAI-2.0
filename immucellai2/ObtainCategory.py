#!/usr/bin/python3
import pandas
import os
import re
import numpy
from pathlib import Path
from .myclasses import CLASS_FOR_RUNRESULT

def isExistspath(FileResult):
    try:
       abs_path = os.path.abspath(FileResult)
       Resultpath = os.path.dirname(abs_path)
    except Exception as e:
       print(f"Error in isExistspath: {e}")
       Resultpath = None
    return Resultpath

def ExtractResult(ResultObject, FileResult, ResultIndex = 0):
    isExistspath(FileResult)
    print("Save the result, wait...")
    ResultObject.save_result(FileResult, ResultIndex = ResultIndex )

def ExtractResult2(ResultObject, FileResult):
   print("Under achieved, and wait for update...")
   isExistspath(FileResult)
   pass

def ExtractResult3(ResultObject, FileResult):
   print("Under achieved, and wait for update...")
   isExistspath(FileResult)

class CLASS_CELLTYPENODE(object):
   def __init__(self):
      pass
   def mycopy(self):
      return {"CellType":"",
         "AlsoKnownAs": [],
         "RankedLayer": 10,
         "RelatedNode":{"HisParentNode":"", "HisChidNode":[]},
         "Description": "...",
         "CellTypeRatio": numpy.NAN } 

def ObtainInformation2(line, SearchString, DictValue):
   obcontent = line.replace(re.findall("\s*%s:\s*"%SearchString, line)[0], "")
   if isinstance(DictValue[SearchString], list):
      obcontent = obcontent.replace("[", "").replace("]","").replace("'","")
      while( ", " in obcontent): obcontent = obcontent.replace(", ", ",")
      DictValue[SearchString] += obcontent.split(",")
   else:
      DictValue[SearchString] =  obcontent

def ObtainCellTypeCategory(CellTypeCateogry):
    if CellTypeCateogry is None or not os.path.isfile(CellTypeCateogry):
        try:
            from importlib.resources import files
            CellTypeCateogry = str(files("immucellai2") / "myconfig" / "Celltype.category")
        except ImportError:
            import pkg_resources
            CellTypeCateogry = pkg_resources.resource_filename("immucellai2", "myconfig/Celltype.category")
    if not os.path.isfile(CellTypeCateogry):
        raise FileNotFoundError(
            f"Package resource 'myconfig/Celltype.category' not found at {CellTypeCateogry}"
        )
    CellTypeCateogryContent = {}
    with open(CellTypeCateogry, "r+") as CellTypeCateogryFile:
        oneCellTypeNode = ( CLASS_CELLTYPENODE() ).mycopy()
        for line in CellTypeCateogryFile.readlines():
            if "#[" not in line and "#" in line: line = line.split("#")[0]
            line = line.strip()
            if re.findall("^#\[NewNode\]", line) != []:
                if oneCellTypeNode["CellType"] != "":
                    CellTypeCateogryContent[oneCellTypeNode["CellType"]] = oneCellTypeNode
                oneCellTypeNode = (  CLASS_CELLTYPENODE() ).mycopy()
            for key1, value1 in oneCellTypeNode.items():
                if isinstance(value1, dict):
                    Break = False
                    for key2, value2 in value1.items():
                        if key2 in line:
                            ObtainInformation2(line, key2, value1)
                            Break = True
                            break
                    if Break: break
                else:
                    if key1 in line:
                       ObtainInformation2(line, key1, oneCellTypeNode)
                       break
        CellTypeCateogryContent[oneCellTypeNode["CellType"]] = oneCellTypeNode
        CellTypeCateogryFile.close()
    return CellTypeCateogryContent

def FindSummcelltyperatio(onesamplecelltyperatio,  
    currentNode, 
    knowratioCelltype, 
    CellTypeCateogryContent, 
    ):
    if currentNode in knowratioCelltype:
       return onesamplecelltyperatio[ currentNode ]
    elif not pandas.isnull(onesamplecelltyperatio[ currentNode ]):
       return onesamplecelltyperatio[ currentNode ]
    else:
       try:
          ChildNodes = CellTypeCateogryContent[currentNode]["RelatedNode"]["HisChidNode"]
       except:
          ChildNodes = []
       celltyperatiosum = 0.00
       for oneChildNode in ChildNodes:
          if oneChildNode == "": continue
          celltyperatiosum += FindSummcelltyperatio(onesamplecelltyperatio ,
             oneChildNode ,
             knowratioCelltype ,
             CellTypeCateogryContent)
       onesamplecelltyperatio[ currentNode ] = celltyperatiosum
       return celltyperatiosum
   
def SummaryCellRatio(ResultObject,
   #CellTypeCateogry = "myconfig/Celltype.cateogory",
   FileResult = "myresult/SummaryCellRatio.txt",
   ResultIndex = 0):
   print("Obtain the Celltype ratio by summary the Cell subtype...")
   isExistspath(FileResult)
   CellTypeCateogry = ResultObject.FileCellTypeCategory
   #CellTypeCateogry = re.findall('^.*/', os.path.abspath(sys.argv[0]) )[0] + CellTypeCateogry
   CellTypeCateogryContent = ObtainCellTypeCateogry(CellTypeCateogry)
   AllCellType = []
   for key, oneCellTypeNode in CellTypeCateogryContent.items():
       AllCellType += [key,
          oneCellTypeNode["RelatedNode"]["HisParentNode"]] + \
          oneCellTypeNode["RelatedNode"]["HisChidNode"]
   AllCellType = list(set(AllCellType))
   SampleName, CellType = ResultObject.SampleName, ResultObject.CellType
   SumCellTypeRatio = pandas.DataFrame(columns = AllCellType, index= SampleName )
   AllsampleCellTypeRatio = ResultObject.get_result(ResultIndex = ResultIndex )
   knowratioCelltype = []
   for Celltypeone in CellType:
      locatecelltype = ""
      for key, oneCellTypeNode in CellTypeCateogryContent.items():
         if Celltypeone in [ key ] + oneCellTypeNode["AlsoKnownAs"]:
            locatecelltype = key
         elif Celltypeone in oneCellTypeNode["RelatedNode"]["HisChidNode"]:
            locatecelltype = Celltypeone
         if locatecelltype != "":
            knowratioCelltype.append(locatecelltype)
            for onesample in SampleName:
               SumCellTypeRatio.loc[onesample, locatecelltype ] = \
                  AllsampleCellTypeRatio.loc[onesample, Celltypeone]
            break
   print("known ratio celltype:", knowratioCelltype)
   for onesample in SampleName:
      print("sample name:",  onesample)
      summtable = []
      for keyword, oneCellTypeNode in CellTypeCateogryContent.items():
         if oneCellTypeNode["RankedLayer"] == '1':
            FindSummcelltyperatio(
               onesamplecelltyperatio = SumCellTypeRatio.loc[onesample, :],
               currentNode = keyword,
               knowratioCelltype = knowratioCelltype,
               CellTypeCateogryContent = CellTypeCateogryContent
            )
            summtable += [ keyword+": %f"%SumCellTypeRatio.loc[onesample, keyword ] ]
      print( summtable )
   ResultObject.SumCellTypeRatio = SumCellTypeRatio
   ResultObject.CellTypeCateogryContent = CellTypeCateogryContent
   SumCellTypeRatio.to_csv(FileResult, sep="\t")
   #SumCellTypeRatio.to_csv(FileResult, index=False, header=False, sep="\t")
   #return CellTypeCateogryContent, SumCellTypeRatio

   
