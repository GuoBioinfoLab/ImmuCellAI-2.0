from .Deconvolution import MainRun # MainRun
from .PrepareData import PrepareData # PrepareDate
from .ObtainCategory import ExtractResult, SummaryCellRatio # ExtractResult()
from .Drawplot import DrawPlotFunction1, DrawPlotFunction2 # visualization()
import os
import pandas
import argparse
import time
import multiprocessing as mp
import numpy as np
import random

def load_tumor_reference_data():
    from importlib.resources import path
    with path("immucellai2.myconfig", "reference_tumorCelltypes.txt") as file_path:
        return pandas.read_table(file_path, sep="\t", header=0, index_col=0)
        
def load_normal_reference_data():
    from importlib.resources import path
    with path("immucellai2.myconfig", "reference_normalCelltypes.txt") as file_path:
        return pandas.read_table(file_path, sep="\t", header=0, index_col=0)

def run_ImmuCellAI2(reference_file, sample_file, output_file, thread_num=16, seed=42):
    run_obj = PrepareData(
        FileReferenceProfile=reference_file,
        FileSampleExpressionProfile=sample_file,
        EnvironmentConfig=thread_num,
        InitialCellTypeRatio=('Normone', 'randn')
    )
    result_obj = MainRun(run_obj, seed=seed)
    ExtractResult(result_obj, output_file, ResultIndex=0)
    return result_obj

#DrawPlotFunction1(result_obj)

def main():
   seed = 43
   np.random.seed(seed) 
   random.seed(seed)
   start_time = time.time()  
   mp.set_start_method('fork', force=True) 
   parser = argparse.ArgumentParser(description='ImmuCellAI2 deconvolution tool')
   parser.add_argument('-f', '--reference', dest = 'reference', default="", help = 'celltype reference experession matrix')
   parser.add_argument('-g', '--genes', dest = 'CoveredGenes', action = 'store_const', help = 'The genes used in the following deconvolution process selected by BayesPrism, temporary variance')
   parser.add_argument('-s', '--sample', dest = 'sample', required = True, help = 'the samples gene expression profile')
   parser.add_argument('-t', '--thread', dest = 'thread', type=int, default = 16, help = "threading numbers for deconvolution")
   parser.add_argument('--sample-type', dest="sample_type", default = 'normal', choices=['tumor','normal'], help='Sample type (tumor/normal, default normal)')
   parser.add_argument('-o', '--output', dest = "output", default = "myresult/ResultDeconvolution.xlsx", help = " the path/filename to save the deconvaluted result.")
   parser.add_argument('--seed', type=int, default=42, help='Random seed')  
   args = parser.parse_args() 
   if args.reference:
      reference_data = pandas.read_csv(args.reference, sep='\t')
   else:
      reference_data = load_tumor_reference_data() if args.sample_type == 'normal' else load_normal_reference_data()
   
   print("### Begin run deconvolution tools, wait....") 
   return run_ImmuCellAI2(
        reference_file=reference_data,
        sample_file=args.sample,
        output_file=args.output,
        thread_num=args.thread,
        seed=args.seed
    )
 
if __name__ == "__main__":
   main() 
