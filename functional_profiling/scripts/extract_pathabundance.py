#!/usr/bin/env python3

#from matplotlib import colors
from numpy import inner
import pandas as pd
import argparse

"""
This code is to extract the pathways from the output of the Humann2
"""

def extract_pathways():
    input_file=pd.read_table(args.input, header=0)
    pa=input_file[~input_file["ID"].str.contains('\|')] #Filtering off the pathways abundance at species level
    pa.columns=pa.columns.str.replace("_pathabundance","") #renaming columnames by replacing "_pathabundance"
    pa=pa.drop(columns=['GA101o','GN101o','NEG','GP101o']) #droping the blank sample(neg), old sequenced samples:'GA101o','GN101o','GP101o'
    
    if args.samples:
        pa=pa.drop(columns=args.samples) 
    
    pa=pa.loc[pa['ID']!='UNMAPPED']
    pa=pa.loc[pa['ID']!='UNINTEGRATED']
    pa.to_csv('pathways_mod.csv',index=False)
    return pa

def topN_pa_across_samples(pa_samples):
    pa_samples.set_index('ID',inplace=True)
    pa_samples=pa_samples.assign(total_abundance=lambda x: x.sum(axis=1))
    pa_samples_sort=pa_samples.sort_values(by=['total_abundance'],ascending=False) #Sort the samples based total abundance in decreasing order
    pa_samples_sort.to_csv(f'{args.output}-sorted.csv')
    pa_samples_sort=pa_samples_sort.drop(columns=['total_abundance']) #Dropping column total_abundance. Since its purpose is served. Not needed in future.
    
    # Select the top N rows from the sorted dataframe
    topNpa_samples_abs=pa_samples_sort.head(args.topN)
    
    #Estimation of others pathways abundance within a sample and adding to the topNpa_samples_abs.
    topNpa_samples_abs=topNpa_samples_abs.T #transpose the topN pathways dataframe
    #Absoulute abundance of pathways in each sample.
    topNpa_samples_abs['Others']=pa_samples_sort.iloc[args.topN::].sum(axis=0)# use the same sorted dataframe to estimate the abundance of others excluding tonp N pathways
    
    #Relative abundance of pathways within each sample.
    topNpa_samples_rel=topNpa_samples_abs.div(topNpa_samples_abs.sum(axis=1), axis=0) #Estimating the relative abundance of each genefamily within a sample.
    
    return topNpa_samples_abs,topNpa_samples_rel


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__ , formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i','--input',type=str, required=True, help="Pathway abundance summary file across multiple samples")
    parser.add_argument('-o','--output',type=str, help="Output file")
    parser.add_argument('-m','--metadata', type=str, help="metadata about the samples in csv format with SampleId containing the sample names")
    parser.add_argument('-s','--samples', nargs='*', help="list of samples(from columns) to discard from table eg -s GA101o GN101o NEG GP101o")
    parser.add_argument('-c','--cellpopn', nargs='*', help="list of samples containing (S and/or N and/or P) to discard from table eg -s A N")
    parser.add_argument('-t','--topN', type=int, help="top N pathways")
    args = parser.parse_args()

pathways=extract_pathways()
topNpa_abs,topNpa_rel=topN_pa_across_samples(pathways)

if args.metadata:
    metadata=pd.read_csv(args.metadata,header=0,index_col="SampleId")
    topNpa_abs_meta=pd.concat([metadata,topNpa_abs],axis=1,join="inner")
    topNpa_rel_meta=pd.concat([metadata,topNpa_rel],axis=1,join="inner")

    #Write absoulute and relative abundance with metadata
    topNpa_abs_meta.to_csv(f'{args.output}-abs.csv')
    topNpa_rel_meta.to_csv(f'{args.output}-rel.csv')

else:
    #Write absoulute and relative abundance without metadata
    topNpa_abs.to_csv(f'{args.output}-abs.csv')
    topNpa_rel.to_csv(f'{args.output}-rel.csv')

