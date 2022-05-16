#!/usr/bin/env python3
#for: humann.nf

import pandas as pd
import argparse

"""
This code is to extract the pathways from the output of the Humann2/3
This script will do the following things:
    1. Extract the pathabundance abundance from the output of the Humann2/3.
    2. Pathway abundance at community level is discarded.
    Note: For detailied description of the input file format please refer to the https://github.com/biobakery/humann#2-pathway-abundance-file
    3. TopN pathways in whole experiment was estimated.
    4. Abundance of TopN pathways across the samples is extracted.
    5. Proportion of TopN pathways contributing each sample is estimated.
"""

def extract_pathways():
    input_file=pd.read_table(args.input, header=0) #Reading input
    if '# Pathway' in input_file.columns[0]:
        pa=input_file[~input_file["# Pathway"].str.contains('\|')] #Filtering off the pathways abundance at species level
        pa.columns=pa.columns.str.replace("_Abundance","") #renaming columnames by replacing "_pathabundance"
    
    if 'ID' in input_file.columns[0]:
        pa=input_file[~input_file["ID"].str.contains('\|')] #Filtering off the pathways abundance at species level
        pa.columns=pa.columns.str.replace("_pathabundance","") #renaming columnames by replacing "_pathabundance"
        pd.columns.values[0]="# Pathway" #Renaming
    
    #Discard samples if provided
    if args.samples:
        pa=pa.drop(columns=args.samples) 
    
    pa=pa.loc[pa['# Pathway']!='UNMAPPED'] #Discard the unmapped reads
    pa=pa.loc[pa['# Pathway']!='UNINTEGRATED'] #Discard the unintegrated means genes that do not belongs to any of the known pathwayys
    pa['# Pathway']=pa['# Pathway'].str.split(":",expand=True)[0]
    #pa.to_csv(f'{args.output}_PWabundance.csv',index=False)
    return pa

def topN_pa_across_samples(pa_samples):
    pa_samples.set_index('# Pathway',inplace=True) #Set index at column 'ID'
    pa_samples=pa_samples.assign(total_abundance=lambda x: x.sum(axis=1)) #Total abundance of each pathway in whole experiment 
    pa_samples_sort=pa_samples.sort_values(by=['total_abundance'],ascending=False) #Sort the samples based total abundance in decreasing order
    #pa_samples_sort.to_csv(f'{args.output}-sorted.csv') 
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
    parser.add_argument('-m','--metadata', type=str, help="Metadata about the samples in csv format with SampleId containing the sample names")
    parser.add_argument('-s','--samples', nargs='*', help="List of samples (from columns) to discard from table eg -s S1 S2 S3")
    parser.add_argument('-t','--topN', type=int, help="Top N pathways")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

#Extracting the Pathways abundance (in RPK) and discarding pathways abundance at the community level.
pathways=extract_pathways()

#Estimating topN pathways in the whole experiment and their abundance across the samples.
topNpa_abs,topNpa_rel=topN_pa_across_samples(pathways)

#If metadata file is provided add it to abundance file
if args.metadata:
    metadata=pd.read_csv(args.metadata,header=0,index_col="SampleId")
    topNpa_abs_meta=pd.concat([metadata,topNpa_abs],axis=1,join="inner")
    topNpa_rel_meta=pd.concat([metadata,topNpa_rel],axis=1,join="inner")

    #Write absoulute and relative abundance with metadata
    topNpa_abs_meta.to_csv(f'{args.output}_PWabundance-abs.csv')
    topNpa_rel_meta.to_csv(f'{args.output}_PWabundance-rel.csv')

else:
    #Write absoulute and relative abundance without metadata
    topNpa_abs.to_csv(f'{args.output}_PWabundance-abs.csv')
    topNpa_rel.to_csv(f'{args.output}_PWabundance-rel.csv')
