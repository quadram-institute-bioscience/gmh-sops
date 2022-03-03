#!/usr/bin/env python3

#from matplotlib import colors
import pandas as pd
import argparse

"""
This code is to extract the genefamilies from the output of the Humann2
"""

def extract_genefamilies():
    input_file=pd.read_table(args.input, header=0)
    gf=input_file[~input_file["ID"].str.contains('\|')] #Filtering off the genefamililes abundance at genus level
    gf.columns=gf.columns.str.replace("_genefamilies","") #renaming columnames by replacing "_genefamilies"
    gf=gf.drop(columns=['GA101o','GN101o','NEG','GP101o']) #droping the blank sample(neg), old sequenced samples:'GA101o','GN101o','GP101o'
    
    #Drop samples
    if args.samples:
        gf=gf.drop(columns=args.samples)

    gf=gf.loc[gf['ID']!='UNMAPPED'] #Single unknown gene of length 1kb consist of all the unmapped reads
    gf.to_csv('genefamilies_mod.csv',index=False)
    return gf

def topN_gf_across_samples(gf_samples):
    gf_samples.set_index('ID',inplace=True)
    gf_samples=gf_samples.assign(total_abundance=lambda x: x.sum(axis=1))
    gf_samples_sort=gf_samples.sort_values(by=['total_abundance'],ascending=False) #Sort the samples based total abundance in decreasing order
    #gf_samples_sort.to_csv(f'{args.output}-sorted.csv')
    gf_samples_sort=gf_samples_sort.drop(columns=['total_abundance']) #Dropping column total_abundance. Since its purpose is served. Not needed in future.
    
    # Select the top N rows from the sorted dataframe
    topNgf_samples_abs=gf_samples_sort.head(args.topN)
    
    #Estimation of others genefamilies abundance within a sample and adding to the topNgf_samples_abs.
    topNgf_samples_abs=topNgf_samples_abs.T #transpose the topN genefamilies dataframe
    #Absoulute abundance of genefamilies in each sample.
    topNgf_samples_abs['Others']=gf_samples_sort.iloc[args.topN::].sum(axis=0)# use the same sorted dataframe to estimate the abundance of others excluding tonp N genefamilies
    
    #Relative abundance of genefamilies within each sample.
    topNgf_samples_rel=topNgf_samples_abs.div(topNgf_samples_abs.sum(axis=1), axis=0) #Estimating the relative abundance of each genefamily within a sample.
    
    return topNgf_samples_abs,topNgf_samples_rel


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__ , formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i','--input',type=str, required=True, help="Gene families summary file across multiple samples")
    parser.add_argument('-o','--output',type=str, help="Output file")
    parser.add_argument('-m','--metadata', type=str, help="metadata about the samples in csv format with SampleId containing the sample names")
    parser.add_argument('-s','--samples', nargs='*', help="list of samples(from columns) to discard from table eg -s GA101o GN101o NEG GP101o")
    parser.add_argument('-c','--cellpopn', nargs='*', help="list of samples containing (S and/or N and/or P) to discard from table eg -s A N")
    parser.add_argument('-t','--topN', type=int, help="top N genefamilies")
    args = parser.parse_args()

print(list(args.samples))
genefamilies=extract_genefamilies()
topNgf_abs,topNgf_rel=topN_gf_across_samples(genefamilies)

if args.metadata:
    metadata=pd.read_csv(args.metadata,header=0,index_col="SampleId")
    topNgf_abs_meta=pd.concat([metadata,topNgf_abs],axis=1,join="inner")
    topNgf_rel_meta=pd.concat([metadata,topNgf_rel],axis=1,join="inner")
    topNgf_abs_meta.to_csv(f'{args.output}-abs.csv')
    topNgf_rel_meta.to_csv(f'{args.output}-rel.csv')

else:
    #Write absoulute and relative abundance.
    topNgf_abs.to_csv(f'{args.output}-abs.csv')
    topNgf_rel.to_csv(f'{args.output}-rel.csv')


    
