#!/usr/bin/env python3
#for: humann.nf
import pandas as pd
import argparse

"""
This script will do the following things:
    1. Extract the genefamilies abundance from the output of the Humann2/3.
    2. Genefamilies abundance at community level is discarded.
    Note: For detailied description of the input file format please refer to the https://github.com/biobakery/humann#1-gene-families-file
    3. TopN Gene families in whole experiment was estimated.
    4. Abundance of TopN Gene families across the samples is extracted.
    5. Proportion of TopN Gene families contributing each sample is estimated.
"""

def extract_genefamilies():
    input_file=pd.read_table(args.input, header=0) #Reading input
    
    if '# Gene Family' in input_file.columns[0]:
        gf=input_file[~input_file["# Gene Family"].str.contains('\|')] #Filtering off the genefamililes abundance at genus level
        gf.columns=gf.columns.str.replace("_Abundance-RPKs","") #renaming columnames by replacing "_genefamilies"
        

    if 'ID' in input_file.columns[0]:
        gf=input_file[~input_file["ID"].str.contains('\|')] #Filtering off the genefamililes abundance at genus level
        gf.columns=gf.columns.str.replace("_genefamilies","") #Renaming the columns
        gf.columns.values[0] = "# Gene Family"

    #Discarding the samples if provided.
    if args.samples:
        gf=gf.drop(columns=args.samples)
        
    unmapped=gf.loc[gf['# Gene Family']=='UNMAPPED'] #Not mapped reads to any of the genefamilies
    gf=gf.loc[gf['# Gene Family']!='UNMAPPED'] #Single unknown gene of length 1kb consist of all the unmapped reads
    #Dicarding the UNMAPPED reads.
    #gf.to_csv("GeneFamilies_abundance-RPKs.csv",index=False) #File consist of the Gene families abundance (in RPK unit) after the above filteration steps.
    return gf, unmapped

def qc_data(unmap):
    stats=pd.read_table(args.stats,header=0,index_col="File")
    unmap.set_index('# Gene Family',inplace=True)
    unmap=unmap.T
    qc_data=pd.concat([stats,unmap],axis=1,join="inner")
    qc_data['Mapped']=qc_data['#Seq']-qc_data['UNMAPPED']
    qc_data=qc_data.rename(columns={"UNMAPPED":"Unmapped","#Seq":"total_reads"})
    qc_data.to_csv(f'{args.qc}_report.csv')
    return

def topN_gf_across_samples(gf_samples):
    gf_samples.set_index('# Gene Family',inplace=True) #Setting index at column name: ID.
    gf_samples=gf_samples.assign(total_abundance=lambda x: x.sum(axis=1)) #Adding up the abudance of each gene family across the samples.
    gf_samples_sort=gf_samples.sort_values(by=['total_abundance'],ascending=False) #Sort the samples based total abundance in decreasing order
    #gf_samples_sort.to_csv(f'{args.output}-sorted.csv')
    gf_samples_sort=gf_samples_sort.drop(columns=['total_abundance']) #Dropping column total_abundance. Since its purpose is served. Not needed in future.
    
    # Select the top N rows from the sorted dataframe
    topNgf_samples_abs=gf_samples_sort.head(args.topN)
    
    #Estimation of others genefamily abundance within a sample and adding to the topNgf_samples_abs.
    topNgf_samples_abs=topNgf_samples_abs.T #transpose the topN genefamilies dataframe
    #Absoulute abundance of top N genefamilies including 'others' in each sample.
    topNgf_samples_abs['Others']=gf_samples_sort.iloc[args.topN::].sum(axis=0) # use the same sorted dataframe to estimate the abundance of others excluding tonp N genefamilies
    
    #Relative abundance of genefamilies within each sample.
    topNgf_samples_rel=topNgf_samples_abs.div(topNgf_samples_abs.sum(axis=1), axis=0) #Estimating the relative abundance of each genefamily within a sample.
    
    return topNgf_samples_abs,topNgf_samples_rel


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__ , formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i','--input', type=str, required=True, help="Gene families summary file across multiple samples")
    parser.add_argument('-o','--output', type=str, required=True, help="Output file")
    parser.add_argument('-st','--stats', type=str, help="Ouptput of the seqfu stats")
    parser.add_argument('-qc','--qc',type=str, help="Save the qc report")
    parser.add_argument('-m','--metadata', type=str, help="Metadata about the samples in csv format with colum name: 'SampleId' containing the sample names")
    parser.add_argument('-s','--samples', nargs='*', help="list of samples(from columns) to discard from table eg -s S1 S2 S3")
    parser.add_argument('-t','--topN', type=int, help="top N genefamilies")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()


#Extracting the genefamilies abundance (in RPK) and discarding Gene family abundance at the community level.
genefamilies,unmapped_reads=extract_genefamilies()

qc_data(unmapped_reads)
#Estimating topN Gene families in the whole experiment and their abundance across the samples.
topNgf_abs,topNgf_rel=topN_gf_across_samples(genefamilies)

#Adding metadata of the samples if provided.
if args.metadata:
    metadata=pd.read_csv(args.metadata,header=0,index_col="SampleId")
    topNgf_abs_meta=pd.concat([metadata,topNgf_abs],axis=1,join="inner")
    topNgf_rel_meta=pd.concat([metadata,topNgf_rel],axis=1,join="inner")
    topNgf_abs_meta.to_csv(f'{args.output}_GFabundance-abs.csv')
    topNgf_rel_meta.to_csv(f'{args.output}_GFabundance-rel.csv')

else:
    #Write absoulute and relative abundance.
    topNgf_abs.to_csv(f'{args.output}_GFabundance-abs.csv')
    topNgf_rel.to_csv(f'{args.output}_GFabundance-rel.csv')