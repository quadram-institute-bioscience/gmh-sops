#!/usr/bin/env python
"""
Extract the reads count from FASTP (json) and print sample name and counts
"""

import sys, os
import json

if __name__ == "__main__":
    import argparse
    args = argparse.ArgumentParser(description="Extract the reads count from FASTP (json) and print a Seqfu like output")
    args.add_argument("JSONFILES", nargs="+", help="Input FASTP json file")
    args.add_argument("-o", "--output", help="Output file")
    args.add_argument("--header", help="Add header line", action="store_true")
    args.add_argument("--fields", help="Number of fields to print", type=int, default=2)
    args = args.parse_args()

    if args.output:
        outfile = open(args.output, "w")
    else:
        outfile = sys.stdout
    
    # SeqFu stats header (not printed by default)
    fields = ["File","#Seq","Total bp","Avg", "N50", "N75", "N90", "auN", "Min", "Max"]
    if args.header:
        print("\t".join(fields[:args.fields]), file=outfile)

    # Create a list "values" with as many elements as fields, all set to ""
    values = [ "-1" for i in range(len(fields)) ]

    for jsonfile in args.JSONFILES:
        values[0] = os.path.basename(jsonfile).replace(".fastp","").replace(".json","")
        data = []
        try:
            with open(jsonfile) as f:
                data = json.load(f)
        except Exception as e:
            print("Error reading %s: %s" % (jsonfile, str(e)), file=sys.stderr)
            continue
        tot_reads = 0
        try:
            tot_reads = data["summary"]["after_filtering"]["total_reads"]
            tot_bases = data["summary"]["after_filtering"]["total_bases"]
            avg_len =   data["summary"]["after_filtering"]["read1_mean_length"]
        except Exception as e:
            print("Error: cannot find the required fields in %s: %s" % jsonfile, e, file=sys.stderr)

        values[1] = str(tot_reads)
        values[2] = str(tot_bases)
        values[3] = str(avg_len)

        # Print fields tab separated
        print("\t".join(values[:args.fields]), file=outfile)
