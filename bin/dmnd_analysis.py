#!/usr/bin/env python

from Bio import SeqIO
import collections
import argparse
import pprint
import csv
import sys
import os



class DMNDanalysis(object):

    def __init__(self, infile, reference, outdir, prefix):

          self.dmd_fields = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
          self.infile = infile
          self.fasta_ref = reference
          self.outdir =  outdir
          self.prefix = prefix
          
          
    def  parse_dmnd(self):

         self.refseq_hits = collections.defaultdict(list)
         with open(self.infile)  as csv_fobj:
             csv_reader = csv.DictReader(csv_fobj, self.dmd_fields, delimiter='\t')
             for row in csv_reader:
                 self.refseq_hits[row["sseqid"]].append(row["qseqid"])       

         self.hit_counts = dict([(hit, len(self.refseq_hits[hit]) ) for hit in self.refseq_hits ])
         self.sorted_hits =  sorted(self.hit_counts, key=self.hit_counts.get, reverse=True)


    def parse_ref(self):

        self.ref_data = collections.defaultdict(tuple)
        ref_seqs  =  SeqIO.parse(self.fasta_ref,  "fasta")
        for seq in ref_seqs:
            try:
                desc, organism  =  seq.description.split(" ",1)[-1].split("[", 1)
            except ValueError:
                desc  =  seq.description.split(" ",1)[-1]
                organism =  ''

            self.ref_data[seq.id] = (desc, organism.rstrip("]") )
            

    def get_hits(self):

        self.function = collections.defaultdict(int)
        self.taxonomy = collections.defaultdict(int)
        total = sum(self.hit_counts.values())
        
        fname = os.path.join(self.outdir, self.prefix+'_hits.tsv')
        csv_fobj = csv.writer(open(fname , "w"), delimiter = "\t")
        csv_err_fobj = csv.writer(open(fname+'.err', 'w'), delimiter ="\t" )
        total_fname = os.path.join(self.outdir, self.prefix+'_totals.tsv')
        self.obj_totals = csv.writer(open(total_fname, "w"), delimiter = "\t")
     
        for hit in self.sorted_hits:
            count =  self.hit_counts[hit]
            perc = round(100* count/float(total), 4)
            annotation  = self.ref_data[hit]
            if annotation:
                desc, organism = annotation
                self.function[desc] += count
                self.taxonomy[organism] += count
                row  = [count, perc, hit, organism,desc]
                csv_fobj.writerow(row)
            else:
                csv_err_fobj.writerow([count, perc, hit, None,None])

                
    def analyse(self):
        
        fname1 = os.path.join(self.outdir, self.prefix+'_function.tsv')
        csv1_fobj = csv.writer(open(fname1 , "w"), delimiter = "\t")
        sum_function = sum(self.function.values())
        self.obj_totals.writerow([fname1, sum_function])
        for  func in sorted(self.function, key=self.function.get, reverse=True):
            func_count = self.function[func]
            func_perc =  round(100 * func_count/float(sum_function),4)
            func_row = [func_count, func_perc, func]
            csv1_fobj.writerow(func_row)
     

        
        fname2 = os.path.join(self.outdir, self.prefix+'_organism.tsv')
        csv2_fobj = csv.writer(open(fname2 , "w"), delimiter = "\t")
        sum_taxonomy = sum(self.taxonomy.values())
        self.obj_totals.writerow([fname2, sum_taxonomy])
        for organism in sorted(self.taxonomy,  key=self.taxonomy.get, reverse=True):
            org_count = self.taxonomy[organism]
            org_perc =  round(100 * org_count/float(sum_taxonomy), 4)
            org_row = [org_count, org_perc, organism]
            csv2_fobj.writerow(org_row)



                                       
if __name__ == '__main__':
    parser = argparse.ArgumentParser("Parses a DIAMOND annotation result file (in BLAST m8 format) to generate counts data")
    parser.add_argument('infile', help = "specifies the infile (a DIAMOND results file in m8 format)")
    parser.add_argument('-r','--reference', help = "specifies a reference fasta to search against for results", required=True )
    parser.add_argument('-o','--outdir', help = "output directory. Default is working directory")
    parser.add_argument('-p','--prefix', default="SHB_", help = "output directory. Default is working directory")
    args = parser.parse_args()
    outdir = args.outdir if args.outdir else os.getcwd()
    sys.stderr.write(outdir+"\n")
    dmnd = DMNDanalysis(args.infile, args.reference, outdir, args.prefix)
    dmnd.parse_dmnd()
    dmnd.parse_ref()
    dmnd.get_hits()
    dmnd.analyse()
    
