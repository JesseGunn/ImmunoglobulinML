from __future__ import annotations
import pandas as pd
import twobitreader
import itertools
import time
import warnings

from Bio.Seq import Seq

NUCLEOTIDES = ['A', 'C', 'G', 'T']

def make_kmer_col(orig_df: pd.DataFrame, col_prefix: str, k=3, norm=False) -> pd.DataFrame:
    # Generate all possible k-mer combinations
    combos = [''.join(p) for p in itertools.product(NUCLEOTIDES, repeat=k)]

    def count_kmer(seq, kmer):
        return sum(seq[i:i+len(kmer)] == kmer for i in range(len(seq) - len(kmer) + 1))
    
    kmer_counts = {
    f"{col_prefix}_{combo}": orig_df['seq'].apply(
        lambda seq: count_kmer(seq, combo) / len(seq) if norm else count_kmer(seq, combo)
    )
    for combo in combos
}
    
    return pd.concat([orig_df, pd.DataFrame(kmer_counts, index=orig_df.index)], axis=1)

def create_kmer_df(df, output_file, kmer_types = [2, 3, 6]) :
    # Generate k-mer columns for each specified k-mer type
    kmer_df = df.copy()
    for kmer in kmer_types:
        print(f"Processing {kmer}-mers for {output_file}_kmers...")
        start_time = time.time()
        kmer_df = make_kmer_col(kmer_df, f"{kmer}mer", k= kmer, norm=True)
        elapsed_time = time.time() - start_time
        print(f"Time taken for {kmer}-mers: {elapsed_time:.2f} seconds")

    # Write to CSV
    kmer_df.to_csv(output_file + f"_{str(kmer_types).replace('[', '').replace(']', '').replace(',', '_')}" + "_kmers.csv", index=False)
    return kmer_df    

def reduce_columns(df, output_file):
    df["source"] = output_file
    df['length'] = df['seq'].str.len()  # Add sequence length as a feature
    # Using str.extract with regex
    df[['gene_type', 'gene_function']] = df['feature'].str.extract(r'(\S+)\s+\(([^)]+)\)')
    return df[['contig', 'source', 'feature', 'start', 'end', 'length', 'seq']]  # Placeholder - adjust based on actual needs

def get_seqs(twobit, df, output_name):
    # Create a function to fetch sequences based on the coordinates and strand information
    def fetch_seq(row):
        seq = tb[row['contig']][row['start']-1:row['end']]
        if row['strand'] == '-': # Handle reverse complement for negative strand
            return str(Seq(seq).reverse_complement())
        return seq
    
    with twobitreader.TwoBitFile(twobit) as tb: # Do this so we close the file after we're done
        df['seq'] = df.apply(fetch_seq, axis=1)
    
    df = reduce_columns(df, output_name) # Reduce columns to only those needed for k-mer analysis
    return df

def beautify_gff(twobit, gff_file, output_name):
    gff_df = pd.read_csv(gff_file, sep='\t', header=None,
                     names=['contig', 'source', 'feature', 'start', 'end',
                            'score', 'strand', 'frame', 'attributes'])
    
    gff_df = get_seqs(twobit, gff_df, output_name)
    gff_df = create_kmer_df(gff_df, output_name)

    return gff_df