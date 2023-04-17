from Bio import SeqIO
import os

INPUTDIR = "."
list_of_files=[ fasta for fasta in os.listdir(INPUTDIR) if fasta.endswith(".ffn") ]

def filter_fasta(input,output):
    records = SeqIO.parse(input, 'fasta')
    filtered = (rec for rec in records if 'G-binding' in rec.description)
    with open (output, "w") as out:
        SeqIO.write(filtered, out, 'fasta')

if __name__ == "__main__":
    for input in list_of_files:
        output="".join(str(input).split(".")[:-1])+"_filtered.fasta"
        filter_fasta(input,output)

