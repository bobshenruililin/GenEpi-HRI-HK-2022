#python extract_fasta.py multi.fasta extract.list extract.fasta

import sys

def extract_fasta(fasta_file, txt_file, output_file):
    fasta_dict = {}
    current_seq = ""
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if current_seq:
                    fasta_dict[current_seq] = sequence
                current_seq = line.strip()[1:]
                sequence = ""
            else:
                sequence += line.strip()
        if current_seq:
            fasta_dict[current_seq] = sequence

    extracted_seqs = []
    with open(txt_file, 'r') as f:
        for line in f:
            seq_name = line.strip()
            if seq_name in fasta_dict:
                extracted_seqs.append(">" + seq_name)
                extracted_seqs.append(fasta_dict[seq_name])

    with open(output_file, 'w') as f:
        f.write("\n".join(extracted_seqs))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py fasta_file txt_file output_file")
        sys.exit(1)

    fasta_file = sys.argv[1]
    txt_file = sys.argv[2]
    output_file = sys.argv[3]

    extract_fasta(fasta_file, txt_file, output_file)
