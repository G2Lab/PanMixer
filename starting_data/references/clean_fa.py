def remove_underscore_sequences(input_fasta, output_fasta):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        keep = False
        for line in infile:
            if line.startswith('>'):
                # Only keep headers with no underscores
                if '_' not in line:
                    keep = True
                    outfile.write(line)
                else:
                    keep = False
            else:
                if keep:
                    outfile.write(line)

# Example usage:
remove_underscore_sequences('hg38.fa', 'hg38_cleaned.fa')
