from collections import Counter
import argparse

def read_fasta_file(file):
    genome = []
    with open(file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                genome.append('')
            else:
                genome[-1] += line
    return genome

def consensus_sequence(genome):
    result = []
    length = len(genome[0])
    for i in range(length):
        chars = [g[i] for g in genome]
        counter = Counter(chars)
        if '-' in counter:
            del counter['-']
        most_common = counter.most_common(1)[0]
        result.append(most_common[0])
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--aligned_fasta", required=True)
    parser.add_argument("-c", "--consensus_name", required=False,default='consensus')
    parser.add_argument("-o", "--output_file", required=True)
    args = parser.parse_args()
    fasta_file = args.aligned_fasta
    name = args.consensus_name
    output_file = args.output_file
    genome = read_fasta_file(fasta_file)
    result = consensus_sequence(genome)
    out = ''.join(result)
    with open(output_file, 'w') as o:
        o.writelines('>' + name + '\n')
        o.writelines(out)
        o.writelines('\n')