import os

def reverse(kmer):
    base_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'a': 'T', 't': 'A', 'c': 'G', 'g': 'C', 'n': 'N'}
    sequence_list = []
    for i in kmer[::-1]:
        sequence_list.append(base_map[i])
    new_sequence = ''.join(sequence_list)
    return new_sequence

def main():
    pattern_file = '/data/home/user/home/project/all_human/workdir/multi_pattern_tree/target_pattern.txt'
    calculate_consensus = '/data/home/user/home/project/workdir/script/consensus_sequence.py'
    pattern_list = []
    with open(pattern_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            seq = line[2].split('_')
            pattern_list.append([line[0], line[1], seq])

    for pattern in pattern_list:
        cmd = 'mkdir -p /data/home/user/home/project/all_human/workdir/multi_pattern_tree/chr' \
              + pattern[0] + '_' + pattern[1]
        os.system(cmd)
        cmd = 'mkdir -p /data/home/user/home/project/all_human/workdir/multi_pattern_tree/chr' \
              + pattern[0] + '_' + pattern[1] + '/samples'
        os.system(cmd)

        pattern_length = len(pattern[2])
        outdir = '/data/home/user/home/project/all_human/workdir/multi_pattern_tree/chr' \
                 + pattern[0] + '_' + pattern[1]

        input_sample_file = '/data/home/user/home/project/all_human/workdir/pattern_cluster/samples/' \
                            'sample_list_chr' + pattern[0] + '.txt'
        input_sample = []
        with open(input_sample_file, 'r') as f:
            for line in f:
                if line.startswith('CHM'):
                    input_sample.append([line.strip().split('\t')[0]])
                else:
                    line = line.strip().split('\t')
                    input_sample.append([line[0], line[1]])

        reverse_pattern = list(reversed(pattern[2]))
        for sample in input_sample:
            if sample[0].startswith('CHM'):
                workdir = '/data/home/user/home/project/all_human/assembly_part/' + sample[0] + \
                          '/workdir/censeq/chr' + pattern[0]
                monomer_file = workdir + '/out_monomer_seq.xls'
                position_file = workdir + '/final_decomposition.tsv'
                fasta_file = workdir + '/input_fasta.1.fa'
                sample_name = sample[0] + '_H'
            else:
                workdir = '/data/home/user/home/project/all_human/assembly_part/' + sample[0] + '/' + \
                          sample[1] + '/workdir/censeq/chr' + pattern[0]
                monomer_file = workdir + '/out_monomer_seq.xls'
                position_file = workdir + '/final_decomposition.tsv'
                fasta_file = workdir + '/input_fasta.1.fa'
                sample_name = sample[0] + '_' + sample[1]

            with open(monomer_file, 'r') as f:
                line = f.readline()
                mon_seq = line.strip().split('mon_seq\t')[1].split(' ')

            target_pattern = {}
            for i in range(len(mon_seq) - pattern_length + 1):
                curr_seq = mon_seq[i: i + pattern_length]
                if curr_seq == pattern[2] or curr_seq == reverse_pattern:
                    target_pattern[i] = ''

            monomer_position = []
            with open(position_file, 'r') as f:
                for line in f:
                    line = line.strip().split('\t')
                    left = line[2]
                    right = line[3]
                    if line[1].endswith('\''):
                        status = '-'
                    else:
                        status = '+'
                    monomer_position.append([left, right, status])

            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        continue
                    else:
                        cen_fa = line.strip()

            for i in target_pattern.keys():
                start = int(monomer_position[int(i)][0])
                end = int(monomer_position[int(i)+pattern_length-1][1])
                count = 0
                for j in range(int(i), int(i)+pattern_length):
                    if monomer_position[j][2] == '+':
                        count += 1
                if count == 0:
                    final_status = '-'
                elif count == pattern_length:
                    final_status = '+'
                else:
                    print('error')
                    print(left)
                    print(right)
                    return
                if final_status == '+':
                    target_pattern[i] = [start, end, cen_fa[start:(end+1)]]
                else:
                    target_pattern[i] = [start, end, reverse(cen_fa[start:(end+1)])]

            with open(outdir + '/samples/' + sample_name + '_' + pattern[0] + '_' + pattern[1] + '.fa', 'w') as f:
                if sample[0].startswith('CHM'):
                    for info in target_pattern.values():
                        f.writelines('>' + sample[0] + '_H_' + pattern[1] + '_' + str(info[0]) + '-' + str(info[1])
                                     + '_' + '\n')
                        f.writelines(str(info[2]) + '\n')
                else:
                    for info in target_pattern.values():
                        f.writelines('>' + sample[0] + '_' + sample[1] + '_' + pattern[1] + '_' + str(info[0]) + '-'
                                     + str(info[1]) + '_' + '\n')
                        f.writelines(str(info[2]) + '\n')

        cmd = 'cat ' + outdir + '/samples/*_' + pattern[0] + '_' + pattern[1] + '.fa >' + outdir + '/' + pattern[1] + \
              '.fa'
        print(cmd)
        os.system(cmd)
        cmd = 'kalign -i ' + outdir + '/' + pattern[1] + '.fa -o ' + outdir + '/' + pattern[1] + '.align.fa'
        print(cmd)
        os.system(cmd)
        cmd = 'python ' + calculate_consensus + ' -i ' + outdir + '/' + pattern[1] + '.align.fa -o ' + outdir + \
              '/' + pattern[1] + '.consensus.fa'
        print(cmd)
        os.system(cmd)
        cmd = 'needle -asequence ' + outdir + '/' + pattern[1] + '.consensus.fa -bsequence ' + outdir + '/' + \
              pattern[1] + '.fa -gapopen 10.0 -gapextend 0.5 -outfile ' + outdir + '/' + pattern[1] + \
              '.pairwise.align.fa'
        print(cmd)
        os.system(cmd)

        for n in range(2, 10):
            process_needle = '/data/home/user/home/project/workdir/script/processNeedle.py'
            recent_expand = '/data/home/user/home/project/workdir/script/recentExpand.py'
            plot = '/data/home/user/home/project/workdir/script/Plot.py'

            cmd = 'python ' + process_needle + ' -i ' + outdir + '/' + pattern[1] + '.pairwise.align.fa -n ' + \
                  str(n) + ' -o ' + outdir + '/'
            print(cmd)
            os.system(cmd)

            cmd = 'python ' + recent_expand + ' -i ' + outdir + ' -f ' + '/' + pattern[1] + \
                  '.fa -n ' + str(n)
            print(cmd)
            os.system(cmd)

            cmd = 'python ' + plot + ' -i ' + input_sample_file + ' -c ' + pattern[0] + ' -n ' + str(n) + ' -o ' + \
                  outdir
            print(cmd)
            os.system(cmd)


if __name__ == '__main__':
    main()