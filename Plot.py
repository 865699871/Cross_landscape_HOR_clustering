import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpathes
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import argparse


def Plot(sample_state, sequence_len,label6file,label2file,windowfile,outdir,n):
    fig, ax = plt.subplots(figsize=(10, 10))
    max_sample_monomer_len = -1
    for sample in sequence_len.keys():
        if sequence_len[sample] > max_sample_monomer_len:
            max_sample_monomer_len = sequence_len[sample]
    color_list = [ '#646464', '#78C7FF','#0CFDFF' , '#0096FF', '#FFB0D8', '#C70000','#FF2E92']

    color_list_6 = [ '#646464', '#78C7FF','#0CFDFF' , '#0096FF', '#FFB0D8', '#C70000','#FF2E92']

    # color_list_6 = [ '#43978F', '#9EC4BE', '#ABD0F1', '#DCE9F4', '#E56F5E', '#F19685', '#F6C957', '#FFB77F', '#FBE8D5']

    color_list_2 = ['#939393','#88163F']

    custom_lines = []
    legend_text = []
    color_table_6 = {}
    label6file = np.asarray(pd.read_csv(label6file,sep='\t'))

    count = 0
    for i in label6file:
        if i[1] not in color_table_6.keys():
            color_table_6[i[1]] = color_list_6[count]
            count += 1

    color_table_2 = {}
    label2file = np.asarray(pd.read_csv(label2file, sep='\t'))

    count = 0
    for i in label2file:
        if i[1] not in color_table_2.keys():
            color_table_2[i[1]] = color_list_2[count]
            count += 1

    window_color = {0: '#F4F4F4',1: '#FCFCFC', 2: '#F8F8A1', 3: '#F7E779', 4: '#F6CC51', 5: '#F4A727',
                    6: '#EB7908', 7: '#C34C05', 8: '#9B2A03', 9: '#730F01', 10: '#4B0101'}
    windowfile = np.asarray(pd.read_csv(windowfile,sep='\t',header=None))

    sample_of_window_match = {}
    for i in windowfile:
        if i[0] not in sample_of_window_match.keys():
            if int(i[3]) > 10:
                color = '#4B0101'
            else:
                color = window_color[int(i[3])]
            sample_of_window_match[i[0]] = [[int(i[1]),int(i[2]),color]]
        else:
            if int(i[3]) > 10:
                color = '#4B0101'
            else:
                color = window_color[int(i[3])]
            sample_of_window_match[i[0]].append([int(i[1]),int(i[2]),color])


    sample_label6file = {}
    for i in label6file:
        items = i[0].split('_')
        start = int(items[-2].split('-')[0])
        end = int(items[-2].split('-')[1])
        if items[0].startswith('CHM'):
            sample = items[0]
        else:
            sample = items[0] + items[1]

        if sample not in sample_label6file.keys():
            sample_label6file[sample] = [i]
        else:
            sample_label6file[sample].append(i)

    sample_label2file = {}
    for i in label2file:
        items = i[0].split('_')

        start = int(items[-2].split('-')[0])
        end = int(items[-2].split('-')[1])
        if items[0].startswith('CHM'):
            sample = items[0]
        else:
            sample = items[0] + items[1]

        if sample not in sample_label2file.keys():
            sample_label2file[sample] = [i]
        else:
            sample_label2file[sample].append(i)


    pattern_count = 0

    track_gap = 20
    track_number = 2
    main_track_len = track_gap * 2
    sub_track_len = track_gap * 2 * (track_number + 1)

    for sample in sample_label6file.keys():
        sample_monomer_len = sequence_len[sample]
        # main track
        xy = np.array([0, pattern_count * max_sample_monomer_len / track_gap])
        rect = mpathes.Rectangle(xy, sample_monomer_len, max_sample_monomer_len / main_track_len, color='#F4F4F4')
        ax.add_patch(rect)
        # sub tracks
        for t in range(track_number):
            xy = np.array([0, pattern_count * max_sample_monomer_len / track_gap - (t + 1) * max_sample_monomer_len / sub_track_len])
            rect = mpathes.Rectangle(xy, sample_monomer_len, max_sample_monomer_len / sub_track_len, color='#F4F4F4')
            ax.add_patch(rect)
        # main track
        for i in sample_label6file[sample]:
            items = i[0].split('_')
            start = int(items[-2].split('-')[0])
            end = int(items[-2].split('-')[1])
            xy2 = np.asarray([start, pattern_count * max_sample_monomer_len / track_gap])
            rect = mpathes.Rectangle(xy2, end + 1 - start, max_sample_monomer_len / main_track_len, color=color_table_6[i[1]], lw=0)
            ax.add_patch(rect)

        # track number 1
        sub_track_id = 1
        for i in sample_label2file[sample]:
            items = i[0].split('_')
            start = int(items[-2].split('-')[0])
            end = int(items[-2].split('-')[1])
            xy2 = np.array(
                [start, pattern_count * max_sample_monomer_len / track_gap
                 - sub_track_id * max_sample_monomer_len / sub_track_len])
            rect = mpathes.Rectangle(xy2, end + 1 - start, max_sample_monomer_len / sub_track_len, color=color_table_2[i[1]], lw=0)
            ax.add_patch(rect)

        sub_track_id = 2
        if sample in sample_of_window_match.keys():
            for i in sample_of_window_match[sample]:
                start = i[0]
                end = i[1]
                color = i[2]
                xy2 = np.asarray(
                    [start, pattern_count * max_sample_monomer_len / track_gap
                     - sub_track_id * max_sample_monomer_len / sub_track_len])
                rect = mpathes.Rectangle(xy2, end + 1 - start, max_sample_monomer_len / sub_track_len, color=color, lw=0)
                ax.add_patch(rect)

        plt.text(sample_monomer_len + max_sample_monomer_len / main_track_len, pattern_count * max_sample_monomer_len / track_gap, sample + ' ' + sample_state[sample], fontsize=10)
        pattern_count += 1



    xy3 = np.asarray([0, -max_sample_monomer_len / main_track_len])
    rect = mpathes.Rectangle(xy3, max_sample_monomer_len, max_sample_monomer_len / 1000, color='black')
    ax.add_patch(rect)
    point_bar = int(max_sample_monomer_len / 10)
    for i in range(10):
        xy3 = np.asarray([0 + i * point_bar, -max_sample_monomer_len / main_track_len])
        rect = mpathes.Rectangle(xy3, max_sample_monomer_len / 1000, -max_sample_monomer_len / 100, color='black')
        ax.add_patch(rect)
        plt.text(0 + i * point_bar, -max_sample_monomer_len / main_track_len - max_sample_monomer_len / main_track_len, str(0 + i * point_bar), fontsize=5)

    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    # ax.legend(custom_lines,legend_text)
    plt.xticks([])
    plt.yticks([])
    plt.axis('equal')
    plt.savefig(outdir + '/plot_pattern.' + n + '.pdf')
    plt.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input_sample_list", required=True)
    parser.add_argument("-c", "--chr", required=True)
    parser.add_argument("-n", "--cluster_number", required=True)
    parser.add_argument("-o", "--output_path", required=True)
    args = parser.parse_args()

    input_sample_file = args.input_sample_list
    chr = args.chr
    n = args.cluster_number
    outdir = args.output_path

    input_sample = []
    with open(input_sample_file, 'r') as f:
        for line in f:
            if line.startswith('CHM'):
                line = line.strip().split('\t')
                print(line)
                input_sample.append([line[0], line[1]])
            else:
                line = line.strip().split('\t')
                input_sample.append([line[0], line[1], line[2]])

    samples = []
    sequence_len = {}
    sample_state = {}
    for sample in input_sample:
        if sample[0].startswith('CHM'):
            samples.append(sample[0])
            sample_dir = '/data/home/user/home/project/all_human/assembly_part/' + sample[0] + \
                         '/workdir/censeq/chr' + chr
            if os.path.exists(sample_dir + '/input_fasta.1.fa.fai'):
                pass
            else:
                cmd = 'samtools faidx ' + sample_dir + '/input_fasta.1.fa'
                os.system(cmd)
            with open(sample_dir + '/input_fasta.1.fa.fai', 'r') as f:
                line = f.readline()
                length = line.strip().split('\t')[1]
                sequence_len[sample[0]] = int(length)
                sample_state[sample[0]] = sample[1]
        else:
            samples.append(sample[0] + sample[1])
            sample_dir = '/data/home/user/home/project/all_human/assembly_part/' + sample[0] + '/' + \
                         sample[1] + '/workdir/censeq/chr' + chr
            if os.path.exists(sample_dir + '/input_fasta.1.fa.fai'):
                pass
            else:
                cmd = 'samtools faidx ' + sample_dir + '/input_fasta.1.fa'
                os.system(cmd)
            with open(sample_dir + '/input_fasta.1.fa.fai', 'r') as f:
                line = f.readline()
                length = line.strip().split('\t')[1]
                sequence_len[sample[0] + sample[1]] = int(length)
                sample_state[sample[0] + sample[1]] = sample[2]

    Plot(sample_state ,sequence_len, outdir + '/datamatrix.' + n + '.label.xls',
         outdir + '/datamatrix.2.label.xls',
         outdir + '/match_in_window.xls',
         outdir, n)