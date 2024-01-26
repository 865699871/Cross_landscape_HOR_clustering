import pandas as pd
import argparse
from sklearn.cluster import KMeans

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--pairwise_aligned_fasta", required=True)
    parser.add_argument("-n", "--cluster_number", required=True)
    parser.add_argument("-o", "--output_path", required=True)
    args = parser.parse_args()

    file = args.pairwise_aligned_fasta
    output = args.output_path
    n = args.cluster_number

    n = int(n)

    file = open(file, 'r')
    data_table = {}
    key = ''
    # e_count = 0
    # d_count = 0
    # k_count = 0
    seqlen = -1
    for i in file:
        if i.startswith('# 2: '):
            key = i.split('# 2: ')[-1][:-1]
            data_table[key] = {'consensus':'','self':'','marker':''}
        if i.startswith('consensus'):
            seq = i[21:].split(' ')[0]
            seqlen = len(seq)
            data_table[key]['consensus'] += seq
        if i.startswith('CHM'):
            item = i[21:21 + seqlen]
            data_table[key]['self'] += item
        if i.startswith('HG'):
            item = i[21:21 + seqlen]
            data_table[key]['self'] += item
        if i.startswith('NA'):
            item = i[21:21 + seqlen]
            data_table[key]['self'] += item
        if i.startswith('RY'):
            item = i[21:21 + seqlen]
            data_table[key]['self'] += item
        if i.startswith('                     '):
            item = i[21:21 + seqlen]
            data_table[key]['marker'] += item
    file.close()
    datamatrix = {}
    col = []
    for i in data_table.keys():
        datamatrix[i] = []

        for j in range(len(data_table[i]['consensus'])):
            if data_table[i]['consensus'][j] == '-':
                continue
            marker = data_table[i]['marker'][j]
            if marker == '|':
                datamatrix[i].append(0)
            else:
                datamatrix[i].append(1)

        col = []
        for j in range(len(datamatrix[i])):
            col.append('f' + str(j))

    datamatrix = pd.DataFrame(datamatrix).T

    # datamatrix.columns = col
    # datamatrix.to_csv(output, sep='\t')
    # ll


    kmeans =KMeans(n_clusters=n, random_state=10)

    kmeans.fit(datamatrix)
    predict = kmeans.predict(datamatrix)
    print(predict)

    predict = pd.DataFrame(predict,index=datamatrix.index)
    print(predict)
    datamatrix.to_csv(output + 'datamatrix.'+ str(n) + '.xls', sep='\t')
    predict.to_csv(output + 'datamatrix.' + str(n) + '.label.xls', sep='\t')




if __name__ == '__main__':
    main()
