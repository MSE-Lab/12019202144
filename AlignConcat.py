#!/opt/miniconda3/bin/python
# -*- coding: UTF-8 -*-
# author: Zhang Zhen

import glob
import os
import argparse
import sys


def get_parameters():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', required=True, action='store', help='Input alignment files (*.fasta) path')
    parse.add_argument('-o', required=False, action='store', help='Output file name. Default: concatenate_result.fasta',
                       default='concatenate_result.fasta')
    args = parse.parse_args()
    return args


def remove_return(data):
    """
    Remove returns in sequence file
    :param data: a list including sequences ['>seqid', 'ADFTAD...SDTT']
    :return:
    """
    sequence_str = ''
    for i in data:
        if i.startswith('>'):
            sequence_str += '\n' + i
        else:
            sequence_str += i.strip()
    return sequence_str.strip().split('\n')


def read_files(ifp):
    """
    Read files from target path
    :param ifp: input files path
    :return:
    """
    fastas = glob.glob(os.path.join(ifp, '*.fasta'))
    dp = {}
    for fasta in fastas:
        if os.path.isfile(fasta) and os.path.getsize(fasta) > 0:
            with open(fasta, 'r') as file:
                data = remove_return(file.readlines())
                dp[fasta.split('/')[-1]] = data
    return dp


def write_file(file_name, content):
    if os.path.exists(file_name):
        print("WARNing: %s exists, will be over writed" % file_name)
    with open(file_name, 'w') as file:
        file.writelines(content)


def calc_avg_distance(target, others):
    """
    calculate average distance
    :param target: target sequence
    :param others: other sequences
    :return:
    """
    distant_sum = list()
    for genome_id in others:
        diff_site_num = 0
        for site in range(0, len(target)):
            if target[site] != genome_id[site]:
                diff_site_num += 1
        distant_sum.append(diff_site_num / len(target))
    distant_average = sum(distant_sum) / len(others)
    return distant_average


def make_main_dic(dp):
    """
    Data Structure:
    data_pool:{file_name(gene_family): [sequence_id, sequence...]}
    gene_dic_dic:gene_dic_dic:{file_name:{genome_id:{gene_id:sequence}}}
    """
    gdd = {}  # gdd: gene dic dic
    gms = set()  # gms: genome set
    for data in dp:
        gene_dic = {}
        for line in range(int(len(dp[data]) / 2)):
            genome_id = dp[data][2 * line].split('|')[0].strip()
            seq_id = dp[data][2 * line].split('|')[1].strip()
            sequence = dp[data][line * 2 + 1]
            gms.add(dp[data][2 * line].split('|')[0].strip())
            if genome_id not in gene_dic.keys():
                gene_dic[genome_id] = {seq_id: sequence}
            else:
                gene_dic[genome_id][seq_id] = sequence
        gdd[data] = gene_dic
    return gdd, gms


def concatenate(gdd, gms):
    """
    concatenate:three situations
    first: none, we used '-' to fill the sequence
    second: Paralogous, it has Multiple sequences, we calculate average distant,and select the best sequence
    finally: single_sequence, The simplest case that we select directly
    :param gdd: gene_dic_dic
    :param gms: genome set
    """
    result = []
    for genome_key in gms:
        result.append(genome_key)
        sequence_result = ''
        for gene_key in gdd:
            if genome_key not in gdd[gene_key]:
                value_len = len(list(list(gdd[gene_key].values())[0].values())[0])
                sequence_result += '-' * value_len
                print('WARNING: ' + genome_key + ' not found in ' + gene_key)
            elif len(list(gdd[gene_key][genome_key].keys())) > 1:
                distant_dic = {}
                for seq_key_para in gdd[gene_key][genome_key]:
                    gene_family = gdd[gene_key].copy()
                    del gene_family[genome_key]
                    gene_family_list = []
                    for genome_id in gene_family.values():
                        for sequence_id in genome_id:
                            gene_family_list.append(genome_id[sequence_id])
                    sequence = gdd[gene_key][genome_key][seq_key_para]
                    distant = calc_avg_distance(sequence, gene_family_list)
                    distant_dic[seq_key_para] = '%.4f' % distant
                sequence_para_pool = []
                for seq in distant_dic:
                    if distant_dic[seq] == min(list(distant_dic.values())):
                        sequence_para_pool.append(seq)
                        sequence_para_pool.append(gdd[gene_key][genome_key][seq])
                sequence_result += sequence_para_pool[1]
                print('WARNING: Gene duplications ' + ' in ' + gene_key + ': ' + genome_key + ' select ' + sequence_para_pool[0])
            else:
                sequence_result += list(gdd[gene_key][genome_key].values())[0]
        result.append(sequence_result)
    return result


def reformat_sequences(result):
    result_new = ''
    for line in result:
        if '>' in line:
            result_new += line + '\n'
        else:
            for i in range(0, len(line), 80):
                result_new += line[i:i + 80] + '\n'
    return result_new


if __name__ == '__main__':
    arge = get_parameters()
    input_file_path = arge.i
    output_name = arge.o

    data_pool = read_files(input_file_path)
    gene_dic_dic, genome_set = make_main_dic(data_pool)
    result_final = concatenate(gene_dic_dic, genome_set)
    result_return = reformat_sequences(result_final)
    write_file(output_name, result_return)

    print('Gene Numbers: ' + str(len(data_pool)))
    print('Genome Numbers: ' + str(len(genome_set)))
    print('Site Numbers: ' + str(len(result_final[3])))
