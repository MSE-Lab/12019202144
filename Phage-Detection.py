#!/usr/bin/python3
#author zhen zhang
import argparse
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import time
import numpy as np
import pymongo


def get_parameters():
    parse = argparse.ArgumentParser(description='This is a script for annotating genome sequence and searching sequences with virus database')
    parse.add_argument('-q', required=True, action='store', help='Query sequences diretory')
    parse.add_argument('--type',required=True, action='store', help='query sequences type: prot for protein, nucle for nucleotide')
    parse.add_argument('--strand', required=False, action='store_true', help='Genes arrangement ignore the strands default: True',default=False)
    parse.add_argument('-r', required=False, action='store', help='Minimum genes in each gene cluster',default='3')
    parse.add_argument('-o', required=True, action='store', help='Output directory')
    parse.add_argument('-t', required=False, action='store', help='Thread numbers default: 8', default='8')
    parse.add_argument('-e', required=False, action='store', help='blastp E-value',default='1e-5')
    parse.add_argument('-c', required=False, action='store', help='Sequences Coverage',default='0.4')
    parse.add_argument('-d', required=False, action='store', help='Database name',default='Virusdb')
    parse.add_argument('-s', required=False, action='store', help='Similarty',default='50')
    parse.add_argument('--clear', required=False, action='store_true', help='Restore Cash file default: True',default=False)
    arguments = parse.parse_args()
    return arguments


def inputfile(file_name):
    with open(file_name) as f:
        data = f.readlines()
    return data


def output_file(output_fileName, content):
    with open(output_fileName, 'w') as file:
        file.writelines(content)


def remove_return(seq_file):
    """Remove sequnences return"""
    seq_not_return = ''
    for line in seq_file:
        if line.startswith('>'):
            seq_not_return += '\n' + line.strip().split(' ')[0] + '\n'
        else:
            seq_not_return += line.strip()
    return seq_not_return


def TimeUsedComputation(StartTime, EndTime):
    timeUsed = EndTime - StartTime
    time_str = "time used: {:.0f}h {:.0f}m {:.0f}s"
    TimeUsed = time_str.format(timeUsed // 3600, (timeUsed % 3600) // 60, ((timeUsed % 3600) % 60) % 60)
    return TimeUsed


def seq_id_pair(recode_sequences):
    """Make sequences id pair with sequences"""
    seq_id_dic = {}
    seq_list = recode_sequences.strip().split('\n')
    for index in range(0, int(len(seq_list) / 2)):
        seq_id_dic[seq_list[2 * index][1:]] = seq_list[2 * index + 1]
    return seq_id_dic


def recode_input_locustag(genome_name,sequences):
    """recode id of genes from query genome """
    locustag_pair = {}
    recode_sequences = ''
    locus_sort = 0
    for line in sequences:
        if line.startswith('>'):
            contig = line.strip('>').split(' ')[0].split('|')[-1].split('_')[0]
            site_location = ''
            for iterms in line.strip('\n').split(' '):
                if 'location' in iterms:
                    site_location = iterms[1:-1].split('=')[-1]
            else:
                pass
            locus = '>%s_genes_%d|%s_prot|%s'%(genome_name,locus_sort,contig,site_location)
            locustag_pair[locus[1:]] = line
            recode_sequences +='\n' + locus + '\n'
            locus_sort += 1
        else:
            recode_sequences += line.strip()
    return recode_sequences,locustag_pair


def generate_queries_file(seq_pair,outdir):
    """Generate query genes from input genome sequences"""
    query_paths = []
    for locus,seq in seq_pair.items():
        new_sequences = '>%s'%locus+'\n'+seq
        path1 = os.path.join(outdir, '%s.fasta'%locus.split('|')[0])
        output_file(path1, new_sequences)
        query_paths.append(path1)
    return query_paths


def prepare_genome_file(input_file, working_dir):
    """Read genome files and recode genome files  """
    query_genome_path = os.path.join(working_dir,'query_genomes')
    try:
        os.mkdir(query_genome_path)
    except FileExistsError:
        pass
    genome_locus_pair = {} # {"recode genome id":"orignal genes locutags"}
    recode_genome_path = [] # ["recode genes path used for blast"]
    genome_all_seq = ''
    genome_num = 0
    genome_file_pair = {} #{"recode genome name":"orignal genome name "} for output result
    genome_name = 'genome%d'%genome_num
    if input_file.strip().split('.')[-1] in ['fa','faa','fna','fasta','fas']:
        genome_file = inputfile(input_file)
        recode_sequences,locustag_pair = recode_input_locustag(genome_name,genome_file)
        for k,v in locustag_pair.items():
            genome_locus_pair[k]=v
        genome_all_seq +=recode_sequences
        seq_pair = seq_id_pair(recode_sequences)
        query_paths = generate_queries_file(seq_pair,query_genome_path)
        recode_genome_path.append(query_paths)
        if '/' in input_file:
            genome_file_pair[genome_name] = input_file.split('/')[-1]
        else:
            genome_file_pair[genome_name] = input_file
    genome_sequences_pair = seq_id_pair(genome_all_seq) #{"recode genes locustags": sequences of genes from query genome"}
    return genome_locus_pair,genome_sequences_pair,recode_genome_path,genome_file_pair


def annotation_nucle_genome(fna_genome_file,working_dir,anotation_thread):
    genome_file_pair = {}
    annotation_file_dic ={}
    genome_num = 0
    if fna_genome_file.strip().split('.')[-1] in ['fa','faa','fna','fasta','fas']:
        genome_name = 'genome%d' % genome_num
        outdir = os.path.join(working_dir, genome_name)
        cmd = ' '.join(["prokka", "--outdir", outdir,"--quiet","--force","--noanno","--norrna","--notrna", "--prefix",genome_name,"--locustag","'%s_genes'"%genome_name,"--cpus",anotation_thread,fna_genome_file])
        pro = subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        pro.wait()
        if 'prokka: not found' in pro.stderr.read():
            print('prokka not found !!! please check you conda environments')
            sys.exit(0)
        annotation_file_dic[genome_name] = (os.path.join(outdir,'%s.faa'%genome_name),os.path.join(outdir,'%s.gff' % genome_name))
        if '/' in input_file:
            genome_file_pair[genome_name] = fna_genome_file.split('/')[-1]
        else:
            genome_file_pair[genome_name] = fna_genome_file
    return genome_file_pair,annotation_file_dic


def read_annotation_file(annotation_file_dic):
    genome_locus_pair ={}
    genome_sequences_pair ={}
    recode_genome_path = []
    genome_seq_path = os.path.join(working_dir,'query_genes')
    try:
        os.mkdir(genome_seq_path)
    except FileExistsError:
        pass
    for genome,anotation_file in annotation_file_dic.items():
        faa_file = anotation_file[0]
        faa_file_seq = inputfile(faa_file)
        gff_file = inputfile(anotation_file[1])
        seq_not_return = remove_return(faa_file_seq)
        seq_pair = seq_id_pair(seq_not_return)
        seq_pair_recode = {}
        for seq_id,seq in seq_pair.items():
            genome_sequences_pair[seq_id] = seq
            num = 0
            for line in gff_file:
                if seq_id in line:
                    genome = seq_id.strip('>').split('_')[0]
                    if line.split('\t')[6] == '+':
                        recode_seq = '>%s_genes_%d|%s_prot|[location=%s..%s]'%(genome,num,line.split('\t')[0],line.split('\t')[3],line.split('\t')[4])
                        genome_locus_pair[recode_seq[1:]] = seq_id
                        genome_sequences_pair[recode_seq[1:]] = seq
                        seq_pair_recode[recode_seq[1:]] = seq
                    else:
                        recode_seq ='>%s_genes_%d|%s_prot|[location=complement(%s..%s)]'%(genome,num,line.split('\t')[0],line.split('\t')[3],line.split('\t')[4])
                        genome_locus_pair[recode_seq[1:]] = seq_id
                        genome_sequences_pair[recode_seq[1:]] = seq
                        seq_pair_recode[recode_seq[1:]] = seq
                else:
                    pass
                num+=1
        query_path = generate_queries_file(seq_pair_recode,genome_seq_path)
        recode_genome_path.append(query_path)
    return genome_locus_pair,genome_sequences_pair,recode_genome_path


def make_mongo_db_dic(refer_DB):
    """Search database and retrive reference sequnences from database"""
    sequences_inf = {}
    seq_refers = ''
    refer_pair = {}
    path1 = os.path.join(refer_DB,'refer_gene_db')
    try:
        os.mkdir(path1)
    except FileExistsError:
        pass
    client = pymongo.MongoClient('mongodb://localhost:27017/')
    mydb = client['virousDB']
    mycol = mydb['virousDB']
    myquery = {'Rank':'Species','Host Source':{'$in':['bacteria']}}
    search_result = mycol.find(myquery)
    for document in search_result:
        # print(document[1]['Host Source'])
        # taxa_dic[document['Name']] = {}
        for ncle_record in document['Nucleotide Data']:
            # taxa_dic[document['Name']][ncle_record['isolate name']] = {}
            for genome in ncle_record['genomic']:
                for cds in genome['features']:
                    locus = '>%s'%cds['protein id']
                    sequences = cds['translation']
                    seq = '%s\n%s'%(locus,sequences)
                    seq_refers += seq+'\n'
                    refer_pair[locus] = sequences
                    sequences_inf[cds['protein id']] = {'speceis name':document['Name'],'isolate name':ncle_record['isolate name'],'accessions':genome['accessions'],'product':cds['product']}
    path2 = os.path.join(path1,'refer.fasta')
    output_file(path2,seq_refers)
    species_inf = {inf['accessions']: [inf['speceis name'],inf['isolate name']] for pid, inf in sequences_inf.items()}
    return sequences_inf,path2,refer_pair, species_inf


def make_blast_db(references_db_path):
    db = '%sdb' % references_db_path
    cmd = ' '.join(['makeblastdb', '-dbtype', 'prot', '-in', references_db_path, '-out', db, '>', '/dev/null'])
    pro = subprocess.Popen(cmd,shell=True)
    pro.wait()
    return db


def blast_search(query,db,evalue,blast_path):
    blast_out = os.path.join(blast_path,'%s.out'%(query.split('/')[-1]))
    cmd = ' '.join(['blastp', '-outfmt', '6', '-query', query, '-db', db , '-max_target_seqs','10000','-evalue',evalue,'-out', blast_out])
    blast = subprocess.Popen(cmd,shell=True)
    blast.wait()


def muti_blast_search(query_paths_genome,db,evalue,working_dir,thread):
    pros = mp.Pool(processes=int(thread))
    path = os.path.join(working_dir,'blast_out')
    try:
        os.mkdir(path)
    except FileExistsError:
        pass
    for query_path in query_paths_genome:
        for query_path_gene in query_path:
            pros.apply_async(func=blast_search, args=(query_path_gene,db,evalue,path,))
    pros.close()
    pros.join()


def handle_blast_out(blast_outdir, sequences_inf, genome_sequences_pair, locus_pair, refer_pair, coverage, similarity):
    """ handle result of sequences searching based on the identity and coverage thread set for user"""
    blast_out = {}
    blast_out_list = os.listdir(blast_outdir)
    print(len(list(genome_sequences_pair.keys())))
    for blast_out_file in blast_out_list:
        try:
            path = os.path.join(blast_outdir,blast_out_file)
            blast = inputfile(path)
            for line in blast:
                q = line.strip('\n').split('\t')[0]
                q_genome = q.lstrip('>').split('_')[0]
                r = line.strip('\n').split('\t')[1]
                match_length = int(line.strip().split('\t')[3])
                identity = line.strip('\n').split('\t')[2]
                evalue = line.strip('\n').split('\t')[-2]
                refer_coverage = '%.4f' % (float(match_length) / len(refer_pair['>%s'%r]))
                query_coverage = '%.4f' % (float(match_length) / len(genome_sequences_pair[q]))
                if float(refer_coverage) >= float(coverage) and float(query_coverage) >= float(coverage):
                    if float(identity) >= float(similarity):
                        mapping_speices = str(sequences_inf[r]['speceis name'])
                        mapping_islate = str(sequences_inf[r]['isolate name'])
                        mapping_accession = str(sequences_inf[r]['accessions'])
                        mapping_product = str(sequences_inf[r]['product'])
                        hit_inf = [q,r,identity,mapping_speices,mapping_product,mapping_islate,mapping_accession,evalue]
                        if q_genome not in blast_out.keys():
                            blast_out[q_genome] = {q: [hit_inf]}
                        else:
                            if q not in blast_out[q_genome].keys():
                                blast_out[q_genome][q] = [hit_inf]
                            else:
                                if len(blast_out[q_genome][q]) < 5:
                                    blast_out[q_genome][q].append(hit_inf)
                                else:
                                    minHit = sorted(blast_out[q_genome][q], key=lambda X: float(X[2]))[0]
                                    if float(hit_inf[2]) > float(minHit[2]):
                                        blast_out[q_genome][q].remove(minHit)
                                        blast_out[q_genome][q].append(hit_inf)
                                    else:
                                        pass
        except UnicodeDecodeError:
            print(path)
    return blast_out


def read_compelement(location):
    if'complement' in location:
        return '-'
    else:
        return '+'


def product(seq):
    if not seq:
        yield []
    else:
        for element in seq[0]:
            for rest in product(seq[1:]):
                yield [element] + rest


def sorted_hits(hit_candidates_list):
    best_hit = []
    phage_hit = {}
    for candidates in hit_candidates_list:
        best = sorted(candidates, key=lambda X: float(X[2]), reverse=True)[0]
        best_hit.append(best)
        for hit_inf in candidates:
            phage = hit_inf[6]
            q = hit_inf[0]
            s = float(hit_inf[2])
            if phage not in phage_hit.keys():
                phage_hit[phage] = {q: (s, hit_inf)}
            else:
                if q not in phage_hit[phage].keys():
                    phage_hit[phage][q] = (s, hit_inf)
                else:
                    print(phage_hit[phage][q])
                    if s > phage_hit[phage][q][0]:
                        phage_hit[phage][q] = (s, hit_inf)
                    else:
                        pass
    phage_hit_sorted, max_numbers, mean_similarity = sorted([
        (phage, len(list(query_inf.keys())), float(np.mean([v for v, h in list(query_inf.values())])))
        for phage, query_inf in phage_hit.items()], key=lambda X: (X[1], X[2]), reverse=True)[0]
    proportion = max_numbers / len(hit_candidates_list)
    return phage_hit_sorted, proportion, best_hit


def genes_arrangement_out(blast_final_out, genes_num, path, region_num, species_inf_dic):
    region_name = 'region%d' % region_num
    if len(blast_final_out) >= int(genes_num):
        print('\n{0:-^100}'.format(region_name))
        predict_phage, max_pro, best_hit = sorted_hits(blast_final_out)
        if predict_phage is not None:
            blast_final_str = '\n'.join(['\t'.join(hit) + '\t' + '\t'.join(species_inf_dic[predict_phage]) + '\t' + '%.4f' % max_pro for hit in best_hit])
            output_file(os.path.join(path, 'region%d.out' % region_num), blast_final_str)
            print('\n'.join([hit[0] + '\t' + '\t'.join(species_inf_dic[predict_phage]) + '\t' +'%.4f' % max_pro for hit in best_hit]))
        else:
            pass
    else:
        pass


def blast_output(filter_blast_out_dic, locus_pair, genome_file_pair, strand, genes_num,outdir, species_inf_dic):
    """final output for searching result based on the genes location"""
    for genome, blast_out in filter_blast_out_dic.items():
        print('\n{0:*^70}'.format(genome_file_pair[genome]))
        path = os.path.join(outdir, genome_file_pair[genome])
        blast_final_out = []
        try:
            os.mkdir(path)
        except FileExistsError:
            pass
        genes_list = list(blast_out.keys())
        blast_out_sorted = sorted(genes_list, key=lambda X: int(X.split('|')[0].split('_')[-1]))
        region_num = 1
        for line_index, query_gene in enumerate(blast_out_sorted):
            hit_condidate_list = blast_out[query_gene]
            contigs = query_gene.split('|')[1].split('_')[0]
            query_strand = read_compelement(query_gene)
            if line_index == 0:
                blast_final_out.append(hit_condidate_list)
            else:
                forward_q = blast_out_sorted[line_index-1]
                forward_contig = forward_q.split('|')[1].split('_')[0]
                forward_site_strand = read_compelement(forward_q)
                forword_query = forward_q.split('|')[0].split('_')[-1]
                query_inf = query_gene.strip('>').split('|')[0].split('_')[-1]
                if contigs == forward_contig:
                    if (int(query_inf) - int(forword_query)) <= 1:
                        if strand:
                            if forward_site_strand == query_strand:
                                blast_final_out.append(hit_condidate_list)
                            else:
                                if len(blast_final_out) >= int(genes_num):
                                    genes_arrangement_out(blast_final_out,genes_num,path,region_num, species_inf_dic)
                                blast_final_out = []
                                blast_final_out.append(hit_condidate_list)
                        else:
                            blast_final_out.append(hit_condidate_list)
                    else:
                        if len(blast_final_out) >= int(genes_num):
                            genes_arrangement_out(blast_final_out,genes_num,path,region_num, species_inf_dic)
                            region_num += 1
                        blast_final_out = []
                        blast_final_out.append(hit_condidate_list)

                else:
                    if len(blast_final_out) >= int(genes_num):
                        genes_arrangement_out(blast_final_out,genes_num,path,region_num, species_inf_dic)
                        region_num +=1
                    blast_final_out = []
                    blast_final_out.append(hit_condidate_list)


if __name__ == '__main__':
    process_start_time = time.time()
    time_string = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(process_start_time))
    print('[%s]Finding Prophage....' % time_string)
    print('Step1:Reading parameters:')
    args = get_parameters()
    input_file = args.q
    evalue = args.e
    coverage = args.c
    thread = args.t
    output_dir = args.o
    similarity = args.s
    input_type = args.type
    strand = args.strand
    genes_num = args.r
    remove_cash = args.clear
    db_path = os.path.join('/home/zhangzhen/Virusdb/',args.d)
    print('Step2:prepare sequences')
    try:
        os.mkdir(db_path)
    except FileExistsError:
        pass
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
    working_dir = os.path.join(output_dir,'working_dir')
    try:
        os.mkdir(working_dir)
    except FileExistsError:
        pass
    sequences_inf,refer_path,refer_pair, species_inf_dic = make_mongo_db_dic(db_path)
    if input_type == 'prot':
        genome_locus_pair,genome_sequences_pair,recode_genome_path,genome_file_pair = prepare_genome_file(input_file,working_dir)
    elif input_type == 'nucle':
        genome_file_pair,annotation_dic = annotation_nucle_genome(input_file,working_dir,thread)
        genome_locus_pair,genome_sequences_pair,recode_genome_path = read_annotation_file(annotation_dic)
    else:
        print('please check your type input')
        sys.exit(0)
    print('Step3:make sequences search')
    # ###make blast search
    DB = make_blast_db(refer_path)
    muti_blast_search(recode_genome_path,DB,evalue,working_dir,thread)
    print('Step4:search result analysis')
    ### handle blast output file
    blast_outdir = os.path.join(working_dir, 'blast_out')
    filter_blast_out_dic = handle_blast_out(blast_outdir,sequences_inf,genome_sequences_pair,genome_locus_pair,refer_pair,coverage,similarity)
    blast_output(filter_blast_out_dic,genome_locus_pair,genome_file_pair,strand,genes_num,output_dir, species_inf_dic)
    process_end_time = time.time()
    print('Whole Process %s' % TimeUsedComputation(process_start_time, process_end_time))
    if remove_cash:
        if os.path.isdir(working_dir):
            shutil.rmtree(working_dir)



