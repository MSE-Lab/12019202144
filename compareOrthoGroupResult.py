import sys
import os

def inputfile(fileName):
	with open(fileName) as f:
		data = f.readlines()
	return data 


def outputfile(fileNames,content):
	with open(fileNames,'w') as f:
		f.writelines(content)


def GetSOG(OrthoGroupFile,Source,GenomeNums=80):
	SOG = []
	OGfile = inputfile(OrthoGroupFile)
	for line in OGfile:
		OG = line.strip('\n').split('\t')[0] 
		geneList = line.strip('\n').split('\t')[1].strip(' ').split(' ')
		genomes = set()
		for gene in geneList:
			genome_id = gene.split('|')[0]
			genomes.add(genome_id)
			if len(genomes) == GenomeNums and len(geneList) == GenomeNums:
				SOG.append(set(geneList))
			else:
				pass
	return SOG


def ComparisonOrthoGroups(OGsetQ,OGsetR):
	overlap = []
	OGSets = []
	for q in OGsetQ:
		matchs = []
		for r in OGsetR:
			overlap.append(len(q&r)/80)
			if	len(q & r) >= 50:
				matchs.append(r)
		if len(matchs) > 0:
			OGSets.append(q)
	print(len(OGSets))
	return OGSets, overlap


def SortedCompare(fileSet):
	compareSet = []
	for qlabels in fileSet:
		for rlabels  in fileSet:
			if qlabels != rlabels:
				if (rlabels, qlabels) not in compareSet:
					compareSet.append((qlabels,rlabels))
	return compareSet


def main(inputdir):
	fileSet = SortedCompare(os.listdir(inputdir))
	CompareResults = {}
	overlapDic = {}
	for query,refer in fileSet:
		queryOG = GetSOG(os.path.join(inputdir,query),query)
		referOG = GetSOG(os.path.join(inputdir,refer),refer)
		OGsetComparsion,overlaps = ComparisonOrthoGroups(queryOG,referOG)
		CompareResults['%s-%s'%(query,refer)] = OGsetComparsion
		overlapDic['%s-%s'%(query,refer)] = overlaps
	fileSet1 = os.listdir(inputdir)
	query1 = fileSet1[0]
	refer1 = fileSet1[1]
	try:
		pairwiseSet = CompareResults['%s-%s'%(query1,refer1)] 	
	except KeyError:
		pairwiseSet = CompareResults['%s-%s'%(refer1,query1)]
	Third = GetSOG(os.path.join(inputdir,fileSet1[2]),fileSet1[2])
	triple = ComparisonOrthoGroups(pairwiseSet,Third)[0]
	outputResults(CompareResults,triple)
	string1 = ''
	for num,k in enumerate(overlapDic):
		for s in overlapDic[k]:
			string1+='%s%d\t%.4f\n'%(k,num,s)
	outputfile('overlap.txt',string1)
	

def outputResults(pairWise,Triple):
	Static = ''
	for k,v in pairWise.items():
		Static+='%s\t%d\n'%(k,len(v))
	Static = Static+'Triple'+'\t'+'%d'%(len(Triple))
	outputfile('Stat.txt',Static)
		
			
if __name__ == '__main__':
	inputDir = sys.argv[1:][0]
	main(inputDir)

			
			
