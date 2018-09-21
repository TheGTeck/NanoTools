# coding=utf-8
#############################################################################
#                                                                           #
# Find the metareads from a mapping over reference transcriptome            #
# from ENSEMBL database                                                     #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#       Version: 0.2                                                        #
#       Author: Quentin Bonenfant                                           #
#               quentin.bonenfant@gmail.com                                 #
#############################################################################


from Read import Read
from Align import Align
from GraphAlign import GraphAlign as ga
import pysam
import re
import sys
import os

SEP = "\t"  # separator
# Rule to detect reference gene for alternative transcripts and
# their position in the reference genome
geneIdRule = re.compile(
    r"\>(ENSMUST\d+\.\d+) cdna chromosome\:\w+\:\d+\:\d+\:\d+\:\-?\d gene\:(\w+\.\d+) gene_biotype\:\.*")


def parseFasta(fastaRef):
    """Associate an alternative transcript to it's original gene,
    and  an id / sequence dict"""

    backReference = {}
    fastaDict = {}
    seq = ""
    prevId = ""
    with open(fastaRef, 'r') as f:
        for line in f:
            if ">" == line[0]:
                found = geneIdRule.search(line)
                alternate = found.group(1)
                ref = found.group(2)
                backReference[alternate] = ref
                if(prevId and seq):
                    fastaDict[prevId] = seq
                    seq = ""
                prevId = alternate
            else:
                seq += line.rstrip("\n")
    fastaDict[prevId] = seq
    return(backReference, fastaDict)


def sam2graph(sam):
    """ Convert a  BAM / SAM file into a alignement set dictionary"""

    refDict = {}
    G = {}
    with pysam.AlignmentFile(sam) as samfile:
        for read in samfile.fetch():
            # Read and Reference data
            target = Read(read.query_name)
            target.sequence = read.query_sequence

            refName = read.reference_name

            if(refName not in refDict):
                reference = Read(refName)
                reference.sequence = fastaDict[refName]
                refDict[refName] = reference
                G[refName] = set()
            else:
                reference = refDict[refName]

            # Alignment data
            s1 = read.query_alignment_start
            e1 = read.query_alignment_end
            s2 = read.reference_start
            e2 = read.reference_end
            cigar = read.cigarstring
            relPos = "-" if read.is_reverse else "+"
            align = Align(target, reference, s1, e1, s2, e2, cigar, relPos)

            G[refName].add(align)
    return(G)


def mappedAltToGene(backReference, graph):
    """Make a reverse version of the reference
    where genes index the alternative transcript,
    if there is a mapping on any isoform."""

    reference = {}
    for alt in graph.keys():
        gene = backReference[alt]
        if gene not in reference.keys():
            reference[gene] = set()
        reference[gene].add(alt)

    return(reference)


dFasta = "/home/cube/Documents/DATA/Références/"
dFasta += "Mus_musculus.GRCm38.cdna.all_ch14.fasta"

samFile = sys.argv[1] if len(sys.argv) > 1 else "./test.bam"
fastaFile = sys.argv[2] if len(sys.argv) > 1 else dFasta
global fastaDict
backReference, fastaDict = parseFasta(fastaFile)
graph = sam2graph(samFile)
references = mappedAltToGene(backReference, graph)


if "tmp" not in os.listdir(os.getcwd()):
    os.system("mkdir tmp")


for gene, alternates in references.items():
    # print("\nGENE: " + gene)
    if gene not in os.listdir(os.getcwd() + "/tmp"):
        os.system("mkdir ./tmp/" + gene)
    for alt in alternates:

        altRef = next(iter(graph[alt])).ref
        data = sorted(graph[alt], key=lambda x: (
            x.read.size, x.startRef))[::-1]
        ga.drawData(data, altRef, True, False, alt, "tmp/" + gene + "/")
        # ga.drawData(data, altRef, False, True, alt, "tmp/" + gene + "/")
