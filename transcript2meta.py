#coding=utf-8
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
from EnsemblFasta import EnsemblFasta as ef
import pysam
import sys
import os

SEP = "\t"  # separator


def sam2graph(sam, fasta):
    """ Convert a  BAM / SAM file into a alignement set dictionary"""

    refDict = {}
    targetDict = {}
    G = {}
    with pysam.AlignmentFile(sam) as samfile:
        for read in samfile.fetch():
            if(not read.flag & 4  ):
                # Fetching read and reference names
                targetName = read.query_name
                refName = read.reference_name

                # Checking existence and creating / fetching objects
                if(targetName not in targetDict):
                    target = Read(read.query_name)
                    targetDict[targetName] = target
                else:
                    target = targetDict[targetName]

                # Since one read can be mapped on multiple transcripts,
                # sequences are only stored on "prefered" mapping.
                # Need to account for this when filling target objects
                if(read.query_sequence):
                        target.sequence = read.query_sequence

                if(refName not in refDict):
                    reference = Read(refName)
                    try:
                        refGene = fasta.genes[fasta.transcriptMap[refName]]
                    except KeyError:
                        print(refName)
                        exit()
                    reference.sequence = refGene.transcripts[refName]
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


dFasta = "/home/cube/Documents/DATA/Références/"
dFasta += "Mus_musculus.GRCm38.cdna.all_ch14.fasta"

fastaFile = sys.argv[1] if len(sys.argv) > 1 else dFasta
samFile = sys.argv[2] if len(sys.argv) > 2 else "./test.bam"

# parsing fasta
fasta = ef(fastaFile)
graph = sam2graph(samFile, fasta)


if "tmp" not in os.listdir(os.getcwd()):
    os.mkdir("tmp")

filePath = os.path.basename(samFile)[:-4]
filePath += "_gene2readlist.txt"

out = open( filePath , "w")
print("OUTPUT FILE:")
print(filePath)
for gene in fasta.genes:
    # if gene not in os.listdir(os.getcwd() + "/tmp"):
    #     os.mkdir("./tmp/" + gene)
    printed = 0
    block = "G" + SEP + gene + "\n"
    for i, alt in enumerate(fasta.genes[gene].transcripts.keys()):
        if(alt in graph.keys()):
            block += "T" + SEP + alt + "\n"
            printed += 1
            altRef = next(iter(graph[alt])).ref

            data = sorted(graph[alt], key=lambda x: (
                x.read.size, x.startRef))[::-1]
            # ga.drawData(data, altRef, True, False, alt, "tmp/" + gene + "/")
            # ga.drawData(data, altRef, False, True, alt, "tmp/" + gene + "/")
            for align in data:
                block += "R" + SEP + align.read.name + "\n"
    if(printed != 0):
        out.write(block)

out.close()
