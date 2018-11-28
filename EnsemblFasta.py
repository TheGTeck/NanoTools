#coding=utf-8
#############################################################################
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#       Version: 0.2                                                        #
#       Author: Quentin Bonenfant                                           #
#               quentin.bonenfant@gmail.com                                 #
#############################################################################

from Gene import Gene
from collections import defaultdict
import re


class EnsemblFasta():

    # CONSTANTS
    global SEP, GENEIDRULE

    # separator
    SEP = "\t"

    # Rule to detect reference gene for alternative transcripts and
    # their position in the reference genome for ENSEMBL fasta
    RULE = r"\>(ENSMUST\d+\.?\d*) .+ gene\:(ENSMUSG\d+\.?\d*) .+"
    GENEIDRULE = re.compile(RULE)

    def __init__(self, path):

        self._genes = {}
        self._transcripts = defaultdict()
        self.parseFasta(path)

    def parseFasta(self, fastaRef):
        """Associate an alternative transcript to it's original gene.
        Sequences are also stored for futur use"""

        seq = ""
        prevId = ""
        with open(fastaRef, 'r') as f:

            for line in f:
                if ">" == line[0]:
                    # asserting the regex don't fail...
                    found = GENEIDRULE.search(line)
                    if(found):
                        alternate = found.group(1)
                        geneName = found.group(2)
                        self._transcripts[alternate] = geneName
                    else:
                        print("EnsemblFasta: NOT FOUND")
                        print(line)
                        exit()

                    if(prevId and seq):
                        geneName = self._transcripts[prevId]
                        if geneName in self._genes:
                            gene = self._genes[geneName]
                        else:
                            gene = Gene(geneName)
                            self._genes[geneName] = gene

                        gene.addTranscripts(prevId, seq)
                        seq = ""
                    prevId = alternate
                else:
                    seq += line.rstrip("\n")
            gene.addTranscripts(prevId, seq)

    @property
    def transcripts(self):
        return(self._transcripts.keys())

    @transcripts.setter
    def transcripts(self):
        raise SyntaxError("FASTA parsing is unmutable")

    @property
    def transcriptMap(self):
        return(self._transcripts)

    @transcriptMap.setter
    def transcriptMap(self):
        raise SyntaxError("FASTA parsing is unmutable")

    @property
    def genes(self):
        return(self._genes)

    @genes.setter
    def genes(self):
        raise SyntaxError("FASTA parsing is unmutable")