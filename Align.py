# coding=utf-8
#############################################################################
#                                                                           #
#    Container for alignement between two sequences (reads or reference)    #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#       Version: 0.2                                                        #
#       Author: Quentin Bonenfant                                           #
#               quentin.bonenfant@gmail.com                                 #
#############################################################################


import re


class Align():

    def __init__(self, r1, r2, st1, end1, st2, end2, cig=None, direction="+"):

        self._read = r1         # Read object aligned
        self._startRead = st1   # Start of the mapping on the read
        self._endRead = end1    # End of the mapping on the
        self._ref = r2          # Reference read
        self._startRef = st2    # Starting point of the mapping on reference
        self._endRef = end2     # End point of the mapping on reference
        self._direction = direction  # Relative orientation (read to reference)
        self._cigar = cig       # CIGAR format of the mapping
        self.pattern = re.compile('([MIDNSHPX=])')

    # CIGAR
    @property
    def cigar(self):
        return(self._cigar)

    @cigar.setter
    def cigar(self, cig):
        self._cigar = cig

    # read
    @property
    def read(self):
        return(self._read)

    # ref
    @property
    def ref(self):
        return(self._ref)

    # startRead
    @property
    def startRead(self):
        return(self._startRead)

    # endRead
    @property
    def endRead(self):
        return(self._endRead)

    # startRef
    @property
    def startRef(self):
        return(self._startRef)

    # endRef
    @property
    def endRef(self):
        return(self._endRef)

    # Direction
    @property
    def direction(self):
        return(self._direction)

    @property
    def cigarGen(self):

        # turn cigar into tuple of values
        values = self.pattern.split(self._cigar)[:-1]
        posR1 = self._startRef
        posR2 = 0
        mapping = (None, None, None)  # alignement tuple
        for i in range(0, len(values), 2):
            n = int(values[i])
            op = values[i + 1]
            for j in range(n):
                if(op in "MX="):
                    mapping = (
                        op,
                        self._ref.sequence[posR1],
                        self._read.sequence[posR2])
                    posR1 += 1
                    posR2 += 1
                elif(op in "S"):
                    mapping = (
                        op,
                        "",
                        self.read.sequence[posR2])
                    posR2 += 1
                elif(op == "I"):
                    mapping = (
                        op,
                        "-",
                        self.read.sequence[posR2])
                    posR2 += 1

                elif(op == "D"):
                    mapping = (
                        op,
                        self._ref.sequence[posR1],
                        "-")
                    posR1 += 1
                else:
                    msg = "Operation " + op
                    msg += " is not implemented yet."
                    raise NotImplemented(msg)

                yield(mapping)

    @property
    def paddedSequence(self):
        """ Return the sequences as it has been mapped to the reference """

        sequence = self._read.sequence
        padded = ""

        # turn cigar into tuple of values
        values = self.pattern.split(self._cigar)[:-1]
        pos = 0
        for i in range(0, len(values), 2):
            # /!\ Dirty version, maybe a cleaner version will come
            n = int(values[i])
            op = values[i + 1]
            if(op in "SMX="):
                padded += sequence[pos:pos + n]
                pos += n
            elif(op == "I"):
                padded += "+" * n
                pos += n
            elif(op == "D"):
                padded += "-" * n
            elif(op == "S"):
                # padded += "_" * n
                # pos += n
                pass
            elif(op == "H"):
                padded += "|"
        return(padded)
