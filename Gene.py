# coding=utf-8
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


class Gene():

    def __init__(self, name):
        self._name = name
        self._fullSequence = ""
        self._size = 0
        self._transcripts = {}
        self._nbTranscripts = 0

    @property
    def name(self):
        return(self._name)

    @name.setter
    def name(self):
        raise SyntaxError("Name should not be modified")

    @property
    def transcripts(self):
        return(self._transcripts)

    def addTranscripts(self, name, sequence):
        if name not in self._transcripts:
            self._transcripts[name] = sequence
            self._nbTranscripts += 1
        else:
            msg = "Duplicate entry in gene " + self._name
            msg += "\n value = " + name
            raise ValueError(msg)

    @property
    def tSize(self):
        return(self._nbTranscripts)

    @tSize.setter
    def tSize(self):
        raise SyntaxError("Size should not be set manually")

    @property
    def length(self):
        return(self._size)

    @length.setter
    def length(self):
        raise SyntaxError("Length should not be set manually")

    @property
    def sequence(self):
        return(self._fullSequence)

    @sequence.setter
    def sequence(self, seq):
        self._sequence = seq
        self._size = len(seq)
