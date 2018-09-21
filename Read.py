# coding=utf-8
#############################################################################
#                                                                           #
#    General Read description                                               #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#       Version: 0.2                                                        #
#       Author: Quentin Bonenfant                                           #
#               quentin.bonenfant@gmail.com                                 #
#############################################################################


class Read():

    def __init__(self, name):

        self._name = name
        self._sequence = None
        self._size = 0

    # Read name
    @property
    def name(self):
        return(self._name)

    @name.setter
    def name(self, n):
        self._name = n

    # SIZE
    @property
    def size(self):
        return(self._size)

    @size.setter
    def size(self, s):
        raise(RuntimeError(
            """\nYou can't explicitly set read size.\n
                 Try changing sequence instead.\n"""))

    # SEQUENCE
    @property
    def sequence(self):
        return(self._sequence)

    @sequence.setter
    def sequence(self, seq):
        self._sequence = seq
        self._size = len(seq)
