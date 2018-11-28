# coding=utf-8
#############################################################################
#                                                                           #
#    Generate images showing alignement between reference and / or reads.   #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#                                                                           #
#       Version: 0.2                                                        #
#       Author: Quentin Bonenfant                                           #
#               quentin.bonenfant@gmail.com                                 #
#############################################################################

from PIL import Image, ImageDraw
import re


class GraphAlign():

    # Initiating Globals
    global BOXHEIGHT, HRMARGIN, BOTTOMMARGIN, SPACER, OFFSET
    global GRADSIZE, ABSCISSE, ORDINATE, CIGARPATTERN
    global BLOCKSIZE, BLOCKSPACE

    # Height of the boxes
    BOXHEIGHT = 15
    # Space between boxes
    SPACER = 2
    # Left margin
    HRMARGIN = 10
    # Bottom margin
    BOTTOMMARGIN = 10
    # Abcisse height from bottom margin
    ABSCISSE = 20
    # Ordonee position from the left margin
    ORDINATE = 325
    # Graduation OFFSET
    OFFSET = 50
    # Size of the graduation
    GRADSIZE = 4
    # CIGAR parttern for splicing
    CIGARPATTERN = re.compile('([MIDNSHPX=])')
    # Sequence block size
    BLOCKSIZE = 20
    # Space between blocks
    BLOCKSPACE = 10

    
    @staticmethod
    def drawCigarAlign(draw, aln, height):
        """ Draw the aligned sequence of a read to a reference"""
        # Setting 'constants' for this function
        center = height + BOXHEIGHT // 2 + 1  # centerline height
        col = "BLUE" if aln.direction == "+" else "RED"
        charWidth = draw.textsize(" ")[0]

        drawPos = aln.startRef * charWidth + HRMARGIN + ORDINATE
        blockPos = aln.startRef % BLOCKSIZE
        spaces = BLOCKSPACE

        # if blocksize is set, add spaces before first nucleotide
        if(BLOCKSIZE != 0):
            drawPos += (int((aln.startRef) / BLOCKSIZE) * BLOCKSPACE) * charWidth
        ncCol = "DARKGREY"
        insertFirst = False
        for readPos, mapping in enumerate(aln.cigarGen):
            # unpacking mapping data
            op, refChar, readChar = mapping

            # avoinding to print soft clipping
            if(op != "S"):
                nc = readChar
                if(op in "MX="):
                    if( refChar == readChar):
                        ncCol = "DARKGREEN"
                    else:
                        ncCol = "RED"
                elif(op == "D"):
                    ncCol = "BLUE"
                elif(op == "I"):
                    ncCol = "ORANGE"
                    if(blockPos % BLOCKSIZE == 0):
                        insertFirst = True
                    spaces -= 1
                    if(spaces == 0):
                        nc = "+"
                    elif(spaces < 0):
                        nc = ""
                        spaces = 0

                if(blockPos == 0 and readPos != aln.startRead):
                    if(not insertFirst):
                        drawPos += spaces * charWidth
                        spaces = BLOCKSPACE
                    else:
                        insertFirst = False

                draw.text((drawPos, center),
                          nc,
                          fill=ncCol)
                if(nc != ""):
                    drawPos += charWidth
                if(op != "I"):
                    blockPos += 1
                blockPos %= BLOCKSIZE

        # Read name
        draw.text((HRMARGIN, center), aln.read.name,
                  fill=col)


    @staticmethod
    def drawBoxes(draw,aln,height):
        """ Draw the reads on a reference as boxes"""
        # Setting 'constants' for this function
        center = height + BOXHEIGHT // 2 + 1  # centerline height
        col = "BLUE" if aln.direction == "+" else "RED"
        charWidth = draw.textsize(" ")[0]

        drawPos = aln.startRef * charWidth + HRMARGIN + ORDINATE
        blockPos = aln.startRef % BLOCKSIZE
        spaces = BLOCKSPACE

        # if blocksize is set, add spaces before first nucleotide
        if(BLOCKSIZE != 0):
            drawPos += (int((aln.startRef) / BLOCKSIZE) * BLOCKSPACE) * charWidth
        ncCol = "DARKGREY"
        insertFirst = False
        for readPos, mapping in enumerate(aln.cigarGen):
            pass


    @staticmethod
    def drawData(data, reference, save=False, show=True, name="unknown", path="./", cigar = False):
        """Create and save/display the reads and overlaps."""

        length = 1

        # Height of the image, depend on the number of data point and htop
        height = len(data) * (BOXHEIGHT + SPACER) + ABSCISSE + BOTTOMMARGIN
        height += BOXHEIGHT + SPACER  # adding one line to avoid collisions
        if(reference):
            height += BOXHEIGHT + SPACER

        # Creating image
        im = Image.new("RGB", (length, height), "WHITE")

        # Creating drawing layer
        draw = ImageDraw.Draw(im)

        # char width
        charWidth = draw.textsize(" ")[0]

        # Finding number of spaces
        nbSpaces = int(reference.size / BLOCKSIZE) + 1

        # Resizing caneva
        length = (reference.size + nbSpaces * BLOCKSPACE) * charWidth
        length += ORDINATE + HRMARGIN * 2
        im = im.resize((length, height))
        draw = ImageDraw.Draw(im)

        # Drawing axis
        draw.line((HRMARGIN + ORDINATE, height - (BOTTOMMARGIN + ABSCISSE),
                   length + HRMARGIN + ORDINATE,
                   height - (BOTTOMMARGIN + ABSCISSE)), fill=0)

        draw.line((HRMARGIN + ORDINATE, 0, HRMARGIN + ORDINATE,
                   height - (BOTTOMMARGIN + ABSCISSE)), fill=0)
        for i in range(0, (length + 1) // charWidth, OFFSET):
            # Current graduation position
            hPos = HRMARGIN + ORDINATE + i * charWidth
            # Drawing graduations
            draw.line((hPos,
                       height - (BOTTOMMARGIN + ABSCISSE + GRADSIZE),
                       hPos,
                       height - (BOTTOMMARGIN + ABSCISSE) + GRADSIZE), fill=0)

            draw.text((hPos - 5,
                       height - (BOTTOMMARGIN + ABSCISSE) + 5),
                      str(i),
                      fill=0)

        #  Drawing Reference if any
        refSpace = 0  # keeping if reference is printed or not
        if(reference):
            refSpace = 1
            rs = ""
            for i in range(0, reference.size, BLOCKSIZE):
                rs += reference.sequence[i: i + BLOCKSIZE] + " " * BLOCKSPACE

            draw.text((HRMARGIN + ORDINATE,
                       SPACER + BOXHEIGHT // 2 + 1),
                      rs, fill="DARKGREY")
            draw.text((HRMARGIN,
                       SPACER + BOXHEIGHT // 2 + 1),
                      reference.name,
                      fill="DARKGREY")

        # Drawing data
        for i, d in enumerate(data):
            off = BOXHEIGHT + SPACER
            GraphAlign.drawAlign(draw, d, (i + refSpace) * off + SPACER)
        del draw
        if(show):
            im.show()
        if(save):
            im.save(path + str(name) + ".png")
