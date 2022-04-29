def reverse_complement_strand(strand):
    """
    A Helper function to reverse complement the strand

    :param strand: a sequence of genome

    :return: a string that reverse-complement the strand
    """

    # initialize an empty string
    reversed = ""
    # Since reverse-complement reads forward strand from the end to the beginning,
    # set an iterator as the length of the forward strand
    i = len(strand) - 1

    # A while loop to start iteration
    while i >= 0:
        # If the letter is "A", add "T" to the reversed string
        if strand[i] == "A":
            reversed += "T"
            # If the letter is "T", add "A" to the reversed string
        elif strand[i] == "T":
            reversed += "A"
        # If the letter is "G", add "C" to the reversed string
        elif strand[i] == "G":
            reversed += "C"
        # If the letter is "C", add "G" to the reversed string
        elif strand[i] == "C":
            reversed += "G"

            # decreasing iterator every time
        i = i - 1
    # Return reversed string
    return reversed


def reverse_strand(strand):
    """
    A Helper function that reverse the strand

    :param strand: a sequence of genome to search for gene inside

    :return: a string that reverse the strand
    """

    # initialize an empty string
    reversed = ""

    # A for loop to iterate the forward-strand
    for i in range(len(strand)):
        # If the letter is "A", add "T" to the reversed string
        if strand[i] == "A":
            reversed += "T"
        # If the letter is "T", add "A" to the reversed string
        elif strand[i] == "T":
            reversed += "A"
        # If the letter is "G", add "C" to the reversed string
        elif strand[i] == "G":
            reversed += "C"
        # If the letter is "C", add "G" to the reversed string
        elif strand[i] == "C":
            reversed += "G"
    # return reversed strand
    return reversed


# A function for transcribe DNA to mRNA
# The function takes an input gene as a parameter (a string)
# It returns a string of mRNA
def mRNA(gene):
    """
    A Helper function to find mRNA

    :param gene: to find mRNA

    :return: a string of mRNA
    """

    # An empty string for mRNA
    mRNA = ""

    # A loop to go through the input gene
    for i in range(len(gene)):
        # If the letter is "A", add "U" to mRNA string
        if gene[i] == "A":
            mRNA += "U"
        # elif the letter is "T", add "A" to mRNA string
        elif gene[i] == "T":
            mRNA += "A"
        # elif the letter is "G", add "C" to mRNA string
        elif gene[i] == "G":
            mRNA += "C"
        # elif the letter is "C", add "G" to mRNA string
        elif gene[i] == "C":
            mRNA += "G"

    # Return the mRNA result
    return mRNA


def BeginEnd(mRNA):
    """
    A Helper function finds the indexes of start position and stop position in mRNA to start translating
    and return the indexes of the start position and stop position

    :param mRNA: a sequence of genome to search for gene inside

    :return: indexes of the start position and stop position
    """

    # A counter for start position
    begin = mRNA.find("AUG")

    # After the start position was found, start at the start position to find stop codons
    end = begin
    # The length of mRNA we are going to look through
    j = len(mRNA) - 2

    # A while loop to find one of three stop codons
    while end < j:
        # If one of three codons is found, stop loop
        if mRNA[end:end + 3] == "UGA" or mRNA[end:end + 3] == "UAG" or mRNA[end:end + 3] == "UAA":
            break
            # Else keep searching and add three each time to end index variable
        else:
            end += 3

        # Return the indexes of start position and stop position
    return begin, end


# A function for finding gene in the genome
# with two bonus
# The function takes three inputs, one is genome, one is
# The function outputs a string of gene, and it starts place and stop place in the input genome
def gene_detection(genome, enable_terminator, enable_pribnow):
    """
    A Helper function to detect gene in the genome and if found,
    return a list of genomes and the information of each genome as a dictionary

    :param genome: a sequence of genome to search for gene inside
    :param enable_terminator: check whether enableTerminator is enabled or not
    :param enable_pribnow: check whether Pribnow is enabled or not
    :return: a list of genomes, and a dictionary containing the information of each genome
    """
    # initialize a list of genes, beginGenes, and dictionary InformationGene
    gene = []
    beginGenes = []
    informationOnGene = {}

    # initialize probabilistic matrices to score potential promoter sequences
    scores = {
        "A": [-2.6, 1.86, -5.02, 1.86, 1, 46, -4.98],
        "C": [-1.08, -9999999999, -3.28, 0, 0, -3.24],
        "G": [-2.43, -5.6, -5.6, -4.28, 0, -3.98],
        "T": [1.67, -1.47, 1.94, -1.7, 0.32, 1.93]
    }

    # find the start of the gene by finding the pribnow box in the genome
    # a while loop to find the start of the gene
    i = 0
    while i < len(genome) - 5:
        score = scores[genome[i]][0] + scores[genome[i + 1]][1] + scores[genome[i + 2]][2] \
                + scores[genome[i + 3]][3] + scores[genome[i + 4]][4] + + scores[genome[i + 5]][5]

        # check if the score > 7 and if enablePribnow is enabled
        if score > 7 and enable_pribnow:
            i += 13
            # put the start of gene inside the list
            beginGenes.append(i)

        # if enablePrinow is not enabled, check the genome if it's "TATAAT"
        elif genome[i:i + 6] == "TATAAT" and not enable_pribnow:
            i += 13
            # put the start of gene inside the list
            beginGenes.append(i)
        # else i += 1
        else:
            i += 1

    # Find the end of gene by finding the terminator inside the genome
    for num in beginGenes:
        RNA = mRNA(reverse_strand(genome))
        # if enableTerminator is not enabled
        if not enable_terminator:
            # Find the terminating poly-U sequence
            for i in range(num, len(RNA) - 9):
                # if poly-U sequence is found
                if RNA[i:i + 10] == "UUUUUUUUUU":
                    sequence1 = RNA[i - 20:i - 12]
                    sequence2 = RNA[i - 8:i]
                    # if 2 rich GC hairpin sequences are reverse-complement
                    if sequence1 == reverse_complement_strand(sequence2):
                        terminator = i + 10
                        # if found, put them in the list gene, and informatino about gene in dictionary
                        gene.append(genome[num:terminator])
                        informationOnGene[len(gene)] = [num, terminator]
                        break

        # if enableTerminator is not enabled
        else:
            # Find the terminating poly-U sequence
            for i in range(num, len(RNA) - 14):
                mismatch = 0
                polyU = RNA[i:i + 15]
                countU = 0
                # check the sequence polyU if they contains enough U inside
                for U in polyU:
                    if U == "U":
                        countU += 1
                if countU >= 10:
                    """
                    a while loop to check 2 sequences before polyU
                    because there exists polyU terminating sequence,
                    check if the 2 hairpin loop sequences before exists
                    with length varies from 4 to 12, and 1 mismatch at most
                    """
                    a = 4
                    while a <= 12:
                        sequence1 = RNA[i - a:i]
                        sequence2 = RNA[i - a - a - 4:i - a-4]
                        # check if the nucleotide in sequence1 and sequence2 are reverse complement
                        # if not, mismatch += 1
                        # if mismatch > 1 then stop the loop
                        for n in range(len(sequence1)):
                            if sequence1[n] == "C":
                                if sequence2[-n - 1] != "G":
                                    if mismatch == 1:
                                        break
                                    mismatch += 1
                            elif sequence1[n] == "G":
                                if sequence2[-n - 1] != "C":
                                    if mismatch == 1:
                                        break
                                    mismatch += 1
                            else:
                                if mismatch == 1:
                                    break
                                mismatch += 1
                        # if mismatch == 1 then break the loop
                        if mismatch == 1:
                            break
                        # else increase the length
                        a += 1
                    # if mismatch is <= 1, put gene in genome in the list of gene
                    # put information of the gene, position of start and end into informationOnGene
                    terminator = i + 15
                    gene.append(genome[num:terminator])
                    informationOnGene[len(gene)] = [num, terminator]
                    break
    # return list of genes, and informationOnGene
    return gene, informationOnGene

def protein(RNA):
    """
    a function for translating protein from RNA
    :param RNA: a string
    :return: a string of protein after translation
    """
    # A dictionary contains the translation codes
    translation = {"A": {"A": {"A": "K","G": "K","C": "N","U": "N"},
                         "C": {"A": "T","G": "T","C": "T","U": "T"},
                         "G": {"A": "R","G": "R","C": "S","U": "S"},
                         "U": {"A": "I","G": "M","C": "I","U": "I"}},
                   "G": {"A": { "A": "E","G": "E","C": "D","U": "D"},
                         "C": {"A": "A","G": "A","C": "A","U": "A"},
                         "G": {"A": "G","G": "G","C": "G","U": "G"},
                         "U": {"A": "V","G": "V","C": "V","U": "V"}},
                   "C": {"A": {"A": "Q","G": "Q","C": "H","U": "H"},
                         "C": {"A": "P","G": "P","C": "P","U": "P"},
                         "G": {"A": "R","G": "R", "C": "R","U": "R"},
                         "U": {"A": "L","G": "L","C": "L", "U": "L"}},
                   "U": {"A": {"A": "*","G": "*","C": "Y","U": "Y"},
                         "C": {"A": "S","G": "S","C": "S","U": "S"},
                         "G": {"A": "*","G": "W","C": "C","U": "C"},
                         "U": {"A": "L","G": "L","C": "F","U": "F"}}
                 }

    # call the BeginEnd helper function to get the start position and stop position
    start_point, end_point = BeginEnd(RNA)
    # initialize an empty string
    protein = ""

    # A while loop to start translating
    while start_point < end_point - 2:
        # three letters are one group
        group = RNA[start_point:start_point + 3]
        # based on the letters in the group, access the dictionary that is initialized to get the letter for protein
        translate = translation[group[0]][group[1]][group[2]]

        # If the loop reaches one of three stop codons, stop the loop
        if translate == "*":
            break

        # otherwise, append the letter to protein and increase start_point by 3 each times
        else:
            protein += translate
            start_point += 3

    # return result
    return protein


# A helper function for reading genome File
# The function takes a filename as a parameter
# The function returns a string of the genome in the file
def genomeStrand(filename):
    """
    a helper function for reading genome file
    :param filename: a filename
    :return: a string of the genome in the file
    """

    # initialize an empty string
    genome = ""

    # open file
    with open(filename, "r") as f:
        # read each line
        for line in f:
            # skip description line
            if line[0] != ">":
                # append genome into the string
                genome += line.strip("\n")
        # return the result
    return genome


# A helper function to output file
# The function takes 5 parameters: input file, forwardgenome, reversegenome, info1 from forwardgenome and info2 from reversegenome
# The function outputs three files that contain gene info, mRNA info and protein info
def outputfile(file, ForwardGenome, ReverseGenome, info1, info2):
    """
    a function to output file
    :param file: read in a filename
    :param ForwardGenome: a sequence of genome that is read forward
    :param ReverseGenome: a sequence of genome that is read from reverse
    :param info1: gene name
    :param info2: gene name
    :return:
    """

    # initialize file names for three output files
    gene = file + "g.fasta"
    rna = file + "r.fasta"
    translate = file + "p.fasta"

    # create 3 output files for gene, mRNA and protein
    gene_f = open(gene, "w")
    rna_f = open(rna, "w")
    trans_f = open(translate, "w")

    # write info into the file from the forwardgenome
    for i in range(len(ForwardGenome)):
        # Call ReverseStrand helper function to get mRNA strand
        RNA = mRNA(reverse_strand(ForwardGenome[i]))
        # Call BegianEnd helper function to get the start and stop position on mRNA for protein file
        first, second = BeginEnd(RNA)
        # store the gene postions
        firstG, secondG = info1[i][0], info1[i + 1][1]

        # write info into files
        trans_f.write("> aaaAp " + str(first + 1) + ":" + str(second) + ":+" + "\n")
        gene_f.write("> aaaA " + str(firstG) + ":" + str(secondG) + ":+" +
                     "\n")
        rna_f.write("> aaaAr " + str(1) + ":" + str(secondG - firstG) +
                    ":+" + "\n")

        # splitting the forward gene strand into 80 characters a line
        num1 = len(ForwardGenome[i]) // 80
        # The rest of characters after splitting the strings into 80 characters each line
        num2 = len(ForwardGenome[i]) % 80

        # write gene strand and mRNA strand into the file
        # each line has 80 characters
        for j in range(num1):
            gene_f.write(ForwardGenome[i][j * 80:(j + 1) * 80] + "\n")
            rna_f.write(RNA[j * 80:(j + 1) * 80] + "\n")

        # splitting the protein strand into 80 characters a line
        num3 = len(protein(RNA)) // 80
        # The rest of characters after splitting the strings into 80 characters each line
        num4 = len(protein(RNA)) % 80

        # write protein strand into the file
        # each line has 80 characters
        for k in range(num3):
            trans_f.write(protein(RNA)[k * 80:(k + 1) * 80] + "\n")

        # After finish each line with 80 characters,
        # write the rest of the letters into the file
        trans_f.write(protein(RNA)[len(protein(RNA)) - num4:] + "\n")
        gene_f.write(ForwardGenome[i][len(ForwardGenome[i]) - num2:] + "\n")
        rna_f.write(RNA[len(RNA) - num2:] + "\n")

    # write info into the file from the reversegenome
    for i in range(len(ReverseGenome)):
        # Call ReverseStrand helper function to get mRNA strand
        RNA = mRNA(reverse_strand(ReverseGenome[i]))
        # Call BegianEnd helper function to get the start and stop position on mRNA for protein file
        first, second = BeginEnd(RNA)
        # store the gene postions
        firstG, secondG = info2[i + 1][0], info2[i + 1][1]

        # write info into files
        trans_f.write("> aaaAp " + str(first + 1) + ":" + str(second + 1) +
                      ":+" + "\n")
        gene_f.write("> aaaA " + str(firstG) + ":" + str(secondG + 1) + ":-" +
                     "\n")
        rna_f.write("> aaaAr " + str(1) + ":" + str(firstG - secondG + 1) +
                    ":+" + "\n")

        # splitting the reverse strand into 80 characters a line
        num1 = len(ReverseGenome[i]) // 80
        # The rest of characters after splitting the strings into 80 characters each line
        num2 = len(ReverseGenome[i]) % 80

        # write gene strand and mRNA strand into the file
        # each line has 80 characters
        for j in range(num1):
            gene_f.write(ReverseGenome[i][j * 80:(j + 1) * 80] + "\n")
            rna_f.write(RNA[j * 80:(j + 1) * 80] + "\n")

        # splitting the protein strand into 80 characters a line
        num3 = len(protein(RNA)) // 80
        # The rest of characters after splitting the strings into 80 characters each line
        num4 = len(protein(RNA)) % 80

        # write protein strand into the file
        # each line has 80 characters
        for k in range(num3):
            trans_f.write(protein(RNA)[k * 80:(k + 1) * 80] + "\n")

        # After finish each line with 80 characters,
        # write the rest of the letters into the file
        trans_f.write(protein(RNA)[len(protein(RNA)) - num4:] + "\n")
        gene_f.write(ReverseGenome[i][len(ReverseGenome[i]) - num2:] + "\n")
        rna_f.write(RNA[len(RNA) - num2:] + "\n")


# A real function to combine all the helper functions to predict protein from genome
# The function takes a string of file name and two boolean values of enableTerminator and enablePribow
# The function returns the output into three files (gene, mRNA, protein) with their infos
def predict(filename, enableTerminator, enablePribnow):
    """
    a function to combine all the helper functions to predict proteins from genome

    :param filename: a filename
    :param enableTerminator:  boolean expression True or False if want to enable to find Terminator or use default Terminator
    :param enablePribnow: boolean expression True or False if want to enable to find pribnow or use default pribnow
    :return:
    """
    # read the genome from the input file
    genome = genomeStrand(filename)

    # get the forward strand
    ForwardGene, info1 = gene_detection(genome, enableTerminator, enablePribnow)
    # get the reverse strand
    ReverseGene, info2 = gene_detection(reverse_complement_strand(genome),
                                        enableTerminator, enablePribnow)

    # for reverse strand,
    # get the position backward
    if len(ReverseGene) != 0:
        for i in range(1, len(ReverseGene) + 1):
            info2[i][0] = len(genome) - info2[i][0]
            info2[i][1] = len(genome) - info2[i][1]

    # get a new file name from the input file
    string = filename[0:len(filename) - 5]
    # call output function to have output file
    outputfile(string, ForwardGene, ReverseGene, info1, info2)


predict("mixed_overlapping_1.r/mixed_overlapping_1.fasta", False, False)
