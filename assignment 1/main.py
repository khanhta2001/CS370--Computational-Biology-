"""
This file reads fasta and matrix files to check similarity of 2 sequences
and give sequence alignment based on global alignment Needleman-Wunsch
This file will return a file with end .g.out that output sequence alignment,
scores, matches, mismatches, and gap and also the similarity of between 2 sequences
"""


def read_fasta(filename):
    """
    this function reads a fasta file that contains all sequences information
    and return a list of all sequences

    :param:  filename        : read a fasta filename
    :return: data_collection : return a nested list containing all sequences
    """

    data_collection = []
    # open and read file inside
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] != ">":
                # push the sequence into the list
                list1 = [b for b in line.strip('\n')]
                # push the entire list into the list to turn into a nested list
                data_collection.append(list1)
    return data_collection


def read_matrix(filename):
    """
    This function reads a matrix files that contains all scoring matches
    and return a nested dictionary of all possible matches, as well as mismatches

    :param: filename          : read a matrix file
    :return nested dictionary : a nested dictionary containing all matching and mismatching

    """
    list1 = []
    # open and read file
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] != "#":
                # get all possible sequence A C G T, etc
                list1.append(line.strip(',').split())
    list1[0].insert(0, ' ')
    matrix1 = {}
    # nested for loop to push the matches, mismatches between 2 sequence info
    for i in range(1, len(list1)):
        matrix2 = {}
        for j in range(1, len(list1[i])):
            matrix2[list1[0][j]] = int(list1[i][j])
        matrix1[list1[i][0]] = matrix2
    return matrix1


def sequenceAligningInfo(list1, list2, score, sequence1, sequence2):
    """
    This is a helper function to print out all information about sequence alignment
    To use for global and local alignment

    :param list1    : contains all sequence alignment of 1st sequence
           list2    : contains all sequence alignment of 2nd sequence
           score    : scores from aligning sequences
           sequence1: the first sequence
           sequence2: the second sequence
    :return string  : information on all sequence alignment

    """

    str1 = str2 = ""
    match = mismatch = 0
    # to push all possible matches, mismatches, and gap into the string to print it out
    for i in range(1, len(list1) + 1):
        str2 += list1[-i]
        str1 += list2[-i]
        # calculate the match
        if list1[-i] == list2[-i]:
            match += 1
        # calculate the mismatch
        elif list1[-i] != list2[-i] and list1[-i] != "-" and list2[-i] != "-":
            mismatch += 1
    # check the length of sequence to find the shortest sequence so that we can find similarity
    if len(sequence1) < len(sequence2):
        minLength = len(sequence1)
    else:
        minLength = len(sequence2)
    info = str1 + "\n" + str2 + "\n" + "score: " + str(score) + \
           "\n" + "Match: " + str(match) + "\n" + "Mismatch: " + str(mismatch) + "\n" \
           + "Gap: " + str(len(list1) - match - mismatch) + "\n" + "Similarity: " \
           + str((match / minLength) * 100)
    return info


def optionFile(filename):
    """
    To read in an option file and output appropriate scores, sequence and scoring

    :param: filename                    : read an option file
    :return: sequence,scores, scoring  :  sequence contains all sequences for alignment
                                          scores contains all match, mismatch scores of alignment
                                          scoring contains gap score and if scores == 0, match, mismatch scores
    """
    sequence, scores, scoring = [], [], {}
    with open(filename, 'r') as fh:
        for line in fh:

            # push the sequences into sequence as a list
            if line[-7:-1] == '.fasta':
                sequence = read_fasta(line.strip("\n"))

            # push the scoring matrix into a nested dictionary to easily find matching
            if line[-5:-1] == ".mat":
                scoring = read_matrix(line.strip("\n"))

            # push all scores to a list, for gapPenalty as well, and in case there is no scoring matrix
            if line != '\n' and line[-7:-1] != '.fasta' and line[-5:-1] != ".mat" and line[0] != "#":
                scores.append(line.strip('\n'))

    # Push the gap extension penalty and affine gap penalties for easier access
    # if they are available
    if len(scores[-1]) != 1:
        list1 = scores[-1].split()
        scores[-1] = list1

    # if there is no scoring matrix then use scores to create a scoring matrix
    if len(scoring) == 0:
        scoring_matrix = {}
        list1 = []
        # push every info about sequence into a distinct list
        for n in sequence[0]:
            for h in sequence[1]:
                if h not in list1:
                    list1.append(h)
                if n not in list1:
                    list1.append(n)
        # create a scoring matrix with different sequence info
        for a in range(len(list1)):
            dict1 = {}
            for b in range(len(list1)):
                if list1[a] == list1[b]:
                    dict1[list1[b]] = int(scores[0])
                else:
                    dict1[list1[b]] = int(scores[1])
            scoring_matrix[list1[a]] = dict1
        scoring = scoring_matrix

    return sequence, scores, scoring


def needleman_wunsch(filename):
    """
    needleman wunsch function takes in a filename(which get processed by optionFile function)
    and align sequences with each other using global alignment technique

    :param: filename : name of option file
    :return: string  : information of sequence alignment
    """

    # read and create variables for use
    optionfile = optionFile(filename)
    sequence1, sequence2 = optionfile[0][1], optionfile[0][0]
    scores = optionfile[1]
    scoring = optionfile[2]
    GapExPenalty = int(scores[2][0])

    # create list and make the list with gapPenalty
    list3 = [1 for x in range(len(sequence2) + 1)]
    for i in range(1, len(list3) + 1):
        list3[i - 1] = i * GapExPenalty

    # create a matrix with scores from the gapPenalty
    matrix = [[1 for a in range(1)] for b in range(len(sequence1) + 1)]
    for i in range(1, len(matrix) + 1):
        matrix[i - 1][0] = i * GapExPenalty
    matrix[0] = list3

    # algorithm to run and aligning sequences
    for i in range(1, len(matrix)):
        nuc1 = sequence1[i - 1]
        for j in range(1, len(matrix[0])):
            nuc2 = sequence2[j - 1]

            # three options to check for possible scores
            opt1 = matrix[i - 1][j - 1]
            opt2 = matrix[i][j - 1]
            opt3 = matrix[i - 1][j]
            """
            find the maximum value from 3 options with:
            option 1 plus the matching state of nuc1 and nuc2
            option 2 and 3 plus the gapPenalty
            """
            # add scores into the matrix
            num = max(opt1 + scoring[nuc1][nuc2], opt2 + GapExPenalty, opt3 + GapExPenalty)
            matrix[i].append(num)

    """
    This is a traceback using the matrix we just created earlier
    This can only find 1 optimal sequence alignment
    The traceback starts at bottom right corner and then move up 
    according to the sequence alignment
    """
    i, j = len(matrix) - 1, len(matrix[0]) - 1
    list1, list2 = [], []

    while i != 0 and j != 0:
        # start at bottom right corner
        start = matrix[i][j]
        # check if it's a match or not
        if sequence1[i - 1] == sequence2[j - 1]:
            list1.append(sequence1[i - 1])
            list2.append(sequence2[j - 1])
            i = i - 1
            j = j - 1
        else:
            # check the mismatch score with current score
            if start == matrix[i - 1][j - 1] + scoring[sequence1[i - 1]][sequence2[j - 1]]:
                list1.append(sequence1[i - 1])
                list2.append(sequence2[j - 1])
                i = i - 1
                j = j - 1
            # check the gap score with starting score
            elif start == matrix[i - 1][j] + GapExPenalty:
                list1.append(sequence1[i - 1])
                list2.append("-")
                i = i - 1
            # check the gap score with starting score
            elif start == matrix[i][j - 1] + GapExPenalty:
                list1.append("-")
                list2.append(sequence2[j - 1])
                j = j - 1

    # put everything remaining into the sequence, if we reach the top
    for a in range(i):
        list1.append(sequence1[a])
        list2.append("-")
    for b in range(j):
        list1.append("-")
        list2.append(sequence2[b])
    score = 0
    for x in range(len(matrix)):
        for y in range(len(matrix[0])):
            if score < matrix[x][y]:
                score = matrix[x][y]

    return "Global Alignment: \n" + sequenceAligningInfo(list1, list2, score, sequence2, sequence1)


def smith_waterman(filename):
    """
    Smith_waterman function takes in a file name(which will be processed by optionFile function)
    and align the sequences with each other using local sequence alignment technique

    :param  :filename
    :return :string (information about sequence alignment)

    """
    optionfile = optionFile(filename)
    sequence1, sequence2 = optionfile[0][1], optionfile[0][0]
    scores = optionfile[1]
    scoring = optionfile[2]
    GapPenalty = int(scores[2][0])

    # create list and make the list with gapPenalty
    list3 = [0 for x in range(len(sequence2) + 1)]

    # create a matrix with scores from the gapPenalty
    matrix = [[0 for a in range(1)] for b in range(len(sequence1) + 1)]
    matrix[0] = list3
    largest_score = 0
    positions = []
    for i in range(1, len(matrix)):
        nuc1 = sequence1[i - 1]
        for j in range(1, len(matrix[0])):
            nuc2 = sequence2[j - 1]
            # three options to check for possbile scores
            opt1 = matrix[i - 1][j - 1]
            opt2 = matrix[i][j - 1]
            opt3 = matrix[i - 1][j]
            """
            find the maximum value from 3 options with:
            option 1 plus the matching state of nuc1 and nuc2
            option 2 and 3 plus the gapPenalty
            """
            # if else statement to see if it's a match or not
            num = max(0, opt1 + scoring[nuc1][nuc2], opt2 + GapPenalty, opt3 + GapPenalty)
            matrix[i].append(num)
            if num > largest_score:
                largest_score = num

    for x in range(1, len(matrix)):
        for y in range(1, len(matrix[0])):
            if matrix[x][y] == largest_score:
                position = (x, y)
                if position not in positions:
                    positions.append(position)
    list10 = []
    position = 0

    # runs all possible locations of largest scores so that we can get all possible sequence alignments
    while position < len(positions):
        position1, position2 = positions[position][0], positions[position][1]
        list1, list2 = [], []
        check = True
        while check:
            # start at bottom right corner
            start = matrix[position1][position2]
            if start == 0:
                check = False
                break
                # check if it's a match or not
            if sequence1[position1 - 1] == sequence2[position2 - 1]:
                list1.append(sequence1[position1 - 1])
                list2.append(sequence2[position2 - 1])
                position1 -= 1
                position2 -= 1
            else:
                # check the mismatch score with current score
                if start == matrix[position1 - 1][position2 - 1] + scoring[sequence1[position1 - 1]][
                    sequence2[position2 - 1]]:
                    list1.append(sequence1[position1 - 1])
                    list2.append(sequence2[position2 - 1])
                    position1 -= 1
                    position2 -= 1
                # check the gap score with starting score
                elif start == matrix[position1 - 1][position2] + GapPenalty:
                    list1.append(sequence1[position2 - 1])
                    list2.append("-")
                    position1 -= 1
                # check the gap score with starting score
                elif start == matrix[position1][position2 - 1] + GapPenalty:
                    list1.append("-")
                    list2.append(sequence2[position2 - 1])
                    position2 -= 1

        # create 2 subsequences to compare them for similarity

        str1 = sequence1[position1:positions[position][0]]
        str2 = sequence2[position2:positions[position][1]]
        pos = position + 1
        str3 = "Local Alignment " + str(pos) + " of " + str(len(positions)) + ":\n" + sequenceAligningInfo(list1, list2,
                                                                                                           largest_score,
                                                                                                           str2, str1)
        position += 1
        list10.append(str3)
        list10.append("\n")

    return list10


def main():
    """
    a main function to call all available functions
    """
    b = needleman_wunsch("test4.in")
    a = smith_waterman("test4.in")
    print(b)
    print(a)


if __name__ == '__main__':
    main()
