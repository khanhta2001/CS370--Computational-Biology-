"""
This function read the file mi.owl and return a dictionary of parent/child relations
"""
def relations(filename = "mi.owl"):
    """
    read the file mi.owl and return a dictionary of parent/child relations

    :param filename: call file mi.owl
    :return:
    """

    # open and read file and store it into a new variable named lines
    file = open(filename, "r")
    lines = file.readlines()

    # get the length of lines for for-loop
    length = len(lines)

    # initial a dictionary
    # the key of the dictionary will be the children
    # the value of each key will be the parent(s)
    parentchild = {}

    # for-loop to create the dictionary
    for i in range(length):

        # a list to put multiple parent values
        values = []

        # if there is "id:" in the line[i]
        if lines[i].find("id:") != -1:

            # the key of the dictionary will be the strings in the position from 4 to 11
            key = lines[i][4:11]

            # initial an empty parent value for the key
            # if the child value (key) does not have a parent (value of each key), there will be an empty list for that child (key) in the dictionary
            parentchild[key] = []

        # else if there is "is_a" in the line[i]
        elif lines[i].find("is_a") != -1:

            # the string in the postion from 6 to 13 will be the first value
            value_1 = lines[i][6:13]

            # append the first value into the values list
            values.append(value_1)

            # check if there is multiple values
            # append the string in the position from 6 to 13 in the line[i-1] to the value list as well
            if lines[i - 1].find("is_a") != -1:
                value_2 = lines[i - 1][6:13]
                values.append(value_2)

            # assign the values list to the key
            parentchild[key] = values

    # return the dictionary
    return parentchild


"""
This function takes a dictionary of protein interaction as input
This functions return a set of MI codes that represents physical protein-protein interactions
"""
def get_physical_interaction_codes(data_relations = relations()):
    """
    This function takes a dictionary of protein interaction as input
    This functions return a set of MI codes that represents physical protein-protein interactions

    :param data_relations: call in function relations()
    :return:
    """

    # Start at MI:0407
    to_look_at = ["MI:0407"]

    # initial a set to return
    looked = set()

    # initial an iterator
    i = 0

    # add MI:0407 to the set
    looked.update(["MI:0407"])

    # start a while-loop
    # there are 71 MI:0407 descendant codes, which belongs to physical interactions
    while True:
        # start with the first element of to_look_at
        start_code = to_look_at[0]
        # start the loop with the key in the input dictionary
        for key in data_relations:
            # another loop with the value in the key
            for value in data_relations[key]:
                # check if the value equals start_code
                # if the value equals to the start_code, means that the key is a physical interaction
                # append the value to to_look_at for the next round of loop
                if value == start_code:
                    to_look_at.append(key)
                    looked.update([key])
        # after finding all the values connecting with the key, pop that key out of the set
        to_look_at.pop(0)
        if len(to_look_at) == 0:
            break

        # update iterator
        i += 1

    # return the set
    return looked


"""
This function takes a dictionary of protein interaction as input
This functions return a set of MI codes that represents experimentally detected interactions
"""


def get_experimental_codes(data_relations = relations()):
    # Start at MI:0045
    to_look_at = ["MI:0045"]
    # initial a set
    looked = set()
    # initial an iterator
    i = 0
    # add MI:04045 to the set
    looked.update(["MI:0045"])
    # start a while-loop
    # there are 311 MI:0045 descendant codes, which belongs to experimentally detected interactions
    while True:
        # start with the first element of to_look_at
        start_code = to_look_at[0]
        # a nested for-loop to loop through the key and the value in the data_relations dictionary
        for key in data_relations:
            for value in data_relations[key]:
                # check if the value equals start_code
                # if the value equals to the start_code, means that the key is experimentally detected interactions
                # append the value to to_look_at for the next round of loop
                if value == start_code:
                    to_look_at.append(key)
                    looked.update([key])
        # after finding all the values connecting with the key, pop that key out of the set
        to_look_at.pop(0)
        if len(to_look_at) == 0:
            break

        # update iterator
        i += 1

    # return the set
    return looked


def process_biogrid_file(biogrid="BIOGRID-ORGANISM-Homo_sapiens-4.4.207.mitab.txt",
                         networkfile="H.sapiens.biogrid.4.4.207.network"):
    # get physical interactions and experimental interactions to check human genes
    physical_interactions = get_physical_interaction_codes()
    experimental_interactions = get_experimental_codes()

    # open the file
    file = open(biogrid, "r")
    file.readline()

    # output all the interaction-interaction in the BIOGRID file with output[key] = list()
    output = {}

    # total interactions in the BIOGRID file
    total_interactions_processed = 0

    # total interactions in the network file
    total_interactions_network = 0

    # total non-human interactions in the BIOGRID file
    Interspecies = 0

    # total physical interactions discarded in the BIOGRID file
    physical_interaction = 0

    # total experimental interactions discarded in the BIOGRID file
    experimental_interaction = 0

    # total of how many interactions that repeats itself
    genes_loop = 0

    # totla of how many interactions that are not duplicate of other interaction
    genes_dup = 0

    # a loop to read each line in BIOGRID file
    for lines in file:
        # calculate the total interactions in the BIOGRID file
        total_interactions_processed += 1

        # split the line into a list for easier finding of the protein
        lists = lines.split("biogrid")

        # get the proteins in the line
        protein_1 = lists[1][lists[1].find('locuslink:') + 10:lists[1].find('|', lists[1].find('locuslink:'))]
        protein_2 = lists[2][lists[2].find('locuslink:') + 10:lists[2].find('|', lists[2].find('locuslink:'))]
        protein_1 = protein_1.split()[0]
        protein_2 = protein_2.split()[0]

        # get experimental interaction
        interactions_experimental = lists[2][lists[2].find("psi-mi:") + 8:lists[2].find("psi-mi:") + 15]

        # get physical interaction
        indexofPSI = lists[2].find("psi-mi:", lists[2].find("psi-mi:") + 15)
        interactions_physical = lists[2][indexofPSI+ 8:indexofPSI + 15]

        # get the first taxid
        taxid1 = lists[2][lists[2].find("taxid:") + 6:lists[2].find("taxid:") + 10]
        
        # get the second taxid
        indexofT2 = lists[2].find("taxid:", lists[2].find("taxid:") + 10)
        taxid2 = lists[2][indexofT2 + 6:indexofT2 + 10]

        # calculate the total non-human interactions
        if taxid1 != "9606" or taxid2 != "9606":
            Interspecies += 1

        # calculate the experimental interactions discarded
        if interactions_experimental not in experimental_interactions and taxid1 == taxid2 == "9606":
            experimental_interaction += 1

        # calculate physical interactions discarded
        if interactions_physical not in physical_interactions and taxid1 == taxid2 == "9606":
            physical_interaction += 1

        # calculate how many interactions that repeats itself
        if protein_1 == protein_2 and taxid1 == taxid2 == "9606" and interactions_experimental in experimental_interactions and interactions_physical in physical_interactions:
            genes_loop += 1

        # check if protein_1 is alphabetically smaller than protein_2
        if protein_1 > protein_2:
            med = protein_1
            protein_1 = protein_2
            protein_2 = med

        # how many interactions that are not duplicate of other interaction

        # check if it is protein-protein interaction
        # if it is, put it into the output file
        if taxid1 == taxid2 == "9606" and protein_1 != protein_2 and interactions_experimental in experimental_interactions and interactions_physical in physical_interactions:
            # check if protein_1 and protein_2 in the output file
            if protein_1 not in output:
                if protein_2 not in output:
                # if not in the output then create a key-value in output like protein_1 = protein_2
                    output[protein_1] = [protein_2]
                    # add 1 to the interaction in the network
                    total_interactions_network += 1
                else:
                    output[protein_1] = [protein_2]
                    # add 1 to the interaction in the network
                    total_interactions_network += 1
            else:
                if protein_2 not in output:
                # if protein_1 in the output then check if protein_2 already in the output[protein_1]
                    if protein_2 not in output[protein_1]:
                        # if protein_2 not in the output then put it in the value list
                        output[protein_1].append(protein_2)
                        # add 1 to the interaction in the network
                        total_interactions_network += 1
                    else:
                        genes_dup += 1
                else:
                    # if protein_1 and protein_2 already in the output then check if protein_2 already in the output[protein_1]
                    if protein_2 not in output[protein_1]:
                        # if protein_2 not in the output then put it in the value list
                        output[protein_1].append(protein_2)
                        # add 1 to the interaction in the network
                        total_interactions_network += 1
                    else:
                        genes_dup += 1

    # open the network file to put the protein-protein interaction inside
    network = open(networkfile, 'w')

    # get the total proteins inside the network
    list1 = sorted(list(output.keys()))
    values = sorted(list(output.values()))
    list2 = []
    for value in values:
        if value not in list2:
            list2.append(value)
    for key in list1:
        if key not in list2:
            list2.append(key)

    # Processing statististics at the top of the network file, starting with #
    network.write("#Total interactions processed: " + str(total_interactions_processed) + "\n")
    network.write("#Total interactions in network: " + str(total_interactions_network) + "  " + str(
        100 * (total_interactions_network / total_interactions_processed))[:4] + "%" + "\n")
    network.write("#Total proteins in network: " + str(len(list1)) + "\n")
    network.write("#Interspecies interactions discarded: " + str(Interspecies) + "  " + str(
        100 * (Interspecies / total_interactions_processed))[:4] + "%" + "\n")
    network.write("#Non-physical interactions discarded: " + str(physical_interaction) + "  " + str(
        100 * (physical_interaction / total_interactions_processed))[:4] + "%" + "\n")
    network.write("#Predicted interactions discarded: " + str(experimental_interaction) + "  " + str(
        100 * (experimental_interaction / total_interactions_processed))[:4] + "%" + "\n")
    network.write(
        "#Self-loops discarded: " + str(genes_loop) + "  " + str(100 * (genes_loop / total_interactions_processed))[
                                                             :4] + "%" + "\n")
    network.write("#Duplicate interactions discarded: " + str(genes_dup) + "  " + str(
        100 * (genes_dup / total_interactions_processed))[:4] + "%" + "\n")

    # add the protein-protein in the alphabetical order into network file
    keys = sorted(list(output.keys()))
    for key in keys:
        output[key].sort()
        for value in output[key]:
            printout = key + " " + value + "\n"
            network.write(printout)


process_biogrid_file()
