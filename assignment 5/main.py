from __future__ import annotations

import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
import itertools, doctest
import operator

# https://realpython.com/working-with-files-in-python/#getting-a-directory-listing
# use this website to get all files inside a directory

# use this website for sorting the dictionary in descending order
# w3resource.com/python-exercises/dictionary/python-data-type-dictionary-exercise-1.php

import os


def brca_mutation_counts():
    """
    This function outputs a dictionary in descending order of mutated genes in breast cancer
    """

    # the files are inside a directory so make a list contain all file names
    files = os.listdir("gdac.broadinstitute.org_BRCA.Mutation_Packager_Calls.Level_3.2016012800.0.0/")

    # a dictionary where key = gene name, and value = number of times it appear
    get_gene = {}

    # a for loop to call each file inside the directory
    for file_name in files:
        file = open("gdac.broadinstitute.org_BRCA.Mutation_Packager_Calls.Level_3.2016012800.0.0/" + file_name, "r")
        file.readline()
        # read each line in the file
        for line in file:
            list = line.split()
            gene_name = list[0]
            mutation = list[8]
            # Ignore mutations that are silents
            if mutation != "Silent":
                if gene_name in get_gene:
                    get_gene[gene_name] += 1
                else:
                    get_gene[gene_name] = 1

                    # sort the dictionary in descending order
    sorted_gene = dict(sorted(get_gene.items(), key=operator.itemgetter(1), reverse=True))
    return sorted_gene


def glioma_survival_curves(glioma="GBMLGG.merged_only_clinical_clin_format.txt"):
    """
  This function takes in a clinical file and output the survival graph of patients with cancer 
  in  brain and central nervous system. 
  There will be 3 lines: all patients, patients with brain cancer, patients with central nervous system 
  cancer
  """
    # read the file
    file = open(glioma, "r")

    # ignore the first line of the file
    file.readline()

    # a lists containing the days since they were dead, or the time they last had a check up
    patient_durations = []

    # status of patients, either dead or alive
    patient_status = []

    # a list containing only patients with brain cancer, when they were dead, or the last time they had a check up
    patient_brain_durations = []

    # status of brain cancer patients, either dead or alive
    patient_brain_status = []

    # a list containing only patients with cancer in central nervous system, when they were dead, or the last time they had a check up
    patient_spinal_durations = []

    # status of patients with cancer in central nervous system, either dead or alive
    patient_spinal_status = []

    # read all the lines inside the file
    lines = file.readlines()

    # get the information on the patients when they were dead
    patients_death = lines[19].split()

    # information on when patients last had a check up
    patients_alive = lines[21].split()

    # status of patient
    vitality = lines[1120].split()

    # where their cancer is
    tumor_tissue = lines[1118].split("\t")

    # run through all the patients
    for patient in range(1, len(patients_death)):
        # ignore those with NA status
        if patients_death[patient] == patients_alive[patient] == "NA" or vitality[patient] == "NA":
            continue
        if (patients_death[patient] != "NA" or  patients_alive[patient] != "NA") and vitality[patient] != "NA":
            # if they were dead, put the time when they last had a check up
            if patients_death[patient] != "NA":
                patient_durations.append(int(patients_death[patient]))
                patient_status.append(1)

            # If they went for a checkup, put the time inside the list
            elif patients_death[patient] == "NA" and patients_alive[patient] != "NA":
                patient_durations.append(int(patients_alive[patient]))

                # check if they died while having a checkup or not
                if vitality[patient] == "alive":
                    patient_status.append(0)
                elif vitality[patient] == "dead":
                    patient_status.append(1)

            # check what caner they had
            if tumor_tissue[patient] == "brain":
                patient_brain_durations.append(patient_durations[-1])
                patient_brain_status.append(patient_status[-1])
            if tumor_tissue[patient] == "central nervous system":
                patient_spinal_durations.append(patient_durations[-1])
                patient_spinal_status.append(patient_status[-1])

    # plot the graph
    kmf = KaplanMeierFitter()
    kmf.fit(patient_durations, patient_status, label="all")
    kmf.plot_survival_function()
    kmf.fit(patient_brain_durations, patient_brain_status, label="brain")
    kmf.plot_survival_function()
    kmf.fit(patient_spinal_durations, patient_spinal_status, label="cns")
    kmf.plot_survival_function()
    plt.ylabel("Survival Rate")
    plt.xlabel("Days")
    plt.show()


def glioma_survival_curves_simplified(glioma="GBMLGG.merged_only_clinical_clin_format.simplified.txt"):
    """
  This function takes in a simplified clinical file and output the survival graph of patients with cancer 
  in  brain and central nervous system. 
  There will be 3 lines: all patients, patients with brain cancer, patients with central nervous system 
  cancer
  """
    # read the file
    file = open(glioma, "r")

    # ignore the first line of the file
    file.readline()

    # a lists containing the days since they were dead, or the time they last had a check up
    patient_durations = []

    # status of patients, either dead or alive
    patient_status = []

    # a list containing only patients with brain cancer, when they were dead, or the last time they had a check up
    patient_brain_durations = []

    # status of brain cancer patients, either dead or alive
    patient_brain_status = []

    # patients with cancer in central nervous system, when they were dead, or the last time they had a check up
    patient_spinal_durations = []

    # status of patients with cancer in central nervous system, either dead or alive
    patient_spinal_status = []

    # run through all the lines inside the file
    for line in file:
        list1 = line.split()
        lines = str(line)
        if "dead" in list1 or "alive" in list1:

            # if they were dead, put the time when they last had a check up
            if list1[10] != "NA":
                patient_durations.append(int(list1[10]))
                patient_status.append(1)

            # If they went for a checkup, put the time inside the list
            elif list1[10] == "NA" and list1[11] != "NA":
                patient_durations.append(int(list1[11]))

                # check if they died while having a checkup or not
                if "dead" in list1:
                    patient_status.append(1)
                elif "alive" in list1:
                    patient_status.append(0)

            # check what caner they had
            if "central nervous system" in lines:
                patient_spinal_durations.append(patient_durations[-1])
                patient_spinal_status.append(patient_status[-1])
            elif "brain" in lines:
                patient_brain_durations.append(patient_durations[-1])
                patient_brain_status.append(patient_status[-1])

    # plot the graph
    kmf = KaplanMeierFitter()
    kmf.fit(patient_durations, patient_status, label="all")
    kmf.plot_survival_function()
    kmf.fit(patient_brain_durations, patient_brain_status, label="brain")
    kmf.plot_survival_function()
    kmf.fit(patient_spinal_durations, patient_spinal_status, label="cns")
    kmf.plot_survival_function()
    plt.ylabel("Survival Rate")
    plt.xlabel("Days")
    plt.show()


def brca_survival_curves_simplified(brca="BRCA.merged_only_clinical_clin_format.simplified.txt"):
    """
  this function read a simplified clincal file and output the survival graph of breast cancer patients of 
  all stages.
  There are 5 lines in the graph: all patients, patients with stage 1, patients with stage 2, patients with   stage 3, and patients with stage 4
  """

    # read the file
    file = open(brca, "r")

    # ignore the first line
    file.readline()

    # a lists containing the days since they were dead, or the time they last had a check up
    patient_durations = []

    # status of patients, either dead or alive
    patient_status = []

    # patients with stage 1 breast cancer
    patient_stage1_durations = []
    patient_stage1_status = []

    # patients with stage 2 breast cancer
    patient_stage2_durations = []
    patient_stage2_status = []

    # patient with stage 3 breast cancer
    patient_stage3_durations = []
    patient_stage3_status = []

    # patient with stage 4 breast cancer
    patient_stage4_durations = []
    patient_stage4_status = []

    for line in file:
        list1 = line.split()
        for i in range(len(list1)):
            if list1[i].startswith("-") and list1[i][1:].isdigit():
                break

        # get the time when they died or last had check up
        if i == len(list1) - 1:
            a = list1.index("hispanic")
            death = list1[a - 6]
            alive = list1[a - 5]
        else:
            death = list1[i + 1]
            alive = list1[i + 2]

        if death != "NA":
            patient_durations.append(int(death))
            patient_status.append(1)
        elif death == "NA" and alive != "NA":
            patient_durations.append(int(alive))
            if list1[-3] == "alive":
                patient_status.append(0)
            else:
                patient_status.append(1)

        if "stage" in list1:
            b = list1.index("stage")
            if list1[b + 1][:2] == "iv":
                patient_stage4_durations.append(patient_durations[-1])
                patient_stage4_status.append(patient_status[-1])
            elif list1[b + 1][:3] == "iii":
                patient_stage3_durations.append(patient_durations[-1])
                patient_stage3_status.append(patient_status[-1])
            elif list1[b + 1][:2] == "ii":
                patient_stage2_durations.append(patient_durations[-1])
                patient_stage2_status.append(patient_status[-1])
            elif list1[b + 1][:1] == "i":
                patient_stage1_durations.append(patient_durations[-1])
                patient_stage1_status.append(patient_status[-1])

    # plot the graph
    kmf = KaplanMeierFitter()
    kmf.fit(patient_durations, patient_status, label="all")
    kmf.plot_survival_function()
    kmf.fit(patient_stage1_durations, patient_stage1_status, label="stage1")
    kmf.plot_survival_function()
    kmf.fit(patient_stage2_durations, patient_stage2_status, label="stage2")
    kmf.plot_survival_function()
    kmf.fit(patient_stage3_durations, patient_stage3_status, label="stage3")
    kmf.plot_survival_function()
    kmf.fit(patient_stage4_durations, patient_stage4_status, label="stage4")
    kmf.plot_survival_function()
    plt.ylabel("Survival Rate")
    plt.xlabel("Days")
    plt.show()


def brca_survival_curves(brca="BRCA.merged_only_clinical_clin_format.txt"):
    """
  This function reads a clinical file and output a survival graph of breast cancer of all stages. 
  There will be 5 lines: all the patients, patient of stage 1, patients of stage 2, patients of 
  stage 3, patients of stage 4
  """
    file = open(brca, "r")

    #ignore the first line
    file.readline()

    # a lists containing the days since they were dead, or the time they last had a check up
    patient_durations = []

    # status of patients, either dead or alive
    patient_status = []

    # patients with stage 1 breast cancer
    patient_stage1_durations = []
    patient_stage1_status = []

    # patients with stage 2 breast cancer
    patient_stage2_durations = []
    patient_stage2_status = []

    # patient with stage 3 breast cancer
    patient_stage3_durations = []
    patient_stage3_status = []

    # patient with stage 4 breast cancer
    patient_stage4_durations = []
    patient_stage4_status = []

    #read all lines in the file
    lines = file.readlines()

    #get all patients time of death
    patients_death = lines[34].split()

    #get all patients time of last check up
    patients_alive = lines[36].split()

    #get patients status: death or alive
    vitality = lines[1477].split()

    #get patients stages at breast cancer
    stages = lines[1460].split()

    stage = 2
    patient = 1
    while stage < len(stages):
        if patients_death[patient] == patients_alive[patient] == "NA" or stages[stage] == "x" or vitality[
            patient] == "NA":
            patient += 1
            stage += 2
        else:
            if patients_death[patient] != "NA":
                patient_durations.append(int(patients_death[patient]))
                patient_status.append(1)
            elif patients_death[patient] == "NA" and patients_alive[patient] != "NA":
                patient_durations.append(int(patients_alive[patient]))
                if vitality[patient] == "alive":
                    patient_status.append(0)
                elif vitality[patient] == "dead":
                    patient_status.append(1)
            if stages[stage][:2] == "iv":
                patient_stage4_durations.append(patient_durations[-1])
                patient_stage4_status.append(patient_status[-1])
            elif stages[stage][:3] == "iii":
                patient_stage3_durations.append(patient_durations[-1])
                patient_stage3_status.append(patient_status[-1])
            elif stages[stage][:2] == "ii":
                patient_stage2_durations.append(patient_durations[-1])
                patient_stage2_status.append(patient_status[-1])
            elif stages[stage][:1] == "i":
                patient_stage1_durations.append(patient_durations[-1])
                patient_stage1_status.append(patient_status[-1])
        stage += 2
        patient += 1

    # plot the graph
    kmf = KaplanMeierFitter()
    kmf.fit(patient_durations, patient_status, label="all")
    kmf.plot_survival_function()
    kmf.fit(patient_stage1_durations, patient_stage1_status, label="stage1")
    kmf.plot_survival_function()
    kmf.fit(patient_stage2_durations, patient_stage2_status, label="stage2")
    kmf.plot_survival_function()
    kmf.fit(patient_stage3_durations, patient_stage3_status, label="stage3")
    kmf.plot_survival_function()
    kmf.fit(patient_stage4_durations, patient_stage4_status, label="stage4")
    kmf.plot_survival_function()
    plt.ylabel("Survival Rate")
    plt.xlabel("Days")

    plt.show()


# brca_survival_curves_simplified()
#glioma_survival_curves()

# glioma_survival_curves_simplified()
# glioma_survival_curves()
#brca_survival_curves_simplified()
brca_survival_curves()
