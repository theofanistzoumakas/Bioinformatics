import random
import copy
import numpy as np
import math
from collections import defaultdict
def production():
    pattern1 = 'AATTGA'
    pattern2 = 'CGCTTAT'
    pattern3 = 'GGACTCAT'
    pattern4 = 'TTATTCGTA'

    patterns = [pattern1, pattern2, pattern3, pattern4]

    alphabet = ['A', 'C', 'G', 'T']
    characters_for_replace = ['A', 'C', 'G', 'T', '']

    sequence = ''
    random_times = random.randint(1, 3)

    letter_list = random.sample(alphabet, k=random_times)

    for letter in letter_list:
        sequence = sequence + letter

    for pattern in patterns:

        random_times = random.randint(1, 2)
        chosen_random_position = -1

        for times in range(random_times):
            random_position = random.randint(0, len(pattern) - 1)
            while random_position == chosen_random_position:
                random_position = random.randint(0, len(pattern) - 1)
            chosen_random_position = random_position

            choice = random.randint(0, 1)
            if choice == 0:
                pattern = pattern[:random_position] + '' + pattern[(random_position + 1):len(pattern)]
            elif choice == 1:
                random_letter = random.choice(alphabet)
                pattern = pattern[:random_position] + random_letter + pattern[(random_position + 1):len(pattern)]

        sequence = sequence + pattern

    random_times = random.randint(1, 2)

    letter_list = random.sample(alphabet, k=random_times)

    for letter in letter_list:
        sequence = sequence + letter

    return sequence


def consensus(first_seq, second_seq):
    consensus_seq = ''
    list1 = list(first_seq)
    list2 = list(second_seq)

    for (i, j) in zip(list1, list2):
        if i == j:
            consensus_seq += i
        elif i == '_' and j != '_':
            consensus_seq += j
        elif j == '_' and i != '_':
            consensus_seq += i
        else:
            consensus_seq += random.choice([i, j])
    return consensus_seq


def find_path(arr, arr2_with_values, sequence1, sequence2):
    seq1 = ''
    seq2 = ''
    i = len(sequence2)
    j = len(sequence1)
    while (arr2_with_values[i][j] != None):
        if arr2_with_values[i][j] == [-1, -1]:
            seq1 = sequence1[-1] + seq1
            sequence1 = sequence1[:-1]
            seq2 = sequence2[-1] + seq2
            sequence2 = sequence2[:-1]
        elif arr2_with_values[i][j] == [0, -1]:
            seq1 = '_' + seq1
            seq2 = sequence2[-1] + seq2
            sequence2 = sequence2[:-1]
        else:
            seq1 = sequence1[-1] + seq1
            sequence1 = sequence1[:-1]
            seq2 = '_' + seq2

        [x, y] = arr2_with_values[i][j]
        i = i + y
        j = j + x
        if arr2_with_values[i][j] == None:
            if i == 0 and j != 0:
                seq1 = sequence1 + seq1
                seq2 = '_' * len(sequence1) + seq2
            elif j == 0 and i != 0:
                seq2 = sequence2 + seq2
                seq1 = '_' * len(sequence2) + seq1

    penalty = 0
    for i in seq1:
        if i == '_':
            penalty += 1
    for i in seq2:
        if i == '_':
            penalty += 1
    for (i, j) in zip(seq1, seq2):
        if i != j:
            penalty += 1

    return penalty, seq1, seq2


def array_constructor(sequence1, sequence2):
    arr = []
    arr2 = []
    number_row = -2
    for i in range(len(sequence2) + 1):
        row = []
        row2 = []
        number_column = 0
        for j in range(len(sequence1) + 1):
            if i == 0:
                row.append(number_column)
                number_column = number_column - 2
                row2.append(None)
            else:
                if j == 0:
                    row.append(number_row)
                    number_row = number_row - 2
                    row2.append(None)
                else:
                    row.append(None)
                    row2.append(None)

        arr.append(row)
        arr2.append(row2)

    for i in range(1, len(sequence2) + 1):
        for j in range(1, len(sequence1) + 1):
            a = 0
            if sequence1[j - 1] == sequence2[i - 1]:
                a = arr[i - 1][j - 1] + 1
            else:
                a = arr[i - 1][j - 1] - 1
            b = arr[i - 1][j] - 2
            c = arr[i][j - 1] - 2
            max_value = max(a, b, c)
            arr[i][j] = max_value
            if max_value == a:
                arr2[i][j] = [-1, -1]
            elif max_value == b:
                arr2[i][j] = [0, -1]
            else:
                arr2[i][j] = [-1, 0]

    penalty, seq1, seq2 = find_path(arr, arr2, sequence1, sequence2)
    return penalty, seq1, seq2


def re_align_gaps(seq_for_re_align, original_seq):
    new_seq = []
    original_index = 0
    for char in original_seq:
        if char == '_':
            new_seq.append('_')
        else:

            new_seq.append(seq_for_re_align[original_index])
            original_index += 1
    return (str("".join(new_seq)))


def viterbi_for_hmm(sequence, indexes, conditions, emission_table):
    d = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    main_prob = [{}]
    p = {}



    for counter_v in range(len(indexes)):
        main_prob[0][indexes[counter_v][0]] = conditions[0][counter_v]
        p[indexes[counter_v][0]] = [indexes[counter_v][0]]

    for t in range(1, len(sequence)):

        updated_path = {}
        main_prob.append({})

        for counter_v2 in range(len(indexes)):
            (probability,symbol) = -100, "MF"
            for y0 in range(len(indexes)):
                (prob_1, state_1) = (main_prob[t - 1][indexes[y0][0]] * conditions[y0+1][counter_v2] , indexes[y0][0])
                if (prob_1 >= probability):


                    if (emission_table[counter_v2][d[sequence[t]]] == 0):
                        (probability, symbol) = (0.00000001 * prob_1, state_1)
                    else:
                        (probability, symbol) = (emission_table[counter_v2][d[sequence[t]]] * prob_1, state_1)







            main_prob[t][indexes[counter_v2][0]] = probability
            updated_path[indexes[counter_v2][0]] = p[symbol] + [indexes[counter_v2][0]]

        p = updated_path

    (probability, symbol) = max([(main_prob[-1][indexes[counter_v3][0]], indexes[counter_v3][0]) for counter_v3 in range(len(indexes))])

    return (probability, p[symbol])






list_of_sequences = []

for i in range(50):
    list_of_sequences.append(production())

datasetA = []
datasetB = []

for i in range(15):
    pop_number = random.randint(0, len(list_of_sequences) - 1)
    number = list_of_sequences.pop(pop_number)
    datasetA.append(number)

datasetB = list(list_of_sequences)

penalties_arr = []

for i in range(len(datasetA)):
    penalties_arr_row = []
    for j in range(len(datasetA)):
        penalties_arr_row.append(None)
    penalties_arr.append(penalties_arr_row)

k = 1
for i in range(len(datasetA)):
    for j in range(k, len(datasetA)):
        penalties_arr[i][j], result_seq1, result_seq2 = array_constructor(datasetA[i], datasetA[j])

    k += 1

supp_datasetA = []
for i in datasetA:
    supp_datasetA.append(i)

counter_seq = 1
counter_conseq = 1
identification = []

while (len(supp_datasetA) > 1):

    # Μπορεί να ενεργοποιηθεί αν βγει από τα σχόλια η προβολή του πίνακα Penalty
    """
    for i in penalties_arr:
        print(i)
    """

    min_value = float('inf')
    row_min = float('inf')
    position_min = float('inf')
    for row in penalties_arr:
        for value in row:
            if value != None:
                if value < min_value:
                    min_value = value
                    row_min = penalties_arr.index(row)
                    position_min = row.index(value)

    first_sequence = supp_datasetA[row_min]
    second_sequence = supp_datasetA[position_min]

    penalties, result_seq1, result_seq2 = array_constructor(first_sequence, second_sequence)

    result = consensus(result_seq1, result_seq2)

    if first_sequence in datasetA:
        identification.append([result_seq1, "seq" + str(counter_seq), "conseq" + str(counter_conseq)])
        counter_seq += 1

    if second_sequence in datasetA:
        identification.append([result_seq2, "seq" + str(counter_seq), "conseq" + str(counter_conseq)])
        counter_seq += 1

    indexer = 0
    for [seq, name, apogon] in identification:
        if (first_sequence == seq):
            identification[indexer][0] = result_seq1
            identification[indexer][2] = "conseq" + str(counter_conseq)
        if (second_sequence == seq):
            identification[indexer][0] = result_seq2
            identification[indexer][2] = "conseq" + str(counter_conseq)
        indexer += 1

    identification.append([result, "conseq" + str(counter_conseq), None])

    counter_conseq += 1

    supp_datasetA.append(result)
    supp_datasetA.remove(first_sequence)
    supp_datasetA.remove(second_sequence)

    new_penalties_arr = []
    for i in range(len(supp_datasetA)):
        penalties_arr_row = []
        for j in range(len(supp_datasetA)):
            penalties_arr_row.append(None)
        new_penalties_arr.append(penalties_arr_row)

    k = 1
    for i in range(len(supp_datasetA)):
        for j in range(k, len(supp_datasetA)):
            new_penalties_arr[i][j], result_seq1, result_seq2 = array_constructor(supp_datasetA[i], supp_datasetA[j])
        k += 1

    penalties_arr = copy.deepcopy(new_penalties_arr)

print("\n\n Οι ακολουθίες πριν την αντιστοίχιση: \n")
for i in identification:
    print(i)

counter_conseq = counter_conseq - 1

for counter in range(counter_conseq, 0, -1):
    the_conseq = None
    indexer = 0
    for [seq, name, apogon] in identification:
        if name == ("conseq" + str(counter)):
            the_conseq = seq
            break

    for [seq, name, apogon] in identification:
        if (apogon == ("conseq" + str(counter))) and (len(the_conseq) > len(seq)):
            variable = re_align_gaps(seq, the_conseq)
            identification[indexer][0] = variable
        indexer += 1

print("\n\n\n Οι πλήρως στοιχισμένες αρχικές ακολουθίες είναι: \n")
total_aligned_sequences = []
for [seq, name, apogon] in identification:
    if name[0:3] == "seq":
        total_aligned_sequences.append(seq)
        print("[", seq, ", ", name, ", ", apogon, "]")

print("\n \n")

f = open("aligned_sequences.txt", "w")
for s in total_aligned_sequences:
    f.write(s)
    f.write("\n")
f.close()

threshold = int(len(total_aligned_sequences) / 2)

for i in range(len(total_aligned_sequences)):
    total_aligned_sequences[i] = list(total_aligned_sequences[i])
array_with_aligned_sequences = np.array(total_aligned_sequences)

rows, columns = array_with_aligned_sequences.shape
conditions =[]

indexing = []
indexing_emission = []
auxilliaryI = 0
counterM = 0
counterI = 0
counterD = 0
for i in range(columns):
    if np.count_nonzero(array_with_aligned_sequences[:,i]=='_')>threshold:
        auxilliaryI += 1
    else:
        counterM +=1
        indexing.append(["M"+str(counterM),i])
        if auxilliaryI>0:
            counterI+=1
            indexing.append(["I" + str(counterI),i-auxilliaryI,auxilliaryI])
            auxilliaryI = 0
        if '_' in array_with_aligned_sequences[:,i]:
                counterD +=1
                indexing.append(["D" + str(counterM),i])

if auxilliaryI > 0:
    counterI += 1
    indexing.append(["I" + str(counterI),len(total_aligned_sequences[0])-auxilliaryI,auxilliaryI])
    auxilliaryI = 0


for i in range(counterM+counterI+counterD+1):
    conditions_row = []
    for j in range(counterM+counterI+counterD):
        conditions_row.append(0)
    conditions.append(conditions_row)


counterm=0
for i in range(len(indexing)-1):
    if indexing[i][0][0] == 'M':
        m_number = indexing[i][0][1:]
        position_i=indexing[i][1]
        position_j = None
        next_state = None
        variable = int(indexing[i][0][1:])+1
        for j in indexing:
            if j[0]=='M'+str(variable):
                position_j = j[1]
                next_state = indexing.index(j)
                break
        if position_j!=None:
            if position_j-position_i==1:
                counterm+=1
                counter_letter_to_letter = 0
                counter_letter_to_gap = 0
                counter_gap_to_letter = 0
                counter_gap_to_gap = 0
                for row in range(rows):
                    if array_with_aligned_sequences[row][position_i]!='_' and array_with_aligned_sequences[row][position_j]!='_':
                        counter_letter_to_letter+=1
                    elif array_with_aligned_sequences[row][position_i] != '_' and array_with_aligned_sequences[row][position_j] == '_':
                        counter_letter_to_gap +=1
                    elif array_with_aligned_sequences[row][position_i] == '_' and array_with_aligned_sequences[row][position_j] != '_':
                        counter_gap_to_letter +=1
                    else:
                        counter_gap_to_gap +=1
                counter_from_letter = counter_letter_to_gap+counter_letter_to_letter
                if counter_from_letter!=0:
                    probability_m_to_m = counter_letter_to_letter/counter_from_letter
                    probability_m_to_d = counter_letter_to_gap / counter_from_letter
                else:
                    probability_m_to_m = 0
                    probability_m_to_d = 0


                counter_from_gap = counter_gap_to_gap + counter_gap_to_letter

                if counter_from_gap!=0:
                    probability_d_to_d = counter_gap_to_gap / counter_from_gap
                    probability_d_to_m = counter_gap_to_letter / counter_from_gap
                else:
                    probability_d_to_d=0
                    probability_d_to_m=0


                if next_state!=None:
                    conditions[i+1][next_state] = round(probability_m_to_m,2)

                previous_state_d = None
                for t in indexing:
                    if t[0] == 'D' + str(m_number):
                        previous_state_d = indexing.index(t)
                        break
                if previous_state_d!=None:
                    conditions[previous_state_d+1][next_state] = round(probability_d_to_m,2)

                next_state_gap = None

                for k in indexing:
                    if k[0] == 'D' + str(variable):
                        next_state_gap = indexing.index(k)
                        break

                if next_state_gap != None:
                    conditions[i + 1][next_state_gap] = round(probability_m_to_d,2)

                if previous_state_d != None and next_state_gap != None:
                    conditions[previous_state_d + 1][next_state_gap] = round(probability_d_to_d,2)
            else:
                i_position = position_i + 1
                i_length = 0
                i_position_in_indexing = 0


                for j in indexing:
                    if j[1] == i_position:
                        i_length = j[2]
                        i_position_in_indexing = indexing.index(j)
                        break

                # Εκκίνηση υπολογισμού M και D πριν το I
                empty_rows = 0
                empty_rows_2 = 0
                countering = 0
                countering_2 = 0
                gap_to_gap_counter_2 = 0

                letter_to_gap_counter_2 = 0
                for row in range(rows):
                    if array_with_aligned_sequences[row][position_i] != '_':
                        countering += 1 #Συνολικό πλήθος που το m δεν έχει κενό
                        is_with_underscore = True
                        for changer in range(i_length):
                            if array_with_aligned_sequences[row][i_position + changer] != '_':
                                is_with_underscore = False
                        if is_with_underscore == True:
                            empty_rows += 1 #Πλήθος απευθείας περάσματος απέναντι
                            if array_with_aligned_sequences[row][position_j] == '_':
                                letter_to_gap_counter_2 += 1 #Πλήθος απευθείας περασμάτων με κενό στο άλλο m
                    else:
                        countering_2 += 1 #Συνολικό πλήθος που το m έχει κενό
                        is_with_underscore_2 = True
                        for changer_2 in range(i_length):
                            if array_with_aligned_sequences[row][i_position + changer_2] != '_':
                                is_with_underscore_2 = False
                                #gap_to_gap_counter_2 += 1
                        if is_with_underscore_2 == True:
                            empty_rows_2 += 1
                            if array_with_aligned_sequences[row][position_j] == '_':
                                gap_to_gap_counter_2 += 1

                probability_m_to_d = letter_to_gap_counter_2/countering
                probability_m_to_m = (empty_rows-letter_to_gap_counter_2)/countering
                probability_m_to_i = (countering-empty_rows)/countering

                probability_d_to_d = 0
                probability_d_to_m = 0
                probability_d_to_i = 0
                if (countering_2 > 0):
                    probability_d_to_i = (countering_2-empty_rows_2)/countering_2
                    probability_d_to_m = (empty_rows_2-gap_to_gap_counter_2)/countering_2
                    probability_d_to_d = gap_to_gap_counter_2/countering_2

                previous_state_d = None
                for t in indexing:
                    if t[0] == 'D' + str(m_number):
                        previous_state_d = indexing.index(t)
                        break

                next_state_gap = None
                for k in indexing:
                    if k[0] == 'D' + str(variable):
                        next_state_gap = indexing.index(k)
                        break


                conditions[i+1][next_state] = round(probability_m_to_m,2)

                conditions[i+1][i_position_in_indexing] = round(probability_m_to_i,2)

                if (next_state_gap != None):
                    conditions[i+1][next_state_gap] = round(probability_m_to_d,2)

                if (previous_state_d != None):
                    conditions[previous_state_d + 1][i_position_in_indexing] = round(probability_d_to_i,2)

                if (next_state_gap != None and previous_state_d!= None):
                    conditions[previous_state_d + 1][next_state_gap] = round(probability_d_to_d, 2)

                if (next_state != None and previous_state_d != None):
                    conditions[previous_state_d + 1][next_state] = round(probability_d_to_m, 2)

                # Ολοκλήρωση Υπολογισμού M και D πριν το I

                # Έναρξη υπολογισμού του  I

                valuable_i_rows = 0
                total_counter = 0
                counter_to_back = 0
                i_to_d_counter = 0
                for row in range(rows):
                    is_empty = True
                    for changers in range(i_length):
                        if array_with_aligned_sequences[row][i_position + changers] != '_':
                            is_empty = False
                    if is_empty == False:
                        valuable_i_rows += 1
                        local_counter = 0
                        for changerss in range(i_length):
                            if array_with_aligned_sequences[row][i_position + changerss] != '_':
                                local_counter += 1
                        counter_to_back += local_counter - 1
                        total_counter += local_counter

                        if array_with_aligned_sequences[row][position_j] == '_':
                            i_to_d_counter += 1

                probability_i_to_i = counter_to_back/total_counter
                prabability_i_to_d = i_to_d_counter/total_counter
                probability_i_to_m = (valuable_i_rows - i_to_d_counter) /total_counter

                next_state_gap = None
                for k in indexing:
                    if k[0] == 'D' + str(variable):
                        next_state_gap = indexing.index(k)
                        break


                conditions[i_position_in_indexing + 1][i_position_in_indexing] = round(probability_i_to_i,2)

                if (next_state != None):
                    conditions[i_position_in_indexing + 1][next_state] = round(probability_i_to_m,2)

                if (next_state_gap != None):
                    conditions[i_position_in_indexing + 1][next_state_gap] = round(prabability_i_to_d, 2)

                # Ολοκλήρωση υπολογισμού του  I

for i in range(len(indexing)-1):
    if indexing[i][0][0] == 'I' and indexing[i][1] == 0:

            i_position = 0
            i_length = indexing[i][2]
            i_position_in_indexing = i


            position_j = 0
            state_d = None
            for t in indexing:
                if t[0] == 'D1':
                    state_d = indexing.index(t)
                    break

            state_m = None
            for t in indexing:
                if t[0] == 'M1':
                    state_m = indexing.index(t)
                    position_j = t[1]
                    break

            valuable_i_rows = 0
            total_counter = 0
            counter_to_back = 0
            i_to_d_counter = 0
            for row in range(rows):
                is_empty = True
                for changers in range(i_length):
                    if array_with_aligned_sequences[row][i_position + changers] != '_':
                        is_empty = False
                if is_empty == False:
                    valuable_i_rows += 1
                    local_counter = 0
                    for changerss in range(i_length):
                        if array_with_aligned_sequences[row][i_position + changerss] != '_':
                            local_counter += 1
                    counter_to_back += local_counter - 1
                    total_counter += local_counter

                    if array_with_aligned_sequences[row][position_j] == '_':
                        i_to_d_counter += 1

            if (total_counter > 0):
                probability_i_to_i = counter_to_back / total_counter
                prabability_i_to_d = i_to_d_counter / total_counter
                probability_i_to_m = (valuable_i_rows - i_to_d_counter) / total_counter



            conditions[i_position_in_indexing + 1][i_position_in_indexing] = round(probability_i_to_i, 2)

            if (state_m != None):
                conditions[i_position_in_indexing + 1][state_m] = round(probability_i_to_m, 2)

            if (state_d != None):
                conditions[i_position_in_indexing + 1][state_d] = round(prabability_i_to_d, 2)



for i in indexing:
    timesA = 0
    timesC = 0
    timesG = 0
    timesT = 0
    timesA_C_G_T = 0

    if i[0][0]=='I':
        times1 = 0
        times2 = 0
        times3 = 0
        times4 = 0
        times5 = 0
        times = 0
        if i[2]>1:
            for j in range(i[1],i[2]+1):
                times1 += np.count_nonzero(array_with_aligned_sequences[:, j] == 'A')
                times2 += np.count_nonzero(array_with_aligned_sequences[:, j] == 'C')
                times3 += np.count_nonzero(array_with_aligned_sequences[:, j] == 'G')
                times4 += np.count_nonzero(array_with_aligned_sequences[:, j] == 'T')
                times5 += np.count_nonzero(array_with_aligned_sequences[:, j] == '_')

                times += (times1 + times2 + times3 + times4 + times5)

            times_A = times1
            times_C = times2
            times_G = times3
            times_T = times4

            timesA_C_G_T = times
        else:
            times_A = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == 'A')
            times_C = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == 'C')
            times_G = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == 'G')
            times_T = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == 'T')
            times_gap = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == '_')
            timesA_C_G_T = times_A + times_C + times_G + times_T + times_gap
        if (timesA_C_G_T == 0):
            timesA_C_G_T = 15
        indexing_emission.append(
            [i[0],round(times_A / timesA_C_G_T, 2), round(times_C / timesA_C_G_T, 2), round(times_G / timesA_C_G_T, 2),
             round(times_T / timesA_C_G_T, 2)])
    elif i[0][0]=='D':
        indexing_emission.append([i[0], 0, 0, 0, 0])
    else:
        times_A = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == 'A')
        times_C = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == 'C')
        times_G = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == 'G')
        times_T = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == 'T')
        times_gap = np.count_nonzero(array_with_aligned_sequences[:, i[1]] == '_')


        timesA_C_G_T = times_A + times_C + times_G + times_T + times_gap

        if(timesA_C_G_T == 0):
            timesA_C_G_T = 15

        indexing_emission.append(
            [i[0], round(times_A / timesA_C_G_T, 2), round(times_C / timesA_C_G_T, 2),round(times_G / timesA_C_G_T, 2),
             round(times_T / timesA_C_G_T, 2)])

#Αρχικοποίηση Start με πρόσβαση σε M1, D1 και I1
first_transition = -10
second_transition = -10
third_transition = -10
counter_initialiser = 0
for i in indexing:
    if i[0]=="M1":
        first_transition = indexing.index(i)
        counter_initialiser += 1

    if i[0]=="D1":
        second_transition = indexing.index(i)
        counter_initialiser += 1
    if i[0] == "I1":
        third_transition = indexing.index(i)
        counter_initialiser += 1
if (counter_initialiser == 1):
    if (first_transition >= 0):
        conditions[0][first_transition] = 1.0
    if (second_transition >= 0):
        conditions[0][second_transition] = 1.0
    if (third_transition >= 0):
        conditions[0][third_transition] = 1.0
if (counter_initialiser == 2):
    if (first_transition >= 0):
        conditions[0][first_transition] = 0.95
    if (second_transition >= 0):
        conditions[0][second_transition] = 0.05
    if (third_transition >= 0):
        conditions[0][third_transition] = 0.05
if (counter_initialiser == 3):
    if (first_transition >= 0):
        conditions[0][first_transition] = 0.95
    if (second_transition >= 0):
        conditions[0][second_transition] = 0.025
    if (third_transition >= 0):
        conditions[0][third_transition] = 0.025




print("HMM Hidden Profile Table: \n")
f = open("hmm_profile_table.txt", "w")
print("   ", [indexing[i][0] for i in range(0, len(indexing))])
f.write("       ")
f.write(",  ".join(map(str, [indexing[i][0] for i in range(0, len(indexing))])))
f.write("\n")
print("Sta.", conditions[0])
f.write("Sta. ")
size = len(conditions[0])
for j in range(size):
    if conditions[0][j] == 0:
        conditions[0][j] = 0.0
f.write(", ".join(map(str, conditions[0])))
f.write("\n")




counter = 0
for i in range(1,len(conditions)):
    print(indexing[counter][0], conditions[i])
    f.write(str(indexing[counter][0]))
    if (len(str(indexing[counter][0])) == 2):
        f.write("   ")
    else:
        f.write("  ")
    size = len(conditions[i])
    for j in range(size):
        if conditions[i][j] == 0:
            conditions[i][j] = 0.0
    f.write(", ".join(map(str, conditions[i])))
    f.write("\n")
    counter += 1

f.close()
print("\n \n")

emission_table = []
for i in indexing_emission:
    emission_table.append([i[1], i[2], i[3], i[4]])



hmm_states = []
for i in indexing:
    hmm_states.append(i[0])


print("DatasetB Sequences (Alignment_Score & Alignment_Path) : \n")
f = open("alignment_scores_&_paths.txt", "w")
for i in datasetB:
    val = viterbi_for_hmm(i, indexing,conditions, emission_table)
    f.write(", ".join(map(str, val)))
    f.write("\n")
    print(val)
f.close()

