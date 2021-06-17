from queue import LifoQueue
from concurrent.futures import *
from copy import deepcopy

MI_Matrix = None

from math import log2
from copy import copy
from pandas import *


BASES = ['A', 'U','G','C']


#  m    m mmmmm
#  ##  ##   #
#  # ## #   #
#  # "" #   #
#  #    # mm#mm


def freq_arr(seq_arr):
    len_dict = []
    empty_pos = {'A':0.0, 'U':0.0, 'G':0.0, 'C':0.0, 'SUM':0.0}
    for j in range(len(seq_arr[0])):
        tmp_dict = copy(empty_pos)
        for seq in seq_arr:
            tmp_dict[seq[j]] += 1.0
            tmp_dict['SUM'] += 1.0
        len_dict.append(tmp_dict)
    return len_dict

            
def MI(pos1, pos2,seq_arr, fr_arr):
    pos1 -= 1
    pos2 -= 1
    MI_sum = 0.0
    double_exist = 0.0
    cur_bp = []
    for seq in seq_arr:
        double_exist = 0
        cur_bp = [seq[pos1], seq[pos2]]
        print("For ", cur_bp)
        for seq_inner in seq_arr:
            if seq_inner[pos1] == seq[pos1] and seq_inner[pos2] == seq[pos2]:
                double_exist += 1.0
        double_freq = ((double_exist)/len(seq_arr))
        pos1_freq = (fr_arr[pos1][seq[pos1]] / fr_arr[pos1]['SUM'])
        pos2_freq = (fr_arr[pos2][seq[pos2]] / fr_arr[pos2]['SUM'])
        print("Adding {} times log2 of {} divided by {} times {}".format(double_freq, double_freq, pos1_freq, pos2_freq))
        print("Adding ", (double_freq) * log2((double_freq)/(pos1_freq * pos2_freq)))
        MI_sum += (double_freq) * log2((double_freq)/(pos1_freq * pos2_freq))
    return MI_sum


def create_MI_mat(seq_arr):
    str_len = len(seq_arr[0]) + 1
    m_freq = freq_arr(seq_arr)
    MI_mat = [[0 for i in range(str_len)] for j in range(str_len)]
    for i in range(1, str_len - 1):
        for j in range(i+1, str_len):
            MI_mat[i][j] = MI(i, j, seq_arr, m_freq)
    return MI_mat


#  m    m mmmmm                       mm   m m    m  mmmm   mmmm  mmmmm  mm   m
#  ##  ##   #             m           #"m  # #    # #"   " #"   "   #    #"m  #
#  # ## #   #             #           # #m # #    # "#mmm  "#mmm    #    # #m #
#  # "" #   #          """#"""        #  # # #    #     "#     "#   #    #  # #
#  #    # mm#mm           #           #   ## "mmmm" "mmm#" "mmm#" mm#mm  #   ##



def get_mi_result(x, y, i, j, grade_mat):
    #if i == j-1:
    #    return 0
    return grade_mat.at[i,j]

def MI_compute_matrix(arr, my_seq, mi_df):
    length = len(my_seq)
    count = length-1  # the amount of diagonals
    while count != 0:  # running on length-1 diagonals
        distance = length - count  # num of cells in the diagonal
        for i in range(length-distance):
            j = i+distance
            val = get_mi_result(my_seq[i], my_seq[j], i + 1, j + 1, mi_df)  # check if it is a base pair
            val1 = 0
            val2 = 0
            val3 = 0

            # Got MI fit
            if val != 0:  
                val1 = val + arr[i+1][j-1]
           
            val3 = get_maximal_val(arr, i, j)
            arr[i][j] = max(val1, val2, val3)
        count -= 1
    for i in range(length):
        print(arr[i])

    reconstruct_MI(arr, my_seq, mi_df)

def reconstruct_MI(arr, my_seq, mi_df):
    count = 0
    length = len(my_seq)
    stack = LifoQueue(length)
    stack.put([0, length-1])
    while not stack.empty():
        curr_cell = stack.get()
        i = curr_cell[0]
        j = curr_cell[1]
        if i < j:
            if arr[i][j] == arr[i+1][j-1] + get_mi_result(my_seq[i], my_seq[j], i + 1, j + 1, mi_df):
                print('[' + my_seq[i] + ' , ' + my_seq[j] + ']', i, j)
                count += 2
                stack.put([i + 1, j - 1])
            else:
                k = i
                while k < j:
                    if arr[i][j] == arr[i][k] + arr[k + 1][j]:
                        stack.put([k + 1, j])
                        stack.put([i, k])
                        break
                    k += 1
        elif(i == j):
            print(my_seq[i], i)
            count += 1


#  mm   m m    m  mmmm   mmmm  mmmmm  mm   m  mmmm  m    m
#  #"m  # #    # #"   " #"   "   #    #"m  # m"  "m "m  m"
#  # #m # #    # "#mmm  "#mmm    #    # #m # #    #  #  #
#  #  # # #    #     "#     "#   #    #  # # #    #  "mm"
#  #   ## "mmmm" "mmm#" "mmm#" mm#mm  #   ##  #mm#    ##

def compute_matrix(arr, my_seq):
    length = len(my_seq)
    count = length-1  # the amount of diagonals
    while count != 0:  # running on length-1 diagonals
        distance = length - count  # num of cells in the diagonal
        for i in range(length-distance):
            j = i+distance
            val = is_base_pair(my_seq[i], my_seq[j], i, j)  # check if it is a base pair
            val1 = 0
            val2 = 0
            val3 = 0
            if val != 0:  # if BP
                val1 = val + arr[i+1][j-1]
           
            val3 = get_maximal_val(arr, i, j)
            arr[i][j] = max(val1, val2, val3)
        count -= 1
    for i in range(length):
        print(arr[i])

    reconstruct_list(arr, my_seq)

def reconstruct_orig(arr, my_seq):
    count = 0
    length = len(my_seq)
    stack = LifoQueue(length)
    stack.put([0, length-1])
    while not stack.empty():
        curr_cell = stack.get()
        i = curr_cell[0]
        j = curr_cell[1]
        if i < j:
            if arr[i][j] == arr[i+1][j-1] + is_base_pair(my_seq[i], my_seq[j], i, j) and is_base_pair(my_seq[i], my_seq[j], i, j) == 1:
                print(my_seq[i] + ' , ' + my_seq[j], i, j)
                count += 2
                stack.put([i + 1, j - 1])
            else:
                k = i
                while k < j:
                    if arr[i][j] == arr[i][k] + arr[k + 1][j]:
                        stack.put([k + 1, j])
                        stack.put([i, k])
                        break
                    k += 1
        elif(i == j):
            print(my_seq[i], i)
            count += 1

def get_result_len(result_list):
    count = 0;
    for it in result_list:
        if len(it) == 4:
            count += 2
        if len(it) == 2:
            count += 1
    return count

def reconstruct_list(arr, my_seq, stack=None, count=None, result_stack = None):
    if count is None:
        count = 0
   
    if stack is None:
        length = len(my_seq)
        stack = []
        stack.append([0, length-1])

    if result_stack is None:
        result_stack = []
        
    while  stack:
        curr_cell = stack.pop()
        i = curr_cell[0]
        j = curr_cell[1]
        if i < j:
            # Checking to see if its a diagonal fit
            if arr[i][j] == arr[i+1][j-1] + is_base_pair(my_seq[i], my_seq[j], i, j):      
                # Coping the record stack
                result_cpy = deepcopy(result_stack)
                result_cpy.append([my_seq[i], my_seq[j], i + 1, j + 1])
                count += 2

                # Coping the pairs stack
                stack_cpy = deepcopy(stack)
                stack_cpy.append([i + 1, j - 1])
                reconstruct_list(arr, my_seq, stack_cpy, count, result_cpy)
            #else:
            k = i
            while k < j:
                # Maximal value fit
                if arr[i][j] == arr[i][k] + arr[k + 1][j]:
                    # Coping both stacks
                    stack_cpy = deepcopy(stack)
                    stack_cpy.append([k + 1, j])
                    stack_cpy.append([i, k])
                    result_cpy = deepcopy(result_stack)
                    reconstruct_list(arr, my_seq, stack_cpy, count, result_cpy)
                k += 1
        elif(i == j):
            result_stack.append([my_seq[i], i + 1])
            count += 1

    # Because I didn't in the final recursion of each loop, I got a lot of semi full stacks.
    # removing with an if condition. isn't pretty but it does the job
    if (result_stack) and (len(my_seq) == get_result_len(result_stack)) :
        print("======================================")
        while (result_stack):
            print(result_stack.pop())
        print("=================\n\n\n")

def get_maximal_val(arr, i, j):
    maximum = 0
    for k in range(i, j):
        val = arr[i][k] + arr[k+1][j]
        if val > maximum:
            maximum = val
    return maximum

def is_base_pair(x, y, i, j):
    if i == j:
        return -10000
    elif x == 'C' and y == 'G':
        return 1
    elif x == 'G' and y == 'C':
        return 1
    elif x == 'A' and y == 'U':
        return 1
    elif x == 'U' and y == 'A':
        return 1
    else:
        return -10000

# initialize the matrix with zeros
def get_matrix(my_length):
    arr = [[0 for i in range(my_length)] for j in range(my_length)]
    return arr


# a_MI_mat = create_MI_mat(a)
# df_a = DataFrame(a_MI_mat)
# df_a_short = df_a.drop(0,axis=0).drop(0,axis=1)
# init_mat = get_matrix(len("123456"))
# MI_compute_matrix(init_mat, "123456", df_a_short)

seq = input("Enter the RNA sequence : ")
initialized_matrix = get_matrix(len(seq))
compute_matrix(initialized_matrix, seq)