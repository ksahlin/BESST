'''
Created on Jun 17, 2013

@author: ksahlin
'''



def traceback_parser(score, row, column, traceback_matrix, seq1, seq2):

    s_i = []
    s_j = []

    stepcounter = 0
    while True:
        stepcounter += 1
        if traceback_matrix[row][column] == 1:
            s_i.append(seq1[row - 1])
            s_j.append(seq2[column - 1])
            row = row - 1
            column = column - 1
        elif traceback_matrix[row][column] == 3:
            s_i.append('-')
            s_j.append(seq2[column - 1])
            column = column - 1
        elif traceback_matrix[row][column] == 2:
            s_i.append(seq1[row - 1])
            s_j.append('-')
            row = row - 1
        else:
            break

    s_i.reverse()
    s_j.reverse()

    return -(stepcounter - 1) # s_i,s_j

def score(a, b, match_score, mismatch_penalty):
    if a == b:
        s = match_score
    if a != b:
        s = mismatch_penalty

    return s

def zero_maker(n, m):
    res_row = []
    res_all = []

    for i in range(n):
        temp = []
        for j in range(m):
            temp.append(0)
        res_all.append(temp)
        del(temp)
    return res_all

def SW(seq1, seq2, gap_penalty, mismatch_penalty):
    score_matrix = zero_maker(len(seq1) + 1, len(seq2) + 1)
    traceback_matrix = zero_maker(len(seq1) + 1, len(seq2) + 1)

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            substitution_score = score_matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1], 1, mismatch_penalty)
            i_indel = score_matrix[i - 1][j] + gap_penalty
            j_indel = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(substitution_score, i_indel, j_indel)

            if score_matrix[i][j] == substitution_score:
                traceback_matrix[i][j] = 1
            elif score_matrix[i][j] == i_indel:
                traceback_matrix[i][j] = 2
            elif score_matrix[i][j] == j_indel:
                traceback_matrix[i][j] = 3

    n = len(score_matrix)     #nr of rows
    m = len(score_matrix[0])     #nr or columns
    max_score = 0

    #if in last column
    for i in range(0, n):
        if score_matrix[i][m - 1] >= max_score:
            max_i = i
            max_j = m - 1
            max_score = score_matrix[i][m - 1]

    #if in last row
    for j in range(0, m):
        if score_matrix[n - 1][j] >= max_score:
                max_i = n - 1
                max_j = j
                max_score = score_matrix[n - 1][j]


    #Returnera score har

    return max_score, max_i, max_j


if __name__ == '__main__':
    seq1 = 'TTAGCCACTG' * 20
    seq2 = 'TTTGCTGGAG' * 20
    score, row, column, traceback_matrix = SW(seq1, seq2, -1, -3)
    b = traceback_parser(score, row, column, traceback_matrix, seq1, seq2)
    print b
    #print score_matrix
#print traceback_matrix
