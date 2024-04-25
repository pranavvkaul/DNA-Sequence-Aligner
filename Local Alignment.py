import numpy as np
import pandas as pd

def la(x, y, match=2, mismatch=-1, gap=-3):
    nx = len(x)
    ny = len(y)
    F = np.zeros((nx + 1, ny + 1), dtype=int) 
    F[:, 0] = 0
    F[0, :] = 0
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    # Variables to store the maximum score and its positions.
    max_score = 0           
    max_posn = []
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            if x[i-1] == y[j-1]:
                match_mismatch_score = match
            else:
                match_mismatch_score = mismatch

            diagonal_score = F[i - 1, j - 1] + match_mismatch_score
            top_score = F[i - 1, j] + gap
            left_score = F[i, j - 1] + gap
            # to choose max score
            if diagonal_score >= 0 and diagonal_score >= top_score and diagonal_score >= left_score:
                F[i, j] = diagonal_score
            elif top_score >= 0 and top_score >= diagonal_score and top_score >= left_score:
                F[i, j] = top_score
            elif left_score >= 0 and left_score >= diagonal_score and left_score >= top_score:
                F[i, j] = left_score
            else:
                F[i, j] = 0
            # Update max score and its pos.
            if F[i, j] > max_score:
                max_score = F[i, j]
                max_posn = [(i, j)]
            elif F[i, j] == max_score:
                max_posn.append((i, j))
            # Assign pointers based on the maximum score.
            if F[i, j] == diagonal_score:
                P[i, j] += 2
            if F[i, j] == top_score:
                P[i, j] += 3
            if F[i, j] == left_score:
                P[i, j] += 4
    print("Scoring Matrix:")
    row_labels = [' '] + list(x)
    column_labels = [' '] + list(y)
    df = pd.DataFrame(F, row_labels, column_labels)
    print(df)
    # Traceback from pos of the max score
    store = []
    def traceback(i, j, alignment_x, alignment_y):
        if F[i, j] == 0:
            store.append((alignment_x[::-1], alignment_y[::-1], max_score))
        else:
            if P[i, j] in [2, 5, 6, 9]:
                traceback(i - 1, j - 1, alignment_x + x[i - 1], alignment_y + y[j - 1])
            if P[i, j] in [3, 5, 7, 9]:
                traceback(i - 1, j, alignment_x + x[i - 1], alignment_y + '-')
            if P[i, j] in [4, 6, 7, 9]:
                traceback(i, j - 1, alignment_x + '-', alignment_y + y[j - 1])
    #print(max_posn)
    #print(max_score)
    # Traceback from pos of max score.
    for i in max_posn:
        traceback(i[0], i[1], '', '')
    print("All Possible Alignmnts:")
    for i in store:
        for j in i:
            print(j, end=" ")  
            print() 
    print()
    #print(store)
    #print the best alignment
    best_alignment = max(store, key=lambda x: x[2])
    print("Best Alignment with Score:")
    print("Score:", best_alignment[2])
    return best_alignment[0] + '\n' + best_alignment[1]

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA"
print("Optimal Local Alignmnts:")
print(la(seq1, seq2))