import numpy as np
import pandas as pd

def nw(x, y, match=2, mismatch=-3, gap=-1):  # needleman-wunsch
    nx = len(x) #storing the length of seq1-x    
    ny = len(y) #storing the length of seq2-y  
    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1), dtype=int) 
    F[:, 0]  = np.linspace(0, gap * nx, nx + 1) #will name first col ele from 0 to -9
    F[0, :] = np.linspace(0, gap * ny, ny + 1)  #will name first row ele from 0 to -7
    # Pointers to trace through an optimal alignment.
    # Matrix filling.
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3  #initialize first col with 3
    P[0, :] = 4  #then initialize first row with 4 , [0][0]also become4
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            if x[i-1] == y[j-1]:
                match_mismatch_score = match #diag ele check for match if equal
                #print(x[i-1] + y[j-1])
            else:
                match_mismatch_score = mismatch
                #print(x[i-1] + y[j-1])
            #getting different scores for top, left and diagonal
            diagonal_score = F[i - 1, j - 1] + match_mismatch_score #get val at prev diag and add/sub match mismatch
            top_score = F[i - 1, j] + gap #getting val at position just above the row 
            left_score = F[i, j - 1] + gap #getting val at position just left side of the col 
            # finding the maximum score for F[i][j]
            if diagonal_score >= top_score and diagonal_score >= left_score:
                F[i, j] = diagonal_score
            elif top_score >= diagonal_score and top_score >= left_score:
                F[i, j] = top_score
            else:
                F[i, j] = left_score
            # Assign pointers based on the maximum score.
            if F[i, j] == diagonal_score:
                P[i, j] += 2 
            if F[i, j] == top_score:
                P[i, j] += 3
            if F[i, j] == left_score:
                P[i, j] += 4
    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    row_labels = [''] + list(x) #first place empty
    column_labels = [' '] + list(y) #first place empty
    df = pd.DataFrame(F, row_labels, column_labels)
    print(df)

    # Trace through all optimal alignmnts.
            #2 to the diagonal dir.
            #3 to the vertical dir.
            #4 to the horizontal dir.
            #5 to both diagonal and vertical dirs.
            #6 to both diagonal and horizontal dirs.
            #7 to both vertical and horizontal dirs.
            #9 to diagonal, vertical, and horizontal dirs
    store = [] #to store all optical alignmnts in a list
    def traceback(i, j, alignment_x, alignment_y):
        if i == 0 and j == 0: #to check if both i and j have reached the top left corner of the traceback matrix P
            store.append((alignment_x[::-1], alignment_y[::-1], F[-1][-1])) #reversal as the chars are appended at the end so in order to make them in correct order
        else:
            if P[i, j] in [2, 5, 6, 9]:
                traceback(i - 1, j - 1, alignment_x + x[i - 1], alignment_y + y[j - 1])
            if P[i, j] in [3, 5, 7, 9]:
                traceback(i - 1, j, alignment_x + x[i - 1], alignment_y + '-')
            if P[i, j] in [4, 6, 7, 9]:
                traceback(i, j - 1, alignment_x + '-', alignment_y + y[j - 1])
    traceback(nx, ny, '', '')
    print("All Possible Alignments:")
    for i in store:
        for j in i:
            print(j, end=" ")  
            print() 
    print()

    print("Best Alignment with Score:")
    print("Score=", store[0][2])
    return store[0][0] + '\n' + store[0][1]

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA"
print(nw(seq1, seq2))