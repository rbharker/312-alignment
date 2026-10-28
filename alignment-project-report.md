# Project Report - Alignment

## Baseline

### Design Experience

*For this part of the project, I had to go over the homework
to understand the algorithm and how each value is decided. I showed my husband
how the algorithm works and explained that I wanted to use an array of arrays to
represent the table, but each entry would be a tuple of the value and the cell
it came from in order to be able to trace the optimized alignment path. I also
decided to use the min() function from python to determine which cell to pull
each value from.*

### Theoretical Analysis - Unrestricted Alignment

#### Time 

```commandline
def align(
    seq1: str,
    seq2: str,
    match_award=-3,
    indel_penalty=5,
    sub_penalty=1,
    banded_width=-1,
    gap="-",
) -> tuple[float, str | None, str | None]:
    """
    Align seq1 against seq2 using Needleman-Wunsch
    Put seq1 on left (j) and seq2 on top (i)
    => matrix[i][j]
    :param seq1: the first sequence to align; should be on the "left" of the matrix
    :param seq2: the second sequence to align; should be on the "top" of the matrix
    :param match_award: how many points to award a match
    :param indel_penalty: how many points to award a gap in either sequence
    :param sub_penalty: how many points to award a substitution
    :param banded_width: banded_width * 2 + 1 is the width of the banded alignment; -1 indicates full alignment
    :param gap: the character to use to represent gaps in the alignment strings
    """

    rows = len(seq1) + 1
    cols = len(seq2) + 1
    alignment_cost = [[(0, None) for j in range(rows)] for i in range(cols)]                                            # O(nm) - has to go through n cols and m rows

    # initialize base cases
    for i in range (cols):
        alignment_cost[i][0] = (indel_penalty*i, None)                                                                  # O(n) - has to go through n columns
    for j in range(rows):
        alignment_cost[0][j] = (indel_penalty*j, None)                                                                  # O(m) - has to go through m rows

    # traverse through the matrix until the final cost
    for j in range(1,rows):                                                                                             # O(m) - has to go through m rows
        for i in range(1, cols):                                                                                        # O(n) - has to go through n columns
            d = alignment_cost[i-1][j-1][0] + match_award if seq1[j-1] == seq2[i-1] \                                   # O(1) - constant lookup time
                else alignment_cost[i-1][j-1][0] + sub_penalty
            l = alignment_cost[i-1][j][0] + indel_penalty                                                               # O(1) - constant lookup time
            u = alignment_cost[i][j-1][0] + indel_penalty                                                               # O(1) - constant lookup time

            cost = min(d, l, u)                                                                                         # O(1) - comparison is constnat

            # priority order is diagonal, left, up
            if cost == d:                                                                                               # O(1) - comparison is constant
                alignment_cost[i][j] = (cost, "d")                                                                      # O(1) - declaration is constant
            elif cost == l:                                                                                             # O(1) - comparison is constant
                alignment_cost[i][j] = (cost, "l")                                                                      # O(1) - declaration is constant
            else:
                alignment_cost[i][j] = (cost, "u")                                                                      # O(1) - declaration is constant

    final_cost = alignment_cost[cols-1][rows-1][0]                                                                      # O(1) - lookup time is constant

    # align strings concurrently
    seq1_alignment = ""                                                                                                 # O(1) - declaration is constant
    seq2_alignment = ""                                                                                                 # O(1) - declaration is constant
    i = cols - 1                                                                                                        # O(1) - declaration is constant
    j = rows - 1                                                                                                        # O(1) - declaration is constant
    index1 = 0                                                                                                          # O(1) - declaration is constant
    index2 = 0                                                                                                          # O(1) - declaration is constant
    while alignment_cost[i][j][1] != None or i != 0 or j != 0:                                                          # O(1) - comparison is constant
        current = alignment_cost[i][j]                                                                                  # O(1) - lookup time is constant
        if current[1] == "d":                                                                                           # O(1) - comparison is constant
            index1 -= 1                                                                                                 # O(1) - declaration is constant
            index2 -= 1                                                                                                 # O(1) - declaration is constant
            seq1_alignment = seq1[index1] + seq1_alignment                                                              # O(1) - declaration is constant
            seq2_alignment = seq2[index2] + seq2_alignment                                                              # O(1) - declaration is constant
            i = i-1                                                                                                     # O(1) - declaration is constant
            j = j-1                                                                                                     # O(1) - declaration is constant
        elif current[1] == "l":                                                                                         # O(1) - comparison is constant
            index2 -= 1                                                                                                 # O(1) - declaration is constant
            seq1_alignment = gap + seq1_alignment                                                                       # O(1) - declaration is constant
            seq2_alignment = seq2[index2] + seq2_alignment                                                              # O(1) - declaration is constant
            i = i - 1                                                                                                   # O(1) - declaration is constant
            j = j                                                                                                       # O(1) - declaration is constant
        elif current[1] == "u":                                                                                         # O(1) - comparison is constant
            index1 -= 1                                                                                                 # O(1) - declaration is constant
            seq2_alignment = gap + seq2_alignment                                                                       # O(1) - declaration is constant
            seq1_alignment = seq1[index1] + seq1_alignment                                                              # O(1) - declaration and lookup is constant
            i = i                                                                                                       # O(1) - declaration is constant
            j = j - 1                                                                                                   # O(1) - declaration is constant
        else:
            if j > 0:                                                                                                   # O(1) - comparison is constant
                while j > 0:                                                                                            # O(1) - declaration is constant
                    index1 -= 1                                                                                         # O(1) - declaration is constant
                    seq2_alignment = gap + seq2_alignment                                                               # O(1) - declaration is constant
                    seq1_alignment = seq1[index1] + seq1_alignment                                                      # O(1) - declaration is constant
                    i = i                                                                                               # O(1) - declaration is constant
                    j = j - 1                                                                                           # O(1) - declaration is constant
            if i > 0:                                                                                                   # O(1) - comparison is constant
                while i > 0:                                                                                            # O(1) - comparison is constant
                    index2 -= 1                                                                                         # O(1) - declaration is constant
                    seq1_alignment = gap + seq1_alignment                                                               # O(1) - declaration is constant
                    seq2_alignment = seq2[index2] + seq2_alignment                                                      # O(1) - declaration is constant
                    i = i - 1                                                                                           # O(1) - declaration is constant
                    j = j                                                                                               # O(1) - declaration is constant

    return final_cost, seq1_alignment, seq2_alignment
```

***This algorithm runs with O(nm) time complexity** because at most it must go through n columns and m rows to
calculate the value at each cell.*

#### Space

```commandline
def align(
    seq1: str,
    seq2: str,
    match_award=-3,
    indel_penalty=5,
    sub_penalty=1,
    banded_width=-1,
    gap="-",
) -> tuple[float, str | None, str | None]:
    """
    Align seq1 against seq2 using Needleman-Wunsch
    Put seq1 on left (j) and seq2 on top (i)
    => matrix[i][j]
    :param seq1: the first sequence to align; should be on the "left" of the matrix
    :param seq2: the second sequence to align; should be on the "top" of the matrix
    :param match_award: how many points to award a match
    :param indel_penalty: how many points to award a gap in either sequence
    :param sub_penalty: how many points to award a substitution
    :param banded_width: banded_width * 2 + 1 is the width of the banded alignment; -1 indicates full alignment
    :param gap: the character to use to represent gaps in the alignment strings
    """

    rows = len(seq1) + 1
    cols = len(seq2) + 1
    alignment_cost = [[(0, None) for j in range(rows)] for i in range(cols)]                                            # O(nm) - has to go allocate n cols and m rows

    # initialize base cases
    for i in range (cols):
        alignment_cost[i][0] = (indel_penalty*i, None)                                                                  # O(1) - space already allocated
    for j in range(rows):
        alignment_cost[0][j] = (indel_penalty*j, None)                                                                  # O(1) - space already allocated

    # traverse through the matrix until the final cost
    for j in range(1,rows):                                                                                             
        for i in range(1, cols):                                                                                        
            d = alignment_cost[i-1][j-1][0] + match_award if seq1[j-1] == seq2[i-1] \                                   # O(1) - single variable declaration is constant
                else alignment_cost[i-1][j-1][0] + sub_penalty
            l = alignment_cost[i-1][j][0] + indel_penalty                                                               # O(1) - single variable declaration is constant
            u = alignment_cost[i][j-1][0] + indel_penalty                                                               # O(1) - single variable declaration is constant

            cost = min(d, l, u)                                                                                         # O(1) - single variable declaration is constant

            # priority order is diagonal, left, up
            if cost == d:                                                                                               
                alignment_cost[i][j] = (cost, "d")                                                                      # O(1) - space already allocated
            elif cost == l:                                                                                             
                alignment_cost[i][j] = (cost, "l")                                                                      # O(1) - space already allocated
            else:
                alignment_cost[i][j] = (cost, "u")                                                                      # O(1) - space already allocated

    final_cost = alignment_cost[cols-1][rows-1][0]                                                                      # O(1) - single variable declaration is constant

    # align strings concurrently
    seq1_alignment = ""                                                                                                 # O(1) - declaration is constant
    seq2_alignment = ""                                                                                                 # O(1) - declaration is constant
    i = cols - 1                                                                                                        # O(1) - declaration is constant
    j = rows - 1                                                                                                        # O(1) - declaration is constant
    index1 = 0                                                                                                          # O(1) - declaration is constant
    index2 = 0                                                                                                          # O(1) - declaration is constant
    while alignment_cost[i][j][1] != None or i != 0 or j != 0:                                                          
        current = alignment_cost[i][j]                                                                                  # O(1) - lookup time is constant
        if current[1] == "d":                                                                                           
            index1 -= 1                                                                                                 # O(1) - declaration is constant
            index2 -= 1                                                                                                 # O(1) - declaration is constant
            seq1_alignment = seq1[index1] + seq1_alignment                                                              # O(1) - declaration is constant
            seq2_alignment = seq2[index2] + seq2_alignment                                                              # O(1) - declaration is constant
            i = i-1                                                                                                     # O(1) - declaration is constant
            j = j-1                                                                                                     # O(1) - declaration is constant
        elif current[1] == "l":                                                                                         # O(1) - comparison is constant
            index2 -= 1                                                                                                 # O(1) - declaration is constant
            seq1_alignment = gap + seq1_alignment                                                                       # O(1) - declaration is constant
            seq2_alignment = seq2[index2] + seq2_alignment                                                              # O(1) - declaration is constant
            i = i - 1                                                                                                   # O(1) - declaration is constant
            j = j                                                                                                       # O(1) - declaration is constant
        elif current[1] == "u":                                                                                         # O(1) - comparison is constant
            index1 -= 1                                                                                                 # O(1) - declaration is constant
            seq2_alignment = gap + seq2_alignment                                                                       # O(1) - declaration is constant
            seq1_alignment = seq1[index1] + seq1_alignment                                                              # O(1) - declaration and lookup is constant
            i = i                                                                                                       # O(1) - declaration is constant
            j = j - 1                                                                                                   # O(1) - declaration is constant
        else:
            if j > 0:                                                                                                   # O(1) - comparison is constant
                while j > 0:                                                                                            # O(1) - comparison is constant
                    index1 -= 1                                                                                         # O(1) - declaration is constant
                    seq2_alignment = gap + seq2_alignment                                                               # O(1) - declaration is constant
                    seq1_alignment = seq1[index1] + seq1_alignment                                                      # O(1) - declaration is constant
                    i = i                                                                                               # O(1) - declaration is constant
                    j = j - 1                                                                                           # O(1) - declaration is constant
            if i > 0:                                                                                                   # O(1) - comparison is constant
                while i > 0:                                                                                            # O(1) - comparison is constant
                    index2 -= 1                                                                                         # O(1) - declaration is constant
                    seq1_alignment = gap + seq1_alignment                                                               # O(1) - declaration is constant
                    seq2_alignment = seq2[index2] + seq2_alignment                                                      # O(1) - declaration is constant
                    i = i - 1                                                                                           # O(1) - declaration is constant
                    j = j                                                                                               # O(1) - declaration is constant

    return final_cost, seq1_alignment, seq2_alignment
```

***This algorithm allocates O(nm) space** because it needs one table with n columns and n rows.*

### Empirical Data - Unrestricted Alignment

| Size | Time (sec) |
| ---- | ---------- |
| 500  | 0.087      |
| 1000 | 0.352      |
| 1500 | 0.799      |
| 2000 | 1.474      |
| 2500 | 2.28       |
| 3000 | 3.312      |


### Comparison of Theoretical and Empirical Results - Unrestricted Alignment

- Theoretical order of growth: O(nm)
- Empirical order of growth (if different from theoretical): O(n)

![](/empirical_analysis/unbanded_align.svg)

*The empirical data shows that this algorithm runs faster on average than I expected it to. I predict
that is because many of the calculations needed for each cell are trivial, so it ends up being faster
O(nm).*

## Project Review

*Fill me in*
