from cmath import inf


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

    if banded_width != -1 and (len(seq1) > 2*banded_width+1) ^ (len(seq2) > 2*banded_width+1):
        return float('inf'), None, None

    rows = len(seq1) + 1
    cols = len(seq2) + 1
    alignment_cost = []
    for i in range(cols):
        row = []
        for j in range(rows):
            if banded_width != -1 and abs(i-j) > banded_width:
                continue
            row.append((0, None))
        if len(row) > 0:
            alignment_cost.append(row)

    final_alignment_cost = find_final_cost(seq1, seq2, rows, cols, alignment_cost,
                    match_award, indel_penalty, sub_penalty, banded_width)

    if banded_width != -1:
        final_cost = final_alignment_cost[cols-1][-1][0]
    else:
        final_cost = final_alignment_cost[cols-1][rows-1][0]

    seq1_alignment, seq2_alignment = align_strings(seq1, seq2, rows, cols,
                  final_alignment_cost, banded_width, gap)

    return final_cost, seq1_alignment, seq2_alignment

def find_final_cost(
    seq1: str,
    seq2: str,
    rows: int,
    cols: int,
    alignment_cost: list,
    match_award=-3,
    indel_penalty=5,
    sub_penalty=1,
    banded_width=-1)->list:

    # initialize base cases
    alignment_cost[0][0] = (0, None)
    for i in range(1, cols):
        if banded_width == -1 or i <= banded_width:
            alignment_cost[i][0] = (indel_penalty * i, "l")
    for j in range(1, rows):
        if banded_width == -1 or j <= banded_width:
            alignment_cost[0][j] = (indel_penalty * j, "u")


    # traverse through the matrix until the final cost
    for i in range(1, cols):
        for j in range(1, rows):

            # if you're in the case outside of banded width, skip
            if banded_width != -1 and abs(i-j) > banded_width:
                continue

            # variable for adjusting for blank spots
            if banded_width != -1 and i > banded_width:
                usable_j = j-(i-banded_width)
            else:
                usable_j = j

            # if you've adjusted for blank spots, diagonal coordinates are different
            if (usable_j != j):
                d = alignment_cost[i - 1][usable_j][0] + match_award if seq1[j - 1] == seq2[i - 1] \
                    else alignment_cost[i - 1][usable_j][0] + sub_penalty
            # otherwise use regular diagonal coordinates
            else:
                d = alignment_cost[i - 1][usable_j - 1][0] + match_award if seq1[j - 1] == seq2[i - 1] \
                    else alignment_cost[i - 1][usable_j - 1][0] + sub_penalty

            # if you're not going to look in blank spots
            if banded_width == -1 or abs((i-1)-j) <= banded_width:
                if (usable_j != j):
                    l = alignment_cost[i - 1][usable_j + 1][0] + indel_penalty
                else:
                    l = alignment_cost[i - 1][usable_j][0] + indel_penalty
            # otherwise don't use left value as minimum
            else:
                l = float('inf')

            # if you're not going to look in blank spots
            if (usable_j != 0):
                u = alignment_cost[i][usable_j - 1][0] + indel_penalty
            # otherwise don't use above value as minimum
            else:
                u = float('inf')

            cost = min(d, l, u)
            if cost == d:
                alignment_cost[i][usable_j] = (cost, "d")
            elif cost == l:
                alignment_cost[i][usable_j] = (cost, "l")
            elif cost == u:
                alignment_cost[i][usable_j] = (cost, "u")

    return alignment_cost

def align_strings(
    seq1: str,
    seq2: str,
    rows: int,
    cols: int,
    alignment_cost: list,
    banded_width=-1,
    gap="-")->tuple[str, str]:

    seq1_alignment = ""
    seq2_alignment = ""
    i = -1
    j = -1
    index1 = -1
    index2 = -1

    change = True
    while change == True:
        change= False
        current = alignment_cost[i][j]
        if current[1] == "d":
            if abs(index1) < rows:
                seq1_alignment = seq1[index1] + seq1_alignment
            if abs(index2) < cols:
                seq2_alignment = seq2[index2] + seq2_alignment
            if banded_width != -1:
                if abs(i) <= banded_width:
                    j = j - 1
                else:
                    j = j
            else:
                j = j - 1
            i = i - 1
            index1 -= 1
            index2 -= 1
            change = True
        elif current[1] == "l":
            seq1_alignment = gap + seq1_alignment
            if abs(index2) < cols:
                seq2_alignment = seq2[index2] + seq2_alignment
            if banded_width != -1:
                if abs(i) <= banded_width:
                    j = j
                else:
                    j = j + 1
            else:
                j = j
            i = i - 1
            index2 -= 1
            change = True
        elif current[1] == "u":
            seq2_alignment = gap + seq2_alignment
            if abs(index1) < rows:
                seq1_alignment = seq1[index1] + seq1_alignment
            i = i
            j = j - 1
            index1 -= 1
            change = True
        else:
            if alignment_cost[i][j] != alignment_cost[i][0]:
                while alignment_cost[i][j] != alignment_cost[i][0]:
                    seq2_alignment = gap + seq2_alignment
                    if abs(index1) < rows:
                        seq1_alignment = seq1[index1] + seq1_alignment
                    i = i
                    j = j - 1
                    index1 -= 1
            if abs(i) < cols:
                while abs(i) < cols:
                    seq1_alignment = gap + seq1_alignment
                    if abs(index2) < cols:
                        seq2_alignment = seq2[index2] + seq2_alignment
                    i = i - 1
                    j = j
                    index2 -= 1
            break
    return seq1_alignment, seq2_alignment

def local_align(
    seq1: str,
    seq2: str,
    match_award=-3,
    indel_penalty=5,
    sub_penalty=3,
    banded_width=-1,
    gap="-",
) -> tuple[float, str | None, str | None]:
    """
    Align seq1 against seq2 using Smith-Waterman
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

    return 0, "", ""
