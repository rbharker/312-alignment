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
    alignment_cost = []
    for i in range(cols):
        row = []
        for j in range(rows):
            if banded_width != -1 and abs(j-i) > banded_width:
                continue
            row.append((0, None))
        alignment_cost.append(row)

    final_alignment_cost = find_final_cost(seq1, seq2, rows, cols, alignment_cost,
                    match_award, indel_penalty, sub_penalty, banded_width)

    final_cost = alignment_cost[cols-1][rows-1][0]

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
    for i in range(cols):
        if banded_width == -1 or i < banded_width + 1:
            alignment_cost[i][0] = (indel_penalty * i, None)
    for j in range(rows):
        if banded_width == -1 or j < banded_width + 1:
            alignment_cost[0][j] = (indel_penalty * j, None)


    # traverse through the matrix until the final cost
    for j in range(1, rows):
        for i in range(1, cols):
            if banded_width != -1 and abs(j-i) > banded_width - 1:
                continue
            d = alignment_cost[i - 1][j - 1][0] + match_award if seq1[j - 1] == seq2[i - 1] \
                else alignment_cost[i - 1][j - 1][0] + sub_penalty
            l = alignment_cost[i - 1][j][0] + indel_penalty
            u = alignment_cost[i][j - 1][0] + indel_penalty

            cost = min(d, l, u)
            if cost == d:
                alignment_cost[i][j] = (cost, "d")
            elif cost == l:
                alignment_cost[i][j] = (cost, "l")
            else:
                alignment_cost[i][j] = (cost, "u")

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
    i = cols - 1
    j = rows - 1
    index1 = 0
    index2 = 0

    while ((alignment_cost[i][j][1] is not None or i != 0 or j != 0)
           and (banded_width == -1 or abs(j-i) < banded_width)):
        current = alignment_cost[i][j]
        if current[1] == "d":
            index1 -= 1
            index2 -= 1
            seq1_alignment = seq1[index1] + seq1_alignment
            seq2_alignment = seq2[index2] + seq2_alignment
            i = i - 1
            j = j - 1
        elif current[1] == "l":
            index2 -= 1
            seq1_alignment = gap + seq1_alignment
            seq2_alignment = seq2[index2] + seq2_alignment
            i = i - 1
            j = j
        elif current[1] == "u":
            index1 -= 1
            seq2_alignment = gap + seq2_alignment
            seq1_alignment = seq1[index1] + seq1_alignment
            i = i
            j = j - 1
        else:
            if j > 0:
                while j > 0:
                    index1 -= 1
                    seq2_alignment = gap + seq2_alignment
                    seq1_alignment = seq1[index1] + seq1_alignment
                    i = i
                    j = j - 1
            if i > 0:
                while i > 0:
                    index2 -= 1
                    seq1_alignment = gap + seq1_alignment
                    seq2_alignment = seq2[index2] + seq2_alignment
                    i = i - 1
                    j = j
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
