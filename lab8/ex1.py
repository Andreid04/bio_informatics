def needleman_wunsch(
    s1,
    s2,
    match=1,
    mismatch=-1,
    gap=-2
):

    rows = len(s1) + 1
    cols = len(s2) + 1

    # scoring matrix
    matrix = [
        [0 for _ in range(cols)]
        for _ in range(rows)
    ]

    for i in range(rows):
        matrix[i][0] = i * gap

    for j in range(cols):
        matrix[0][j] = j * gap

    for i in range(1, rows):
        for j in range(1, cols):

            if s1[i - 1] == s2[j - 1]:
                diagonal = matrix[i - 1][j - 1] + match
            else:
                diagonal = matrix[i - 1][j - 1] + mismatch

            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap

            matrix[i][j] = max(
                diagonal,
                up,
                left
            )

    aligned_s1 = []
    aligned_s2 = []

    i = len(s1)
    j = len(s2)

    while i > 0 or j > 0:

        if i > 0 and j > 0:

            if s1[i - 1] == s2[j - 1]:
                score = match
            else:
                score = mismatch

            if matrix[i][j] == matrix[i - 1][j - 1] + score:

                aligned_s1.append(s1[i - 1])
                aligned_s2.append(s2[j - 1])

                i -= 1
                j -= 1

                continue

        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap:

            aligned_s1.append(s1[i - 1])
            aligned_s2.append("-")

            i -= 1

        else:

            aligned_s1.append("-")
            aligned_s2.append(s2[j - 1])

            j -= 1

    # traceback goes backwards
    aligned_s1.reverse()
    aligned_s2.reverse()

    aligned_s1 = "".join(aligned_s1)
    aligned_s2 = "".join(aligned_s2)

    return matrix, aligned_s1, aligned_s2


def display_matrix(matrix, s1, s2):

    print("\nscoring matrix")
    print("-" * 10)

    print("      -", end="")

    for base in s2:
        print(f"{base:>5}", end="")

    print()

    for i, row in enumerate(matrix):

        if i == 0:
            label = "-"
        else:
            label = s1[i - 1]

        print(f"{label:>5}", end="")

        for value in row:
            print(f"{value:>5}", end="")

        print()


def display_alignment(aligned_s1, aligned_s2):

    print("\nalignment")
    print("-" * 10)

    print("S1:", aligned_s1)

    match_line = ""

    for a, b in zip(aligned_s1, aligned_s2):

        if a == b:
            match_line += "|"
        elif a == "-" or b == "-":
            match_line += " "
        else:
            match_line += "."

    print("    ", match_line)
    print("S2:", aligned_s2)

    matches = sum(
        1 for a, b in zip(aligned_s1, aligned_s2)
        if a == b
    )

    mismatches = sum(
        1 for a, b in zip(aligned_s1, aligned_s2)
        if a != b and a != "-" and b != "-"
    )

    gaps = sum(
        1 for a, b in zip(aligned_s1, aligned_s2)
        if a == "-" or b == "-"
    )

    print("\nalignment statistics")
    print("-" * 10)

    print("Alignment length:", len(aligned_s1))
    print("Matches:", matches)
    print("Mismatches:", mismatches)
    print("Gaps:", gaps)

    identity = (matches / len(aligned_s1)) * 100

    print(f"Sequence identity: {identity:.2f}%")


def main():

    print("Needleman-Wunsch DNA Alignment")
    print("-" * 10)

    s1 = "ACCGTGAAGCCAATAC"
    s2 = "AGCGTGCAGCCAATAC"

    print("\nINPUT SEQUENCES")
    print("-" * 10)

    print("S1:", s1)
    print("S2:", s2)

    match = 1
    mismatch = -1
    gap = -2

    print("\nscoring parameters")
    print("-" * 10)

    print("Match:", match)
    print("Mismatch:", mismatch)
    print("Gap:", gap)

    matrix, aligned_s1, aligned_s2 = needleman_wunsch(
        s1,
        s2,
        match,
        mismatch,
        gap
    )

    display_matrix(matrix, s1, s2)

    display_alignment(aligned_s1, aligned_s2)

    print("\noptimal alignment score")
    print("-" * 10)

    print(matrix[len(s1)][len(s2)])


if __name__ == "__main__":
    main()
