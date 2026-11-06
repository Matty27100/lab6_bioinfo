dna = (
    "GTTCTGGAACGCGCTTGACTACGTTAGTGCATCTATTTACATCTCTTAGTGCCTAGGGTA"
    "ACGTGACTCCCGGCATTAGGCGTGCACAAATGTTTACTAGGGCATAGCACGCCCGGCTTA"
    "CACAGGATGCTTCGACCCGGGAGGTCCGGCACTTCTACGGACTGTATTGTCGAAATGCAG"
    "CCGTTCTGACCAGTGCGACCAAGCACATTACCGCTATGCCCTAGCGAGTCCTTTCACAAA"
    "CCTCATACGCACTGTGACTTGTTAATTCCACTACAGAGCCCCTTTAACTTAGCCAGTTCC"
)

def detect_transposons(dna, ir=6, dr=4):
    positions = []
    complement = str.maketrans("ACGT", "TGCA")
    for start in range(dr, len(dna) - ir):
        repeat = dna[start:start + ir].translate(complement)[::-1]
        for end in range(start + 2 * ir + 1, len(dna) - dr + 1):
            if dna[end - ir:end] == repeat and dna[start - dr:start] == dna[end:end + dr]:
                positions.append((start + 1, end))
    return positions

print("DNA length:", len(dna))
print("DNA:", dna)
print("Positions (1-based, inclusive):")
elements = detect_transposons(dna)

for i, (start, end) in enumerate(elements):
    print(f"T{i + 1}: start={start}, end={end}")

for i in range(len(elements)):
    for j in range(i + 1, len(elements)):
        a, b = elements[i]
        c, d = elements[j]
        if b < c or d < a:
            case = "separate"
        elif (a <= c and d <= b) or (c <= a and b <= d):
            case = "nested"
        else:
            case = "partial overlap"
        print(f"T{i + 1} and T{j + 1}: {case}")
