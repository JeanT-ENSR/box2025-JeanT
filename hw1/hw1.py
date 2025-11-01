def parser(filename):
    f = open(filename, "r")
    dico = []
    w_buffer = ""

    for line in f :
        if line[0] == '>':
            dico.append(w_buffer)
            w_buffer = ""
        else:
            w_buffer = ''.join((w_buffer, line[0:(len(line))-1]))
    dico.append(w_buffer)

    f.close()
    return dico

sequences = parser("genome_hw1.fa")

def F_one_word(k, w):
    '''
    F(k) is the number of distinct subwords of length k in w
    '''
    compteur = 0
    dict = {}
    for i in range(len(w) - k + 1):
        sw = w[i:i+k]
        compteur += not(sw in dict)
        dict[sw] = 1
    return (compteur,dict)

def F_many_words(k, words):
    compteur, dict = F_one_word(k, words[0])
    for w in words:
        for i in range(len(w) - k + 1):
            sw = w[i:i+k]
            compteur += not(sw in dict)
            dict[sw] = 1
    return compteur


F_sequences = []
for k in range(2,31):
    F_sequences.append(F_many_words(k,sequences))
print(F_sequences)

def Fibonacci_w (l):
    if (l < 1):
        return ""
    elif (l == 1):
        return "A"
    elif (l == 2):
        return "AB"
    else:
        n0 = 1
        n = 2
        w = ("A", "AB")
        word = "AB"
        while (n < l):
            word = ''.join(w)
            w = (w[1], word)
            n += n0
            n0 = n - n0
        return word

F_fibonacci= []
l_sequences = 0
for s in sequences:
    l_sequences = max(l_sequences, len(s))

for k in range(2,31):
    F_fibonacci.append(F_one_word(k,Fibonacci_w(l_sequences)))
print(F_fibonacci)
