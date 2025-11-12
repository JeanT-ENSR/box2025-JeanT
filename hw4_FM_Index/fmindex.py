# >>===================================================<<
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# ||--------------------|B|O|X|-----------------------|||
# |||B|u|r|r|o|w|s|-|W|h|e|e|l|e|r|-|T|r|a|n|s|f|o|r|m|||
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# >>===================================================<<

# Template by Leo Ackermann (2025)
# Code by ____

from ks import simple_kark_sort
import numpy as np


# >>==========================================================================<<
# >>                              FM-index class                              <<
# >>==========================================================================<<

class FMindex:

    # >>===[ Class constructor ]================================================

    def __init__(self, seq, verbose=False):
        self.sa = None
        self.bwt = None
        self.fm_count = None
        self.fm_rank = None
        self.fm_ranks = None
        self.next_smallest_letter = None

        self.__compute_sa(seq)
        self.__compute_bwt_from_sa(seq)
        self.__compute_fm_count()
        self.__compute_fm_rank()
        self.__compute_fm_ranks()



    # >>===[ Attribute initialisation functions ]===============================
    def cmp_str(x, y):
        '''
        Fonction de comparaison pour les strings qui prend en compte $ comme caractère le plus petit
        '''
        x0=x[0]
        y0=y[0]
        if x0 == '$':
            if y0 == '$':
                return cmp(x[1:], y[1:])
            else:
                return -1
        else:
            if y0 == '$':
                return 1
            else:
                if x0<y0:
                    return -1
                if x0>y0:
                    return 1
                return cmp(x[1:], y[1:])

    def __compute_sa(self, seq):
        '''
        Compute the suffix array of the sequence (with termination symbol), in
        linear time
        '''
        # NOTE. Relies on an external function

        self.sa = simple_kark_sort(seq + '$')


    def __compute_bwt_from_sa(self, seq):
        '''
        Compute the Burrows-Wheeler transform using the SA-based definition, in linear time
        NOTE : self.sa must be initialized
        '''

        self.bwt = ''.join(map(lambda x: seq[x-1] if x>0 else '$', self.sa))

    def get_string__naive(self):
        '''
        Retrieve the string encoded within the BWT, in time O(|S|^3 * log|S|)
        NOTE : self.bwt must be initialized
        '''

        first_column = sorted(self.bwt, key=cmp_str)
        first_cols = []
        l = length(self.bwt)

        # Initialisation de la première colonne de la matrice
        for i in range(l):
            first_cols.append(first_column[i])

        # Construction de la colonne suivante de la matrice à partir des colonnes précédentes
        for _ in range(1,l):
            for i in range(l):
                first_col[i] = self.bwt[i] + first_col[i]
            first_col = sorted(first_col, key=cmp_str)

        # On recherche le mot qui se termine par '$'
        k = 0
        while self.bwt[k] != "$" :
            k=k+1
        return first_cols[k]

    def run_length_encoding(s):
        encoding = []
        l = len(s)
        letter = s[0]
        length = 1
        for i in range(1,l):
            if s[i] == letter:
                length = length+1
            else:
                encoding.append(letter,length)
                length = 1
                letter = s[i]
        return encoding

    def find_l(alpha, l):
        i = 0
        letter = alpha[i]
        while letter != l:
            i=i+1
            letter = alpha[i]
        return i

    def move_to_front_encoding(s):
        alphabet = [('A'+i) for i in range(26)]
        encoding = []
        l = len(s)
        for i in range(l):
            encoding.append((26-find_l(alphabet, s[i])))
            alphabet.remove(s[i])
            alphabet.append(s[i])

    def __compute_fm_count(self):
        '''
        Compute the Cumulative count map
        NOTE : self.bwt must be initialized
        '''
        count = {}
        bwt_sort = sorted(self.bwt)
        l = len(self.bwt)
        for i in range(l):
            if not(bwt_sort[i] in count):
                count[bwt_sort[i]] = i
        self.fm_count = count

    def __compute_fm_rank(self):
        '''
        Compute the Rank array
        NOTE : self.bwt must be initialized
        '''
        rank = []
        l = len(self.bwt)
        for i in range(l):
            nb = 0
            for j in range(i):
                if self.bwt[i] == self.bwt[j]:
                    nb = nb + 1
            rank.append(nb+1)
        self.fm_rank = rank

    def __compute_fm_ranks(self):
        '''
        Compute the rank array but split by letters
        NOTE : self.bwt must be initialized
        '''
        l = len(self.bwt)
        ranks = {"$":np.zeros(l), "A":np.zeros(l), "C":np.zeros(l), "G":np.zeros(l), "T":np.zeros(l)}
        ranks[self.bwt[0]][0] = 1
        for i in range(1,l):
            for key in ranks:
                if key == self.bwt[i]:
                    ranks[key][i] = ranks[key][i-1] + 1
                else:
                    ranks[key][i] = ranks[key][i-1] + 1
        self.fm_ranks = ranks

    def __compute_next_smallest_letter(self):
        '''
        JE NE ME RAPELLE PLUS CE QUE C'EST
        '''
        letter = "A"
        self.next_smallest_letter = letter

    # >>===[ Pattern matching functions ]=======================================



## Test ##
fmi = FMindex("ACAACG")
sa = fmi.sa
print(sa)

