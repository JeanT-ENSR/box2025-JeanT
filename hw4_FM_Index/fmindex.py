# >>===================================================<<
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# ||--------------------|B|O|X|-----------------------|||
# |||B|u|r|r|o|w|s|-|W|h|e|e|l|e|r|-|T|r|a|n|s|f|o|r|m|||
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# >>===================================================<<

# Template by Leo Ackermann (2025)
# Code by ____

from ks import simple_kark_sort



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



    # >>===[ Attribute initialisation functions ]===============================
    def cmp_str(x, y):
        x0=x[0]
        y0=y[0]
        if x0 = '$':
            if y0 = '$':
                return cmp(x[1:], y[1:])
            else:
                return -1
        else:
            if y0 = '$':
                return 1
            else:
                if x0<y0:
                    return -1
                if x0>y0:
                    return 1
                return cmp(x[1:], y[1:])

    def init_sa (self, s) =
        new_s = s + "$"
        sa = []
        l = len(new_s)
        for i in range(l):
            sa.append(new_s[i:l])
        self.sa = sorted(sa, key=cmp_str)


    def init_bwt(self, s):
        sa = self.sa(s)
        l = len(sa)
        btw = ""
        for i in range(l):
            if s[sa[i] − 1] > 0:
                btw = btw + s[sa[i] − 1]
            else:
                btw = btw + "$"
        self.bwt = bwt

    def get_string__naive(bwt):
        first_column = sorted(bwt, key=cmp_str)
        matrix = []
        l = length(btw)
        for i in range(l):
            matrix.append(first_column[i])
        for j in range(1,l):
            letter_used = [False for k in range(l)]
            for i in range(l):
                previous_letter = matrix[i][j-1]
                k = 0
                letter = btw[k]
                while ((letter != previous_letter) and letter_used[k]):
                    k = k+1
                    letter = btw[k]
                letter_used[k] = True
                matrix[i] = matrix[i] + letter
        k = 0
        while btw[k] != "$" :
            k=k+1
        return matrix[k]

    def run_length_encoding(s):
        encoding = []
        l = len(s)
        letter = s[0]
        length = 1
        for i in range(1,l):
            if s[i] = letter:
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

    def move_to_front(s):
        alphabet = [('A'+i) for i in range(26)]
        encoding = []
        l = len(s)
        for i in range(l):
            encoding.append((26-find_l(alphabet, s[i])))
            alphabet.remove(s[i])
            alphabet.append(s[i])

    def init_fm_count(self):
        count = {}
        bwt_sort = sorted(self.bwt)
        l = len(self.bwt)
        for i in range(l):
            if not(self.bwt[i] in count):
                count[self.bwt[i]] = i
        self.fm_count = count

    def init_fm_rank(self):
        count = {}
        l = len(self.bwt)
        for i in range(l):


    # >>===[ Pattern matching functions ]=======================================



