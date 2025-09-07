def F(k, w):
    '''
    F(k) is the number of distinct subwords of length k in w
    '''
    compteur = 0
    dict = {}
    for i in range(len(w) - k + 1):
        sw = w[i:i+k]
        compteur += not(sw in dict)
        dict[sw] = 1
    return compteur

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

print(Fibonacci_w(1))
print(Fibonacci_w(2))
print(Fibonacci_w(3))
print(Fibonacci_w(5))
print(Fibonacci_w(8))
print(Fibonacci_w(30))