import readfa

nb_seq = []

files = []
files.append((open('file1.fa', 'r')))
files.append((open('file2.fa', 'r')))
files.append((open('file3.fa', 'r')))
files.append((open('file4.fa', 'r')))
files.append((open('file5.fa', 'r')))
files.append((open('file6.fa', 'r')))

reads = []
reads_save = []
for i in range(6):
    reads.append((readfa.readfq(files[i])))
    reads_save.append(list ( reads[i] ))

print("num_seqs = ")
for i in range(6):
    print( len ( reads_save[i] ))


print("cumulative seq lengths")
for i in range(6):
    sum = 0
    for j in range(len(reads_save[i])):
        sum += len ((reads_save[i])[j][1])
    print( sum )












