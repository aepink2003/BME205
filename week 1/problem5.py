    #Questionsn 5
import this
i = 1
f = open('rosalind_ini6.txt')
for line in f.readlines():
    if i % 2 == 0 :
        print line
    i += 1