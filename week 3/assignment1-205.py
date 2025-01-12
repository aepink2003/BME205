

#Questionsn 4
a = 4422
b= 8548
s=0

for i in range(a,b+1):
    if i%2 == 1:
        s += i
        
print(s)
    
    #Questionsn 5
import this
i = 1
f = open('rosalind_ini6.txt')
for line in f.readlines():
    if i % 2 == 0 :
        print line
    i += 1
    
#QUESTION 6
s = "When I find myself in times of trouble Mother Mary comes to me Speaking words of wisdom let it be And in my hour of darkness she is standing right in front of me Speaking words of wisdom let it be Let it be let it be let it be let it be Whisper words of wisdom let it be And when the broken hearted people living in the world agree There will be an answer let it be For though they may be parted there is still a chance that they will see There will be an answer let it be Let it be let it be let it be let it be There will be an answer let it be Let it be let it be let it be let it be Whisper words of wisdom let it be Let it be let it be let it be let it be Whisper words of wisdom let it be And when the night is cloudy there is still a light that shines on me Shine until tomorrow let it be I wake up to the sound of music Mother Mary comes to me Speaking words of wisdom let it be Let it be let it be let it be yeah let it be There will be an answer let it be Let it be let it be let it be yeah let it be Whisper words of wisdom let it be"
d = {}
for i in s.split(" "):
    if i in d:
        d[i] +=1
    else:
        d[i] = 1
        

#QUESTION 7
s = "TTGGCAATCCACGCAAGTACATACACCAGAGGCGACGGCACTCGTCTATTTCTATTTTCGCTCACAACCCCCTATAGCTGTTCGCGCATGACCGCATGGAAACGAATGGAGACAGCCTTGACCGTACGTATCGAGTGGAGGTGATATTCCGGAGTTTAGGCTTCACTCCATAACCCTGATGACCTCAGATTCCATCTCGACTGTCAAACAAATCGACTTAGACCCTTGTTCCCATAAATACGATGAAAGCGACTGCCATAGCCATGTAGTCCACATCGATAGATAGGGGCGCACCTCCCCGCTCCTACTCCGCACACTAGTAGGGTGAGTGATTGTGCTATATAATAGCGTCACGGAAAATAAGCCAGAGTCAATCGCCGATGTCTCGTTCAAACCCATGACTACTGTGTCATCACCCGACTGTTAAAACATATTCTCGAAGCTTGAAAGGTCTGTCAGGCATCGTTGAGAACCAGTCCCCTAATCTCCCGGTCGTGAGGTGCACTATCAGACATGTCCCGATTTAAGTGCGGCTCACAGAAGGAGAGCCGGTTACTGACGACCTAATGGCTTAAGTCCCGGTCTATTACTTAAATTACCAGAGATCGCGCGATTACTCAAAAACGATTCCCCTTCATATCGAGGCTTGTCCTAGAAACTGTAGAAGGTCTCAATATTCCGATTCATTGTCGACCATTGACTGGTTAATCGCGTTTAGCAATCCAAATGTGCACGGGCCCAGTTATATGCCAGGAAGGGAAAATATTCATATATCCACTTTACAGCGCAAGAGGGGATACTAGCATGAGCTTCTAGGTTTAAATATGACGACAGACGCCCATATTGGACCGTCTGGCGGTGATAATTTGT"
d = {"A":0,"C":0,"G":0,"T":0}

for char in s:
    if char in d:
        d[char] +=1

for c in d.items():
    print(c)
        