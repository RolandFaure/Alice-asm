#python

import sys

def main():
    seq = sys.argv[1]
    seq = seq.upper()
    #reverse complement
    seq = seq[::-1]
    seq = seq.replace('A','t')
    seq = seq.replace('T','a')
    seq = seq.replace('G','c')
    seq = seq.replace('C','g')
    seq = seq.replace('<', 'w')
    seq = seq.replace('>', '<')
    seq = seq.replace('w', '>')
    seq = seq.upper()
    print(seq)

if __name__ == '__main__':
    main()