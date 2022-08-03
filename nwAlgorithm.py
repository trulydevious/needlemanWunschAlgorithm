import numpy as np
import matplotlib.pyplot as plt 
import sys
import os
import csv


# Program do ilościowego porównania par sekwencji za pomocą algorytmu dopasowania globalnego Needlemana-Wunscha

def readFastaFile(file):
    
    """
    Funkcja wczytująca plik FASTA
    Args:
        file (str): nazwa pliku FASTA
    Returns:
        sequence (str): zczytana sekwencja DNA
    Raise:
        Exception("Empty file")
        Exception("Inorrect file format")
    """
    
    fastaFile = open(file)
    sequence = ''
    
    if file.endswith(('.fasta', '.txt')):                   # Sprawdzenie rozszerzenia pliku

        
        if os.stat(file).st_size == 0:                      # Sprawdzenie, czy plik jest pusty, czy nie
            raise Exception("Empty file")
    
        else:
            for line in fastaFile:                          # Wczytanie pliku z pominięciem pierwszej linii
                if line[0] != '>':
                    sequence = sequence + line.strip()
    else:
        raise Exception("Inorrect file format")
    
    fastaFile.close()
    return sequence


def createScoreMatrix(seq1, seq2):
    
    """
    Funkcja tworząca macierz z punktami, czy sekwencje mają te same nukleotydy
    Args:
        seq1 (str): pierwsza sekwencja DNA
        seq2 (str): druga sekwencja DNA
    Returns:
        matrix: macierz z punktami
    """       
            
    scoreMatrix = np.zeros((len(seq1), len(seq2)))          # Wygenerowanie macierzy wypełnionej zerami
        
    for i in range(len(seq1)):                              # Wypełnienie macierzy punktami, w zależności czy nukleotydy się zgadzają
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                scoreMatrix[i][j]= match
            else:
                scoreMatrix[i][j]= mismatch
    
    return scoreMatrix


def smithWaterman(seq1, seq2, scoreMatrix):
    
    """
    Funkcja tworząca macierz na podstawie sekwencji DNA
    Args:
        seq1 (str): pierwsza sekwencja DNA
        seq2 (str): druga sekwencja DNA
    Returns:
        matrix: macierz z punktami wg algorytmu NW
    """       
    
    right_arrow = "\u2192"
    down_arrow = "\u2193"
    down_right_arrow = "\u2198"
                
    matrix = np.zeros((len(seq1)+1, len(seq2)+1))   # Generowanie macierzy
    pathMatrix = np.zeros((len(seq1)+1, len(seq2)+1), dtype = str)  
        
    for i in range(len(seq1)+1):                    # Wypełnienie pierwszej kolumny i wiersza
        matrix[i][0] = i * 0
        pathMatrix[i][0] = right_arrow
        
    for j in range(len(seq2)+1):
        matrix[0][j] = j * 0
        pathMatrix[0][j] = down_arrow
                
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            diagonal = matrix[i - 1][j - 1] + scoreMatrix[i-1][j-1]
            vertical = matrix[i - 1][j] + gap
            horizontal = matrix[i][j - 1] + gap
            matrix[i][j] = max(diagonal, horizontal, vertical, 0)
            if max(diagonal, horizontal, vertical) == diagonal:
                pathMatrix[i][j] = down_right_arrow
            elif max(diagonal, horizontal, vertical) == horizontal: 
                pathMatrix[i][j] = right_arrow
            else:
                pathMatrix[i][j] = down_arrow
                
    return matrix, pathMatrix


def needlemanWunsch(seq1, seq2, scoreMatrix):
    
    """
    Funkcja tworząca macierz na podstawie sekwencji DNA
    Args:
        seq1 (str): pierwsza sekwencja DNA
        seq2 (str): druga sekwencja DNA
    Returns:
        matrix: macierz z punktami wg algorytmu NW
        pathMatrix: macierz z przedstwioną scieżką
    """       
    
    right_arrow = "\u2192"
    down_arrow = "\u2193"
    down_right_arrow = "\u2198"
                
    matrix = np.zeros((len(seq1)+1, len(seq2)+1))               # Wygenerowanie macierzy wypełnionej zerami
    pathMatrix = np.zeros((len(seq1)+1, len(seq2)+1), dtype = str)    
    
    for i in range(len(seq1)+1):                                # Wypełnienie pierwszej kolumny i wiersza
        matrix[i][0] = i * gap
        pathMatrix[i][0] = right_arrow
    
    for j in range(len(seq2)+1):                                # Wypełnienie środka macierzy
        matrix[0][j] = j * gap
        pathMatrix[0][j] = down_arrow
                
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            diagonal = matrix[i - 1][j - 1] + scoreMatrix[i-1][j-1]         # Po skosie - zgodnie z punktami we wczesniej utworzonej macierzy
            vertical = matrix[i - 1][j] + gap                               # Pionowo - dodając wartość kary za przerwę
            horizontal = matrix[i][j - 1] + gap                             # Poziomo - dodając wartość kary za przerwę 
            matrix[i][j] = max(diagonal, horizontal, vertical)              # Wybranie maksymalnej wartości z trzech możliwych
            if max(diagonal, horizontal, vertical) == diagonal:
                pathMatrix[i][j] = down_right_arrow
            elif max(diagonal, horizontal, vertical) == horizontal: 
                pathMatrix[i][j] = right_arrow
            else:
                pathMatrix[i][j] = down_arrow
                
    return matrix, pathMatrix


def globalAlligment(seq1, seq2, matrix, scoreMatrix):
    
    """
    Funkcja globalne generująca dopasowanie i ścieżkę
    Args:
        seq1 (str): pierwsza sekwencja DNA
        seq2 (str): druga sekwencja DNA
        matrix: macierz z punktami wg algorytmu NW
        scoreMatrix: macierz z punktami
    Returns:
        alignment1 (str): pierwsza sekwencja DNA po dopasowaniu
        alignment2 (str): druga sekwencja DNA po dopasowaniu
    """    
    
    matrix, pathMatrix = needlemanWunsch(seq1, seq2, createScoreMatrix(seq1, seq2))
    alignment1 = ''                                             # Zdefiniowanie danych początkowych                                        
    alignment2 = ''
    i = len(seq1)
    j = len(seq2)

    while (i > 0 and j > 0):
  
        score = matrix[i][j]                                    # Wyniki w zależności od położenia (po skosie, pionowo, poziomo)
        scoreDiag = matrix[i-1][j-1]
        scoreVertic = matrix[i-1][j]
        scoreHoriz = matrix[i][j-1]
    
        if score == scoreDiag + scoreMatrix[i-1][j-1]:          # Jeśli aktualny wynik, to wynik otrzymany po przejściu po skosie,
            alignment2 += seq2[j-1]                             # to pierwsza i druga sekwencja zostają tak samo
            alignment1 += seq1[i-1]
            i -= 1
            j -= 1
        
        elif score == scoreVertic + gap:                        # Jeśli aktualny wynik, to wynik otrzymany po przejściu pionowo,
            alignment2 += '-'                                   # to dodaj "-" do drugiej sekwencji, a pierwsza zostaje tak samo
            alignment1 += seq1[i-1]
            i -= 1
        
        elif score == scoreHoriz + gap:                         # Jeśli aktualny wynik, to wynik otrzymany po przejściu poziomo,
            alignment1 += '-'                                   # to dodaj "-" do pierwszej sekwencji, a druga zostaje tak samo
            alignment2 += seq2[j-1]
            j -= 1

    alignment1 = alignment1[::-1]                               # Odwrócenie wartości, ponieważ algorytm "szedł od tyłu"
    alignment2 = alignment2[::-1]
    
    return alignment1, alignment2


def localAlligment(seq1, seq2, matrix, scoreMatrix):
    
    """
    Funkcja globalne generująca dopasowanie i ścieżkę
    Args:
        seq1 (str): pierwsza sekwencja DNA
        seq2 (str): druga sekwencja DNA
        matrix: macierz z punktami wg algorytmu NW
        scoreMatrix: macierz z punktami
    Returns:
        alignment1 (str): pierwsza sekwencja DNA po dopasowaniu
        alignment2 (str): druga sekwencja DNA po dopasowaniu
    """   
     
    matrix, pathMatrix = smithWaterman(seq1, seq2, createScoreMatrix(seq1, seq2))
    alignment1 = ''
    alignment2 = ''
    i = len(seq1)
    j = len(seq2)

    while (i > 0 and j > 0):
  
        score = matrix[i][j]
        scoreDiag = matrix[i-1][j-1]
        scoreVertic = matrix[i-1][j]
        scoreHoriz = matrix[i][j-1]
        
        if score != 0:
    
            if score == scoreDiag + scoreMatrix[i-1][j-1]:
                alignment1 += seq1[i-1]
                alignment2 += seq2[j-1]
                i -= 1
                j -= 1
            
            elif score == scoreVertic + gap:
                alignment1 += seq1[i-1]
                alignment2 += ' '
                i -= 1
            
            elif score == scoreHoriz + gap:
                alignment1 += ' '
                alignment2 += seq2[j-1]
                j -= 1
        
        else:
            break

    alignment1 = alignment1[::-1]
    alignment2 = alignment2[::-1]
    
    return(alignment1, alignment2)


def makeStat(seq1, seq2):
    
    """
    Funkcja generująca statystyki
    Args:
        seq1 (str): pierwsza sekwencja DNA
        seq2 (str): druga sekwencja DNA
    Returns:
        identicalPercent: procent identycznych pozycji w sekwencjach DNA
        gapsPercent: procent przerw w sekwencjach DNA
    """   
    
    align1, align2 = globalAlligment(seq1, seq2, needlemanWunsch(seq1, seq2, createScoreMatrix(seq1, seq2)), createScoreMatrix(seq1, seq2))
    identical = 0
    gaps = 0
    
    for i in range(len(align1)):                                               # Przejście przez sekwencję 1 i 2 oraz zliczenie przerw oraz tych samych pozycji
        
        if align1[i] == align2[i]:
            identical += 1
        
        elif align1[i] == '-':
            gaps += 1
        
        elif align2[i] == '-':
            gaps += 1
            
    identicalPercent = round(identical * 100 / len(align1), 2)                 # Obliczenie procentu identycznych pozycji
    gapsPercent = round(gaps * 100 / len(align1), 2)                           # Obliczenie procentu przerw
            
    return(identicalPercent, gapsPercent)


def createTextToSave(seq1, seq2):
    
    """
    Funkcja generująca tekst do zapisu
    Args:
        seq1 (str): pierwsza sekwencja DNA
        seq2 (str): druga sekwencja DNA
    Returns:
        text: tablica z tekstem, który ma być zapisany
    """  
    
    matrixToSave, pathMatrix = needlemanWunsch(seq1, seq2, createScoreMatrix(seq1, seq2))
    matrixToSave2, pathMatrix2 = smithWaterman(seq1, seq2, createScoreMatrix(seq1, seq2))
    align1, align2 = globalAlligment(seq1, seq2, needlemanWunsch(seq1, seq2, createScoreMatrix(seq1, seq2)), createScoreMatrix(seq1, seq2))
    align11, align22 = localAlligment(seq1, seq2, smithWaterman(seq1, seq2, createScoreMatrix(seq1, seq2)), createScoreMatrix(seq1, seq2))
    identicalPercent, gapsPercent = makeStat(seq1, seq2)

    text = ["Sequence 1: " + seq1, "Sequence 2: " + seq2, "", "Scoring system: ", "Match: " + str(match), "Mismatch: " + str(mismatch), 
            "Gap: " + str(gap), "", "Needleman-Wunsch algorithm: ", matrixToSave, "", pathMatrix,  "", "Global allignment: ", align1, align2, "", 
            "Smith-Waterman algorithm: ", matrixToSave2,  "", pathMatrix2,  "", "Local allignment: ", align11, align22, "", 
            "Statistics: ", "Identical positions: " + str(identicalPercent) + " %", "Gaps: " + str(gapsPercent) + " %"]
    
    return text


def saveToFile(text):
    
    """
    Funkcja zapisująca dane do pliku .txt
    """
      
    with open("result.txt", "w", encoding='utf-8') as file:                     # Zapis do pliku .txt od nowych linii
        cw = csv.writer(file, delimiter="\n")
        cw.writerow(text)
        


# Podanie punktacji
print(os.listdir())
match = int(input("Podaj punktację dla zgodności:\n"))

if match < -10 or match > 11 :
    print("Podano nieprawidłowe dane.")
    sys.exit()

mismatch = int(input("Podaj punktację dla niezgodności:\n"))

if mismatch < -10 or mismatch > 11 :
    print("Podano nieprawidłowe dane.")
    sys.exit()

gap = int(input("Podaj punktację dla przerwy:\n"))

if gap < -10 or gap > 11 :
    print("Podano nieprawidłowe dane.")
    sys.exit()

# Wprowadzenie i przetworzenie danych

method = int(input("Wybierz nr metody:\n 1 - wprowadzenie danych w konsoli\n 2 - wczytanie pliku typu FATSA\n"))

if method == 1:
    
    dna1 = str(input("Podaj pierwszą sekwencję DNA:\n"))          # Wprowadzenie pierwszej sekwencji i sprawdzenie poprawności danych
    for letter in dna1:
        if letter not in 'ATCG':
            print("Podano nieprawidłowe dane.")
            sys.exit()

        
    dna2 = str(input("Podaj drugą sekwencję DNA:\n"))             # Wprowadzenie drugiej sekwencji i sprawdzenie poprawności danych
    for letter in dna2:
        if letter not in 'ATCG':
            print("Podano nieprawidłowe dane.")
            sys.exit()
            
    saveToFile(createTextToSave(dna1, dna2))
                
elif method == 2:
    
    file1 = str(input("Podaj nazwę pierwszego pliku:\n"))           # Wprowadzenie i wczytanie plików
    dna1 = readFastaFile(file1)
    
    file2 = str(input("Podaj nazwę drugiego pliku:\n"))
    dna2 = readFastaFile(file2)
    
    saveToFile(createTextToSave(dna1, dna2))
    
else: 
    print("Podano nieprawidłowe dane.")
    sys.exit()