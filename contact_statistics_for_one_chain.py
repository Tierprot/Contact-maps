__author__ = 'Bones'
import matplotlib.pylab as plt
import numpy as np
import re

def process_line(string):
    temp_var = string.split()
    return temp_var[3]+temp_var[5]+'_'+temp_var[4],\
           [float(temp_var[6]),float(temp_var[7]),float(temp_var[8])]

def make_library(filname):
    residue = []
    coordinates = []
    with open(filname + ".pdb", 'r') as inputfile:
        for string in inputfile:
            if "ATOM" in string and\
               "CA" in string:
                    residue_curr, coordinates_curr = process_line(string)
                    residue.append(residue_curr)
                    coordinates.append(coordinates_curr)
    return residue, coordinates

def distance(a, b):
    return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5

def print_to_file(contacts,filename):
    with open(filename + ".txt", 'w') as output:
        for record in contacts:
            output.write("{0} {1} {2}\n".format(record[0],record[1], record[2]))

def filter_by_distance(contacts, min_distance = 3, max_distance = 5):
    return list(filter(lambda contact: min_distance <= contact[2] <= max_distance, contacts))

def filter_by_molecule(contacts, chain):
    return list(filter(lambda contact: contact[0][-1] == chain and contact[1][-1] == chain, contacts))

def filter_by_intermolecule(contacts, *chains):

    test = filter(lambda contact: True if contact[0][-1] in chains and contact[1][-1] in chains\
                                          and contact[0][-1] != contact[1][-1] else False, contacts)
    return list(test)


def print_contacts(contacts):
    for contact in contacts:
        print(contact)

if __name__ == '__main__':

    #inputname = input("Enter pdf filename:\n")
    #inputchain = input("Enter chain name (letter):\n")
    #min_distance = input("Enter minimum distance:\n")
    #max_distance = input("Enter maximum distance:\n")
    inputname = "model.000.00"
    inputchain = "A"
    min_distance = 3
    max_distance = 12

    #делаем поостаточную библиотеку
    residues, coordinates = make_library(inputname)

    #считаем контакты
    contacts = []

    for i in range(len(residues)):
        for token in range(i+1, len(residues)):
            delta = distance(coordinates[i], coordinates[token])
            if min_distance <= delta <= max_distance:
                contacts.append([residues[i], residues[token], delta])

    #сортировка по значениям
    contacts = sorted(contacts, key=lambda contacts: contacts[2])

    #пишем данные в файл
    print_to_file(contacts, inputname)

    #фильтр на втутримолекулярные и межмолекулярные остатки
    mol = filter_by_molecule(contacts, inputchain)

    #строим график
    #выделяем столбец с контактами

    a, b, c = zip(*mol)
    x = range(len(mol))
    cont = c[:]

    numseqA = []
    numseqB = []

    #подготовка к картированию контактов,
    #формируем числовые столбцы последовательностей
    pat = re.compile(r'[\d]+', re.I)
    for i in range(len(mol)):
        num1 = re.search(pat, mol[i][0])
        num2 = re.search(pat, mol[i][1])
        numseqA.append(int(num1.group(0)))
        numseqB.append(int(num2.group(0)))

    a, b, z = zip(*mol)


    plt.figure(figsize=(18, 6))
    plt.subplot(121)
    plt.hist(cont, bins=max_distance-min_distance, histtype='stepfilled', color="blue", range=(min_distance, max_distance),label="contacts per 1 A")
    plt.hist(cont, bins=(max_distance-min_distance)*2, histtype='stepfilled', color="red",alpha=0.5,label="contacts per 0.5 A", range=(min_distance, max_distance))
    plt.xlabel('distance between atoms')
    plt.ylabel('number of contacts')

    plt.subplot(122, axisbg='darkblue')
    gridsize = lambda f: max(f) if max(f) >= 100 else 100
    plt.hexbin(numseqA, numseqB, z, gridsize=gridsize(numseqA), cmap=plt.cm.jet_r, vmin=min_distance, vmax=max_distance)
    plt.axis([0, max(numseqA), 0, max(numseqB)])
    plt.xlabel('sequence')
    plt.ylabel('sequence')
    plt.title("Contacts map")
    cb = plt.colorbar()
    cb.set_label('distance, A')
    plt.show()
