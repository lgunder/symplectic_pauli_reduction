import random
import numpy
import math
import os

# computes comm of two paulis
def comm(p1, p2):
    comm = 0
    n = int(len(p1) / 2)
    for i in range(n):
        comm = comm + (p1[i] * p2[i + n] - p1[i + n] * p2[i])
    return comm % 2

# returns the full symp prod matrix
def sympmatrix(paulis):
    m = numpy.zeros([len(paulis), len(paulis)])
    for i in range(len(m)):
        for j in range(len(m)):
            m[i, j] = int(comm(paulis[i], paulis[j]))
    return m

##returns v, non-identical paulis on n registers
def randompaulis(v,n):
    paulis=[]
    while(len(paulis)<v):
        pauli=[]
        for i in range(2*n):
            pauli.append(random.randrange(0,2))
        if(not(pauli in paulis)):
            paulis.append(pauli)
    return paulis

def symptop(psymp):
    pauli=""
    n=int(len(psymp)/2)
    for i in range(n):
        if(psymp[i]==0 and psymp[i+n]==0):
            pauli=pauli+"I"
        elif(psymp[i]==1 and psymp[i+n]==0):
            pauli=pauli+"X"
        elif(psymp[i]==0 and psymp[i+n]==1):
            pauli=pauli+"Z"
        else:
            pauli=pauli+"Y"
    return pauli

def ptosymp(pstr):
    symp=numpy.zeros(2*len(pstr))
    for i in range(len(pstr)):
        if(pstr[i]=='X'):
            symp[i]=1
        if(pstr[i]=='Z'):
            symp[i+len(pstr)]=1
        if(pstr[i]=='Y'):
            symp[i]=1
            symp[i+len(pstr)]=1
    symp=symp.astype(int)
    return symp.tolist()

def pstosymp(pstrs):
    paulis=[]
    for pauli in pstrs:
        paulis.append(ptosymp(pauli))
    return paulis

#note that iso is at end and this is just one option, of course
def construct_min_rep(sympmatr):
    Lrank=symm_gaussian_elimination(sympmatr)
    L=Lrank[0]
    rank=int(Lrank[1])
    dim=int(len(sympmatr))
    regs=int(dim-rank/2)
    fund_p=numpy.zeros((dim,2*regs))
    for i in range(dim-rank):
        fund_p[i+2*int(rank/2)][i+regs+int(rank/2)]=1
    for i in range(int(rank/2)):
        fund_p[2*i][i]=1
        fund_p[2*i+1][i+regs]=1
    Linv=numpy.linalg.inv(L)%2
    minimal_rep=numpy.zeros((dim,2*regs))
    for i in range(dim):
        for j in range(dim):
            minimal_rep[i]=minimal_rep[i]+(Linv[i,j]*fund_p[j])%2
    pauli_forms=[]
    for i in range(dim):
        pauli_forms.append(symptop(minimal_rep[i]))
###note that this commented out line will print the minimal representations with the isotropic registers still attached
#    print(pauli_forms)
    return pauli_forms

def symm_gaussian_elimination(sympmatr):
    initialmat=numpy.copy(sympmatr)
    pos=0 #position number (up to len()-1)
    L=numpy.eye(len(sympmatr))
    for k in range(len(sympmatr)-1):
        pos=k
        pivot=-1
        for j in range(pos,len(sympmatr)):
            if(int(sympmatr[k,j])==1):
                pivot=j
                break
        #whole row is 0, so skip to next pos
        if(not(pivot==-1)):
            sympmatr=s_swap(sympmatr,pos+1,pivot)
            L[[pos+1,pivot]]=L[[pivot,pos+1]]
        else:
            continue
        for i in range(pos,len(sympmatr)):
            if(i==(pos+1)):
                continue
            if(int(sympmatr[k,i])==1):
                s_add(sympmatr,pos+1,i)
                L[i,:]=(L[i,:]+L[pos+1,:])%2
        for i in range(pos+1,len(sympmatr)):
            if((int(sympmatr[i,pos+1])==1) and not(i==k)):
                s_add(sympmatr,pos,i)
                L[i,:]=(L[i,:]+L[pos,:])%2
    print("half-rank (minimal direct sum registers) is:")
    rank=0
    for row in sympmatr:
        for entry in row:
            if(int(entry)==1):
                rank+=1
    print(rank/2)
    return [L,rank]

def s_swap(sympmatr,pos1,pos2):
    temp=numpy.copy(sympmatr[pos1,:])
    sympmatr[pos1,:]=sympmatr[pos2,:]
    sympmatr[pos2,:]=temp
    temp=numpy.copy(sympmatr[:,pos1])
    sympmatr[:,pos1]=sympmatr[:,pos2]
    sympmatr[:,pos2]=temp
    return sympmatr

def s_add(sympmatr,pos1,pos2):
    sympmatr[pos2,:]=(sympmatr[pos2,:]+sympmatr[pos1,:])%2
    sympmatr[:,pos2]=(sympmatr[:,pos2]+sympmatr[:,pos1])%2
    return sympmatr

##takes a single file only, a folder will be done afterwards
##assumes the file is in the form of whatever then Pauli operators, each on a different line. Weights are tossed out by this, but each Pauli will correspond to the original one in the given file
def readin(filename):
    f = open(filename,"r")
    lines = f.readlines()
    paulis=[]
    for line in lines:
        pauliop=""
        for char in line:
            if(char=="I" or char=="X" or char=="Y" or char=="Z"):
                pauliop+=char
        paulis.append(pauliop)
    print("number of operators:")
    print(len(paulis))
    print("original length:")
    print(len(paulis[0]))
    return paulis

##note that this will tack on each file name within the folder
def fullfolder(folderpath):
    dir_list = os.listdir(folderpath)
    for file in dir_list:
        print(file)
        onefile=readin(folderpath+"//"+file)
        construct_min_rep(sympmatrix(pstosymp(onefile)))
    return
     
path = ""#put your path to here

#this runs a full folder giving the data
fullfolder(path+"//Symp_pauli_reduction//example_molecules")
