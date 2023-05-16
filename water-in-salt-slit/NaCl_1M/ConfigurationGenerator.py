#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import random
import copy
from numpy.linalg import norm

def neighborsearch(neighbor,molecule,cptatm, x, y, z, Lx, Ly, Lz):
    '''Search all neighbor to a molecule in a box and return the closest distance'''
    box = np.array([Lx, Ly, Lz])
    minr = 10
    for m in molecule.T:
        x0 = m[0] + x
        y0 = m[1] + y
        z0 = m[2] + z
        dxdydz = np.remainder(neighbor[:cptatm].T - np.array([x0,y0,z0]) + box/2., box) - box/2.
        minr = np.min([minr,np.min(norm(dxdydz,axis=1))])
    return minr

Na = 6.022e23 #constants.Avogadro
Mh2o = 0.018053 # kg/mol - water

N = 4000

# desired concentration in mol/L
c = 1
nion = c*N*Mh2o/(2*(1+Mh2o*c)) # desired number for each ion
nwater = N - 2*nion

# choose the initial box dimensions
dnacl = 2.84
nx = 14
ny = 14
nz = 4
dw = 3.1
layer = nz*dnacl
h = 5
Lx = nx*dnacl
Ly = ny*dnacl

cptH2O = 0
nCl = 0
nNa = 0

txlo, txhi = 0, Lx
tylo, tyhi = 0, Ly

attemps = 0

while cptH2O+nNa+nCl < N:
    
   
    Lz = layer + h + dw*attemps
    tzlo, tzhi = 0,Lz
    
    cptatom = 0
    cptbond = 0
    cptangle = 0
    cptmol = 0
    cptNa = 0
    cptH2O = 0
    cptCl = 0
    cptClf = 0
    cptNaf = 0
    cptres = 1
    nCl = 0
    nNa = 0
    
    box = np.array([Lx, Ly, Lz])

    # allocate memory
    XYZ = np.zeros((1000000,3))
    Typ = ["" for x in range(1000000)]
    ResName = ["" for x in range(1000000)]
    ResNum = np.zeros((1000000,1))

    # Load NaCl positions for the wall crystal structure
    wallNaCl = np.zeros((10000,7))
    file1 = open('NaCl/Position.dat', 'r')
    Lines = file1.readlines()
    count = 0
    for line in Lines:
        wallNaCl[count]=line.strip().split(' ')
        count += 1
    wallNaCl = wallNaCl[0:count]

    # replicate the initial structure
    wallNaClrep = copy.deepcopy(wallNaCl)
    for xx in np.arange(txlo+dnacl/2,txhi,2*dnacl):
        for yy in np.arange(tylo+dnacl/2,tyhi,2*dnacl):
            for zz in np.arange(tzlo+dnacl/2,tzlo+layer,2*dnacl):
                wallNaClrep = np.append(wallNaClrep,wallNaCl+[0,0,0,0,xx,yy,zz], axis=0)
    wallNaClrep = wallNaClrep[8:]
    assert len(wallNaClrep[wallNaClrep.T[2]==1]) == len(wallNaClrep[wallNaClrep.T[2]==2])

    for n in range(len(wallNaClrep)):
        wallNaClrep[n,0] = np.int64(n+1)

    shift = np.min(wallNaClrep.T[6])
    wallNaClrep.T[6] -= shift
    layer = np.max(wallNaClrep.T[6]) - np.min(wallNaClrep.T[6])

    wallCl = wallNaClrep[wallNaClrep.T[2] == 2]

    # choose the Cl to maintain fix all over the production run (just one atom)
    selection = wallCl[(wallCl.T[6]>layer*0.3) & (wallCl.T[6]<layer*0.7) & (wallCl.T[2] == 2)][1]
    assert len(selection)>0

    for n in np.arange(selection[0], selection[0]+1,1):
        id0 = np.where(wallNaClrep.T[0] == n)
        wallNaClrep[id0,2] += 1

    wallCl = []
    wallCl = wallNaClrep[wallNaClrep.T[2] > 1]
    wallNa = wallNaClrep[wallNaClrep.T[2] == 1]  

    # place the Cl of the wall
    for m in wallCl:
        if m[2] == 2:
            x0 = m[4]
            y0 = m[5]
            z0 = m[6]
            XYZ[cptatom] = [x0,y0,z0]
            Typ[cptatom] = 'Cl'
            ResNum[cptatom] = cptres
            ResName[cptatom] = 'Cl'
            cptCl += 1
            cptres += 1
            cptatom += 1

    # add Cl randomly
    fail = 0
    while (nCl < nion) & (fail < 1e4):
        x = random.randint(1,1000)/1000*(txhi-txlo)+txlo
        y = random.randint(1,1000)/1000*(tyhi-tylo)+tylo
        z = random.randint(1,1000)/1000*(tzhi-tzlo)+tzlo
        pos = np.array([x,y,z])
        
        dxdydz = np.remainder((pos - XYZ[0:cptatom]) + box/2., box) - box/2.
        d = np.min(norm(dxdydz,axis=1))

        if d > 3 and z>tzlo+layer+3 and z<tzhi-3:
            fail = 0
            XYZ[cptatom] = [x,y,z]
            Typ[cptatom] = 'Cl'
            ResNum[cptatom] = cptres
            ResName[cptatom] = 'Cl'
            cptatom += 1  
            nCl += 1
            cptCl += 1
            cptres += 1
        else:
            fail += 1

    # place the Cl fixed
    for m in wallCl:
        if m[2]==3:
            x0 = m[4]
            y0 = m[5]
            z0 = m[6]
            XYZ[cptatom] = [x0,y0,z0]
            Typ[cptatom] = 'Cl'
            ResNum[cptatom] = cptres
            ResName[cptatom] = 'Clf'          
            cptClf += 1
            cptres += 1
            cptatom += 1  

    # place the Na
    for m in wallNa:
        x0 = m[4]
        y0 = m[5]
        z0 = m[6]
        XYZ[cptatom] = [x0,y0,z0]
        Typ[cptatom] = 'Na'
        ResNum[cptatom] = cptres
        ResName[cptatom] = 'Naf'
        cptatom += 1
        cptNaf += 1
        cptres += 1

    # add Na randomly
    fail = 0
    while (nNa < nion) & (fail < 1e4):
        x = random.randint(1,1000)/1000*(txhi-txlo)+txlo
        y = random.randint(1,1000)/1000*(tyhi-tylo)+tylo
        z = random.randint(1,1000)/1000*(tzhi-tzlo)+tzlo

        pos = np.array([x,y,z])
        
        dxdydz = np.remainder((pos - XYZ[0:cptatom]) + box/2., box) - box/2.
        d = np.min(norm(dxdydz,axis=1))
         
        if d > 3 and z>tzlo+layer+1.5 and z<tzhi-1.5:
            fail = 0
            XYZ[cptatom] = [x,y,z]
            Typ[cptatom] = 'Na'
            ResNum[cptatom] = cptres
            ResName[cptatom] = 'Na'
            cptatom += 1   
            nNa += 1
            cptres += 1
            cptNa += 1
        else:
            fail += 1

    cptH2O = 0    
    # create water    
    rOH = 0.9572 
    rOM = 0.105 # tip4p/epsilon
    thetaHOH = 104.52    
    PosH2O = np.array([[0, 0, 0],        [rOH*np.cos((thetaHOH/2)*np.pi/180),   rOH*np.sin((thetaHOH/2)*np.pi/180), 0.0],        [rOH*np.cos((thetaHOH/2)*np.pi/180),   -rOH*np.sin((thetaHOH/2)*np.pi/180),  0.0],        [rOM,  0.0, 0.0]])
    
    cption = cptatom

    TypH2O = ['OW', 'HW1', 'HW2', 'MW']
    for x in np.arange(txlo+dw/2,txhi,dw):
        for y in np.arange(tylo+dw/2,tyhi,dw):
            for z in np.arange(tzlo+dw/2,tzhi,dw):

                pos = np.array([x,y,z])
                dxdydz = np.remainder((pos - XYZ[0:cption]) + box/2., box) - box/2.
                d = np.min(norm(dxdydz,axis=1))
                
                if d > dw and z>tzlo+layer+3 and z<tzhi-3 and cptH2O < nwater:
                    for j in range(4):
                        XYZ[cptatom] = [x,y,z]+np.array(PosH2O[j])
                        Typ[cptatom] = TypH2O[j]
                        ResNum[cptatom] = cptres
                        ResName[cptatom] = 'SOL'
                        cptatom += 1    
                    cptH2O += 1
                    cptres += 1

    attemps += 1
    print('attemp number : '+str(attemps)+', Lz = '+str(Lz) + ', molecules + ions = '+str(cptH2O+nNa+nCl))
         
print('Lx = '+str(Lx/10)+' nm, Ly = '+str(Ly/10)+' nm, Lz = '+str(Lz/10)+' nm, hinit = '+str(h/10)+' nm')
print(str(nNa)+' Na ions in electrolyte') 
print(str(nCl)+' Cl ions in electrolyte')
print(str(cptH2O)+' water molecules')

Vwater = cptH2O/6.022e23*0.018 # kg or litter
Naddion = (nCl+nNa)/6.022e23 # mol
cion = Naddion/Vwater
print('The initial ion concentration is '+str(cion)+' M')

# write conf.gro
f = open('conf.gro', 'w')
f.write('SO4Na2 slit\n')
f.write(str(cptatom)+'\n')
for n in range(cptatom):
    f.write("{: >5}".format(str(np.int32(ResNum[n][0])))) # residue number (5 positions, integer) 
    f.write("{: >5}".format(str(ResName[n]))) # residue name (5 characters) 
    f.write("{: >5}".format(str(Typ[n]))) # atom name (5 characters) 
    f.write("{: >5}".format(str(np.int32(n+1)))) # atom number (5 positions, integer)
    f.write("{: >8}".format(str("{:.3f}".format(XYZ[n][0]/10)))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    f.write("{: >8}".format(str("{:.3f}".format(XYZ[n][1]/10)))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
    f.write("{: >8}".format(str("{:.3f}".format(XYZ[n][2]/10)))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
    f.write("\n")
f.write("{: >10}".format(str("{:.5f}".format(Lx/10))))
f.write("{: >10}".format(str("{:.5f}".format(Ly/10))))
f.write("{: >10}".format(str("{:.5f}".format(Lz/10))))
f.close()
# write topol.top
f = open('topol.top', 'w')
f.write('#include "../ff/forcefield.itp"\n')
f.write('#include "../ff/tip4peps.itp"\n')
f.write('#include "../ff/ions.itp"\n\n')
f.write('[ System ]\n')
f.write('SO4Na2 slit\n\n')
f.write('[ Molecules ]\n\n')
f.write('Cl '+ str(cptCl)+'\n')
f.write('Clf '+ str(cptClf)+'\n')
f.write('Naf '+ str(cptNaf)+'\n')
f.write('Na '+ str(cptNa)+'\n')
f.write('SOL '+ str(cptH2O)+'\n')
f.close()


# In[ ]:




