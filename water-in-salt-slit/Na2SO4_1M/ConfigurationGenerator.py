#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import random
import copy


# In[2]:


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


# In[6]:


Na = 6.022e23 #constants.Avogadro
Mh2o = 0.018053 # kg/mol - water
# desired number of water molecule
Ntot = 4000
#Ntot *= (3.99721*4.53025) / (3.97600*3.97600)
if Ntot < 800:
    double = True
else:
    double = False


# ### choose the initial box dimention

# In[7]:


lx = 5.7103 
ly = 9.0605
lz = 7.0275
nx = 7
ny = 5
nz = 4 
if double == True:
    nz *= 2
dw = 3.1
layer = nz*5.7103


# In[8]:


# desired concentration in mol/L
c = 1
nion = c*Ntot*Mh2o/(3*(1+Mh2o*c)) # desired number for each ion
nwater = Ntot - 3*nion
print('Number of ions '+str(np.round(nion,3)))


# ### initial NaCl mesh

# In[9]:


attemp = 0
cptH2O = 0
nNa = 0
nSO4 = 0
while cptH2O+nNa+nSO4 < Ntot:
    h = 60+attemp*3
    Lx = nx*lx
    Ly = ny*ly
    Lz = layer + h
    print('Lx = '+str(np.round(Lx/10,3))+' nm, Ly = '+str(np.round(Ly/10,2))+' nm, Lz = '+str(np.round(Lz/10,3))+' nm, hinit = '+str(np.round(h/10,3))+' nm')
    
    txlo = 0
    txhi = Lx
    tylo = 0 
    tyhi = Ly
    tzlo = 0 
    tzhi = Lz
    cptatom = 0
    cptbond = 0
    cptangle = 0
    cptmol = 0
    cptNa = 0
    cptH2O = 0
    cptSO4 = 0
    cptSO4fi = 0
    cptSO4fp = 0
    cptNafi = 0
    cptres = 1

    XYZ = np.zeros((1000000,3))
    Typ = ["" for x in range(1000000)]
    ResName = ["" for x in range(1000000)]
    ResNum = np.zeros((1000000,1))

    # Load Na2SO4 positions
    wallSO4Na2 = np.zeros((10000,7))
    file1 = open('Na2SO4/Position.dat', 'r')
    Lines = file1.readlines()
    count = 0
    for line in Lines:
        wallSO4Na2[count]=line.strip().split(' ')
        count += 1
    wallSO4Na2 = wallSO4Na2[0:count]


    wallSO4Na2rep = copy.deepcopy(wallSO4Na2)
    for xx in np.arange(txlo+lx/2,txhi,lx):
        for yy in np.arange(tylo+ly/2,tyhi,ly):
            for zz in np.arange(tzlo+lz/2,tzlo+nz*lz/2,lz):
                wallSO4Na2rep = np.append(wallSO4Na2rep,wallSO4Na2+[0,0,0,0,xx,yy,zz], axis=0)
    wallSO4Na2rep = wallSO4Na2rep[28:]
    assert len(wallSO4Na2rep[wallSO4Na2rep.T[2]==2]) == len(wallSO4Na2rep[wallSO4Na2rep.T[2]==3])/4

    for n in range(len(wallSO4Na2rep)):
        wallSO4Na2rep[n,0] = np.int64(n+1)

    shift = np.min(wallSO4Na2rep.T[6])
    wallSO4Na2rep.T[6] -= shift
    layer = np.max(wallSO4Na2rep.T[6]) - np.min(wallSO4Na2rep.T[6])

    wallSO4 = wallSO4Na2rep[wallSO4Na2rep.T[2] > 1]

    # choose the SO4 to maintain fix
    selection = wallSO4[(wallSO4.T[6]>layer*0.35) & (wallSO4.T[6]<layer*0.65) & (wallSO4.T[2] == 2)][1]
    assert len(selection)>0

    for n in np.arange(selection[0]-4, selection[0]+1,1):
        id0 = np.where(wallSO4Na2rep.T[0] == n)
        wallSO4Na2rep[id0,2] += 2

    wallSO4 = []
    wallSO4 = wallSO4Na2rep[wallSO4Na2rep.T[2] > 1]
    TypSO4 = ['O1', 'O2', 'O3', 'O4', 'S']
    wallNa = wallSO4Na2rep[wallSO4Na2rep.T[2] == 1]       

    # place the SO4 
    cpSO = 0
    for m in wallSO4:
        if m[2]<4:
            x0 = m[4]
            y0 = m[5]
            z0 = m[6]
            XYZ[cptatom] = [x0,y0,z0]
            nO = (cpSO)%5+1
            cpSO += 1
            Typ[cptatom] = TypSO4[nO-1]
            ResNum[cptatom] = cptres
            ResName[cptatom] = 'SO4fi' # fixed during the initial stage
            if Typ[cptatom] == 'S':
                cptSO4fi += 1
                cptres += 1
            cptatom += 1

    PosSO4 = np.array([[1.238,   0.587,   1.119],            [0.778,   1.501,  -1.263],            [-0.962,   1.866,   0.623],            [-0.592,  -0.506,  -0.358],           [0.115,   0.862,   0.030]])

    # add SO4 randomly
    nSO4 = 0
    while nSO4 < nion:
        x = random.randint(1,1000)/1000*(txhi-txlo)+txlo
        y = random.randint(1,1000)/1000*(tyhi-tylo)+tylo
        z = random.randint(1,1000)/1000*(tzhi-tzlo)+tzlo

        XYZrep = copy.deepcopy(XYZ[0:cptatom])
        for xx in [-Lx,0,Lx]:
            for yy in [-Ly,0,Ly]:
                for zz in [-Lz,0,Lz]:
                    XYZrep = np.append(XYZrep,XYZ[0:cptatom]+[xx,yy,zz], axis=0)
        d = np.sqrt((x-XYZrep.T[0])**2+(y-XYZrep.T[1])**2+(z-XYZrep.T[2])**2)

        if np.min(d) > 5 and z>tzlo+layer+3 and z<tzhi-3:
            for j in range(5):
                XYZ[cptatom] = [x,y,z]+PosSO4[j]
                Typ[cptatom] = TypSO4[j]
                ResNum[cptatom] = cptres
                ResName[cptatom] = 'SO4' # never fixed
                cptatom += 1  
            nSO4 += 1
            cptSO4 += 1
            cptres += 1

    # place the SO4 fixed
    cpSO = 0
    for m in wallSO4:
        if m[2]>=4:
            x0 = m[4]
            y0 = m[5]
            z0 = m[6]
            XYZ[cptatom] = [x0,y0,z0]
            nO = (cpSO)%5+1
            cpSO += 1
            Typ[cptatom] = TypSO4[nO-1]
            ResNum[cptatom] = cptres
            ResName[cptatom] = 'SO4fp' # fixed during the production run   
            if Typ[cptatom] == 'S':
                cptSO4fp += 1
                cptres += 1
            cptatom += 1  

    print(str(cptSO4fp+cptSO4fi+cptSO4)+' SO4 ions')

    # place the Na
    for m in wallNa:
        x0 = m[4]
        y0 = m[5]
        z0 = m[6]
        XYZ[cptatom] = [x0,y0,z0]
        Typ[cptatom] = 'Na'
        ResNum[cptatom] = cptres
        ResName[cptatom] = 'Nafi' # fixed during the initial stage
        cptatom += 1
        cptNafi += 1
        cptres += 1

    # add Na randomly
    nNa = 0
    while nNa < nion*2:
        x = random.randint(1,1000)/1000*(txhi-txlo)+txlo
        y = random.randint(1,1000)/1000*(tyhi-tylo)+tylo
        z = random.randint(1,1000)/1000*(tzhi-tzlo)+tzlo

        XYZrep = copy.deepcopy(XYZ[0:cptatom])
        for xx in [-Lx,0,Lx]:
            for yy in [-Ly,0,Ly]:
                for zz in [-Lz,0,Lz]:
                    XYZrep = np.append(XYZrep,XYZ[0:cptatom]+[xx,yy,zz], axis=0)
        d = np.sqrt((x-XYZrep.T[0])**2+(y-XYZrep.T[1])**2+(z-XYZrep.T[2])**2)
        if np.min(d) > 3 and z>tzlo+layer+1.5 and z<tzhi-1.5:
            XYZ[cptatom] = [x,y,z]
            Typ[cptatom] = 'Na'
            ResNum[cptatom] = cptres
            ResName[cptatom] = 'Na'
            cptatom += 1   
            nNa += 1
            cptres += 1
            cptNa += 1

    print(str(cptNa+cptNafi)+' Na ions') 

    cptH2O = 0
    # create water
    PosH2O = np.array([[0, 0, 0],            [0.05858,   0.0757, 0.0],            [0.05858,   -0.0757,  0.0],            [0.0104,  0.0, 0.0]])*10

    XYZrep = copy.deepcopy(XYZ[0:cptatom])
    for xx in [-Lx,0,Lx]:
        for yy in [-Ly,0,Ly]:
            for zz in [-Lz,0,Lz]:
                XYZrep = np.append(XYZrep,XYZ[0:cptatom]+[xx,yy,zz], axis=0)

    TypH2O = ['OW', 'HW1', 'HW2', 'MW']
    for x in np.arange(txlo+dw/2,txhi,dw):
        for y in np.arange(tylo+dw/2,tyhi,dw):
            for z in np.arange(tzlo+dw/2,tzhi,dw):

                d = np.sqrt((x-XYZrep.T[0])**2+(y-XYZrep.T[1])**2+(z-XYZrep.T[2])**2)

                if np.min(d) > dw and z>tzlo+layer+3 and z<tzhi-3 and cptH2O < nwater:
                    for j in range(4):
                        XYZ[cptatom] = [x,y,z]+np.array(PosH2O[j])
                        Typ[cptatom] = TypH2O[j]
                        ResNum[cptatom] = cptres
                        ResName[cptatom] = 'SOL'
                        cptatom += 1    
                    cptH2O += 1
                    cptres += 1
    print(str(cptH2O)+' water molecule')

    Vwater = cptH2O/6.022e23*0.018 # kg or litter
    Naddion = (nSO4+nNa)/6.022e23 # mol
    if Vwater>0:
        cion = Naddion/Vwater
        print('The initial ion concentration is '+str(cion)+' M')
    attemp += 1


# ## Displace hafl the NaCl ions of the side

# In[10]:


xyz_s = []
xyz_na = []
for cpt, mytyp in enumerate(Typ):
    if mytyp == 'S':
        xyz_s.append(XYZ[cpt])
    if mytyp == 'Na':
        xyz_na.append(XYZ[cpt])
xyz_s = np.array(xyz_s)
xyz_na = np.array(xyz_na)

S_histo = np.histogram(xyz_s.T[2], bins = 50)
Na_histo = np.histogram(xyz_na.T[2],bins = 50)
plt.plot((S_histo[1][:-1] + S_histo[1][1:])/2, S_histo[0], 'o')
plt.plot((Na_histo[1][:-1] + Na_histo[1][1:])/2, Na_histo[0], '.')


# In[11]:


to_displace = []
for cpt, mytyp in enumerate(Typ):
    if (mytyp == 'Na') & (XYZ[cpt][2] == 0) :
        to_displace.append(1)
    else:
        to_displace.append(0)


# In[12]:


displaced_atoms = 0
cpt_main = 0
for cpt, to_d in enumerate(to_displace):
    if (to_d == 1):
        cpt_main += 1
        if cpt_main%2 == 0:
            if double == True:
                XYZ[cpt][2] += 14.29813*2
            else:
                XYZ[cpt][2] += 14.29813
            displaced_atoms += 1
print(str(displaced_atoms) + ' atom displaced')


# In[13]:


xyz_s = []
xyz_na = []
for cpt, mytyp in enumerate(Typ):
    if mytyp == 'S':
        xyz_s.append(XYZ[cpt])
    if mytyp == 'Na':
        xyz_na.append(XYZ[cpt])
xyz_s = np.array(xyz_s)
xyz_na = np.array(xyz_na)

S_histo = np.histogram(xyz_s.T[2], bins = 50)
Na_histo = np.histogram(xyz_na.T[2],bins = 50)
plt.plot((S_histo[1][:-1] + S_histo[1][1:])/2, S_histo[0], 'o')
plt.plot((Na_histo[1][:-1] + Na_histo[1][1:])/2, Na_histo[0], '.')


# ### log file

# In[14]:


f = open("log.system", "w")
f.write('Lx = '+str(Lx/10)+' nm, Ly = '+str(Ly/10)+' nm, Lz = '+str(Lz/10)+' nm, hinit = '+str(h/10)+' nm')
f.write('They are '+ str(nSO4)+' additional SO4')
f.write('They are '+ str(nNa)+' additional Na')
f.write('They are '+ str(cptH2O)+' water molecules')
f.write('The initial ion concentration is '+str(cion)+' M')
f.close()


# ### write conf.gro

# In[15]:


f = open('conf.gro', 'w')
f.write('SO4Na2 slit pore\n')
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


# ### write topol.top

# In[16]:


f = open('topol.top', 'w')
f.write('#include "../ff/forcefield.itp"\n')
f.write('#include "../ff/tip4peps.itp"\n')
f.write('#include "../ff/ions.itp"\n\n')
f.write('[ System ]\n')
f.write('SO4Na2 slit pore\n\n')
f.write('[ Molecules ]\n\n')
f.write('SO4fi '+ str(cptSO4fi)+'\n')
f.write('SO4 '+ str(cptSO4)+'\n')
f.write('SO4fp '+ str(cptSO4fp)+'\n')
f.write('Nafi '+ str(cptNafi)+'\n')
f.write('Na '+ str(cptNa)+'\n')
f.write('SOL '+ str(cptH2O)+'\n')
f.close()


# In[ ]:





# In[ ]:




