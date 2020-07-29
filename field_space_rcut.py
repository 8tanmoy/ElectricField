#--field calculation for convergence test--
#--tanmoy jun 19, 2020--
#--intra-CN removed upon upload jul 29, 2020--

"""
tasks
_____
for each rcut:
    1. decomposition of the electric field. output just the field due to cation + anion
"""
import shutil
import matplotlib.pyplot as plt
import os
import sys
import subprocess as sp
import numpy as np
import copy
import time
np.set_printoptions(threshold=sys.maxsize)

#--constants--
gro2si = 1.602176 * 9.0 * 0.1
print(f"gro2si = {gro2si}")

root    =   '/projectnb/cui-buchem/tanmoy/projects/ionic_liquids/setup_may2/1_emim_ace/'
drf     =   root + 'field_rcut'
wrk     =   os.getcwd()
r_cut   =  0.1 + 0.05 * _RCUTNO_


#--process gro file--
#--inp : .gro file path--
#--out : 1> timestamp 2> boxinfo 3> n_atoms (debug) 4> coordonate array
def process_gro(coord_path):
    f1 = open(coord_path, "r")
    arr = list()
    for line in f1:
        arr.append(line)
    f1.close()
    tstamp = arr[0]
    nat = int(arr[1])
    boxinfo = arr[-1].split()
    arr = arr[2:-1]
    molinfo = list()        #information array
    #--06/22--
    for jj in arr:
        temp_resid      = jj[0:5].strip()
        temp_resname    = jj[5:10].strip()
        temp_atname     = jj[10:15].strip()
        temp_idx        = jj[15:20].strip()
        molinfo.append([temp_resid, temp_resname, temp_atname, temp_idx])
    #print(molinfo)
    arr = [ii.split() for ii in arr]
    arr = [ii[-3:] for ii in arr]                                           #new
    assert nat == len(arr), "natoms error in process_gro function"
    assert len(arr[0]) == len(arr[-1]), "len error in process_gro"
    assert len(molinfo) == len(arr), "len molinfo != len coordinates"
    #--zip molinfo and coordinates--
    for ii in range(len(molinfo)):
        molinfo[ii].extend(arr[ii])
    return(tstamp, nat, boxinfo, molinfo)

def read_charges(ch_path):
    f1 = open(ch_path, "r")
    arr = list()
    for line in f1:
        arr.append(line.split())
    # --resname--atomanme--charge--
    f1.close()
    return(arr)

def calc_maxmin(arr_coord):
    x_arr   = [ float(row[0]) for row in arr_coord]
    y_arr   = [ float(row[1]) for row in arr_coord]
    z_arr   = [ float(row[2]) for row in arr_coord]
    print(f"xmax = {np.max(x_arr)}\t\
xmin = {np.min(x_arr)}\n\
ymax = {np.max(y_arr)}\t\
ymin = {np.min(y_arr)}\n\
zmax = {np.max(z_arr)}\t\
zmin = {np.min(z_arr)}\n")
    return

def get_ref_idx(fname):
    f1 = open(fname, "r")
    data_idx    = f1.read().replace('\n',' ')
    data_idx    = data_idx.split()
    f1.close()
    print("ref array passed from ", fname)
    return(data_idx)

def strip_gold(withGold):
    noGold = list()
    for item in withGold:
        #--item[2] is atomname--
        if item[2][:2] != 'AU':
            noGold.append(item)
    return(noGold)

def add_charges(noCharge, chargeList):
    for item in noCharge:
        selectAtom = [ float(ii[2]) for ii in chargeList if (ii[0] == item[1] and ii[1] == item[2]) ]
        item.append(selectAtom[0])  #because using list comprehension
    return(noCharge)

def calc_field(refIndex, noGold, boxDim, chargeList, rcut):
    inMol       = copy.deepcopy(noGold)
    ref_atom    = [ item for item in inMol if item[3] == refIndex ]
    ref_atom    = ref_atom[0]
    print(f'using reference: {ref_atom}')
    inMol.remove(ref_atom)                 #inplace removal
    ref_atom[4], ref_atom[5], ref_atom[6] = float(ref_atom[4]), float(ref_atom[5]), float(ref_atom[6])
    boxDim = [float(ii) for ii in boxDim]
    #--calculate vector and use PBC--
    for atom in inMol:
        for cart in [4, 5, 6]:
            atom[cart] = float(atom[cart])
            atom[cart] -= ref_atom[cart]
        #--set pbc--
            if abs(atom[cart]) > boxDim[cart - 4]*0.5:
                if atom[cart] > boxDim[cart - 4]*0.5:
                    atom[cart] -= boxDim[cart - 4]
                else:
                    atom[cart] += boxDim[cart - 4]

    #--get a new list with max distance rcut--
    noGoldCut = list()
    for atom in inMol:
        dr2     = atom[4]*atom[4] + atom[5]*atom[5] + atom[6]*atom[6]
        if dr2 < rcut*rcut:
            noGoldCut.append(atom)
    noGoldCut = add_charges(noGoldCut, chargeList)
    #--loop over atoms and collect field with conditionals--
    fIntraX,  fIntraY,  fIntraZ     = 0.0, 0.0, 0.0
    fLigandX, fLigandY, fLigandZ    = 0.0, 0.0, 0.0
    fCationX, fCationY, fCationZ    = 0.0, 0.0, 0.0
    fAnionX,  fAnionY,  fAnionZ     = 0.0, 0.0, 0.0
    fTotalX,  fTotalY,  fTotalZ     = 0.0, 0.0, 0.0
    countIntra, countLigand, countCation, countAnion, countTotal = 0, 0, 0, 0, 0
    print(f'len inMol : {len(inMol)} \nlen noGoldCut : {len(noGoldCut)} \n')
    print(f'len inMol[0] : {len(inMol[0])}')
    for atom in noGoldCut:
        #print(atom)
        countTotal += 1
        dr2     = atom[4]*atom[4] + atom[5]*atom[5] + atom[6]*atom[6]
        dr3i    = dr2**(-1.5)
        #--check for intra. Match 0resid and 1resname
        if (atom[0] == ref_atom[0] and atom[1] == ref_atom[1]):
            countIntra += 1
            fIntraX    += atom[7] * dr3i * atom[4]
            fIntraY    += atom[7] * dr3i * atom[5]
            fIntraZ    += atom[7] * dr3i * atom[6]
        #--check for ligand
        elif (atom[1] == 'MBN'):
            countLigand += 1
            fLigandX    += atom[7] * dr3i * atom[4]
            fLigandY    += atom[7] * dr3i * atom[5]
            fLigandZ    += atom[7] * dr3i * atom[6]            
        #--check for cation, assuming we are just using EMIM
        elif (atom[1] == 'EMIM'):
            countCation += 1
            fCationX    += atom[7] * dr3i * atom[4]
            fCationY    += atom[7] * dr3i * atom[5]
            fCationZ    += atom[7] * dr3i * atom[6]            
        #--check for anion
        else:
            countAnion += 1
            fAnionX    += atom[7] * dr3i * atom[4]
            fAnionY    += atom[7] * dr3i * atom[5]
            fAnionZ    += atom[7] * dr3i * atom[6]
        #--total field--
        fTotalX    += atom[7] * dr3i * atom[4]
        fTotalY    += atom[7] * dr3i * atom[5]
        fTotalZ    += atom[7] * dr3i * atom[6]
    print(f'countTotal : {countTotal} \ncountIntra : {countIntra} \ncountLigand : {countLigand} \ncountCation : {countCation} \ncountAnion : {countAnion} \n')
    assert countTotal == (countIntra + countLigand + countCation + countAnion), 'error in decomposition'
    #--collect arrays--
    fIntra      = gro2si * np.array([fIntraX, fIntraY, fIntraZ])
    fLigand     = gro2si * np.array([fLigandX, fLigandY, fLigandZ])
    fCation     = gro2si * np.array([fCationX, fCationY, fCationZ])
    fAnion      = gro2si * np.array([fAnionX, fAnionY, fAnionZ])
    fTotal      = gro2si * np.array([fTotalX, fTotalY, fTotalZ])
    fTotNoIntra = fCation + fAnion
    return(fIntra, fLigand, fCation, fAnion, fTotNoIntra)

#--get coordinates as gro file--
def gen_trajout(time):
    gmx_trj_cmd = "echo 0 | gmx trjconv -f " + root + "prod.xtc -s " + root + "prod.tpr"\
                  + " -b " + str(time) + " -e " + str(time) + " -pbc none"\
                  + " -o trajout.gro"
    #center or boxcenter does not make a difference
    print(f'gmx_trj_cmd = {gmx_trj_cmd}')
    sp.check_call(gmx_trj_cmd, shell=True)
    return

simTime        = 100000
gen_trajout(simTime)

#--process gro file--
tstamp, n_atoms, boxInfo, molinfo = process_gro(wrk + '/trajout.gro')
idx_ref      = get_ref_idx( drf + '/index_nz.ndx')
chargeList  = read_charges( drf + '/charges.dat')
noGold      = strip_gold(molinfo)
#print(chargeList)
#print(idx_ref)

#test = add_charges(noGold, chargeList)
#print(test[0])

tbeg    = time.perf_counter()
#--make files for field output--
#f01     = open("field_intra.dat", "w")
#f02     = open("field_ligand.dat", "w")
#f03     = open("field_cation.dat", "w")
#f04     = open("field_anion.dat", "w")
f05     = open("field_nointra.dat", "w")
#--read charges from charges
for index in idx_ref:
    fIntra, fLigand, fCation, fAnion, fTotNoIntra = calc_field(index, noGold, boxInfo, chargeList, r_cut)
    #print(fIntra, fLigand, fCation, fAnion, fTotNoIntra, sep='\n')
    #f01.write(" ".join(str(x) for x in fIntra))
    #f01.write("\n")
    #f02.write(" ".join(str(x) for x in fLigand))
    #f02.write("\n")
    #f03.write(" ".join(str(x) for x in fCation))
    #f03.write("\n")
    #f04.write(" ".join(str(x) for x in fAnion))
    #f04.write("\n")
    f05.write(" ".join(str(x) for x in fTotNoIntra))
    f05.write("\n")
#f01.close()
#f02.close()
#f03.close()
#f04.close()
f05.close()
tend    = time.perf_counter()
print(f'took {tend - tbeg} s for {len(idx_ref)} refs')
