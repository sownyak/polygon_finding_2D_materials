import sys
import re
import numpy as np
import math as m

#Please cite the paper 
# @article{mondal2024quantifying,
#   title={Quantifying defects in graphene oxide structures},
#   author={Mondal, Sownyak and Ghosh, Soumya},
#   journal={Carbon Trends},
#   volume={14},
#   pages={100323},
#   year={2024},
#   publisher={Elsevier}
# }

# purpose: find indices of the vertices of polygons
# usage: polygons.py image_with_connectivity.txt position.xyz output_polygons.txt

file1 = sys.argv[1] # inputfile: connectivity file: 1. with indices arranged in ascending order 2. index of atoms with 1 neighbor removed 3. coordinates of the neighbors, adjusting for periodic images, included. format: [index  neighbor_index neighbor_coordinates ....]
file2 = sys.argv[2] # actual xyz file
file3 = sys.argv[3] # outputfile

numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
rx = re.compile(numeric_const_pattern, re.VERBOSE)

eps=float(1E-5)
# collect both inputfiles in memory (1 snapshot at a time)
fo1 = open(file1)
lfo1 = fo1.readlines()
fo2 = open(file2)
lfo2 = fo2.readlines()
# go through the connectivity file ---> for each index (i) generate sets (at max 6) of 3 indices (jkl) corresponding to neighbours
nindices = len(lfo1)
print('nindices: ', nindices)

# function to compute dihedral
def calc_dih(pos_4points, eps):
    p1 = [x for x in pos_4points[0]]
    p2 = [x for x in pos_4points[1]]
    p3 = [x for x in pos_4points[2]]
    p4 = [x for x in pos_4points[3]]
    vec1 = np.array([x-y for x,y in zip(p1,p2)])
    vec2 = np.array([x-y for x,y in zip(p4,p3)])
    cos12 = float(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
    if (abs(cos12) > eps):
        coscorr12 = [cos12, pos_4points]
    else:
        coscorr12 = [float(0), pos_4points]
    return coscorr12

# function to compute sets of 2 from neighborlist
def create_set(list_neighbors):
    tmp_lneighbors = [x for x in list_neighbors]
    tmp_lneighbors.sort()
    sets = []
    i = 0
    while i < (len(tmp_lneighbors)-1):
        j = i + 1
        while j < len(tmp_lneighbors):
            sets.append([tmp_lneighbors[i],tmp_lneighbors[j]])
            j = j + 1
        i = i + 1
    return sets
        

pos_v_dih = []
list_cycles = []
ind_cycle = []
incomplete_sets = []
#incomplete_sets_trial = []
tmp_incomplete_sets = []
tmp_complete_sets = []
all_indices = []
ind_neighbors = []
pos_neighbors = []
dummy_ind_neighbors = []
dummy_pos_neighbors = []
neighbor_pos = []
list_pos_cycles = []
pos_cycle = []
for i in range(nindices):
    floats = [float(x) for x in rx.findall(lfo1[i])]
    curr_index = int(floats[0])
    all_indices.append(curr_index)
for i in range(nindices):    # nindices can be less than number of atoms if we exclude atoms wth 1 neighbor
    floats = [float(x) for x in rx.findall(lfo1[i])]
    icurr = int(floats[0])      
    #print('current index: ', icurr)
    fneighbors = [x for x in floats[1:]]
    #print('fneighbors: ', fneighbors)
    #print('No. of neighbors: ', len(sneighbors))
    nneighbors = int(len(fneighbors)/4)
    ind_neighbors[:] = []
    pos_neighbors[:] = []
    for j in range(nneighbors):
        neighbor = int(fneighbors[j*4])
        ind_neighbors.append(neighbor)
        neighbor_pos[:] = []
        for k in range(1,4,1):
            jk = k + j*4
            #print(jk)
            coord = fneighbors[jk]
            neighbor_pos.append(coord)
        #print('neighbor_pos: ', neighbor_pos)
        tmp_pos = [x for x in neighbor_pos]
        pos_neighbors.append(tmp_pos)
    #print('ind_neighbors: ', ind_neighbors)
    #print('pos_neighbors: ', pos_neighbors)
    ind_sets = create_set(ind_neighbors)
    #print('sets:')
    #print(ind_sets)
    # start with each set and choose a neighbor
    #tmp_incomplete_sets[:] = []
    # cycle through the different sets
    #if (icurr == 1):
    #    print('icurr, ind_sets: ', icurr, ',', ind_sets)
    for item0 in ind_sets:
        #if (icurr == 1):
        #    print('current set = ', item0)
        icycle = 0
        # if the elements of the sets could not be connected by going in one direction, then try the other direction
        while icycle < 2:
            pos_v_dih[:] = []
            ind_cycle[:] = []
            pos_cycle[:] = []
            if (icycle == 0):
                v1 = item0[0]
                v2 = item0[1]
            else:
                v1 = item0[1]
                v2 = item0[0]
            ind_cycle.append(v2)
            ind_cycle.append(icurr)
            ind_cycle.append(v1)
            ref_ind = v2
            #if (icurr == 1):
            #    print('icycle = ', icycle)
            #    print('ref_ind = ', ref_ind)
            dummy_ind = [v2, icurr, v1]
            #if (icurr == 1):
                #print('icycle, dummy_ind: ', icycle, ',', dummy_ind)
            #dummy_pos[:] = []
            for j in dummy_ind:
                #print(j)
                if (j == icurr):
                    coords = [float(x) for x in rx.findall(lfo2[j+1])]
                    pos_v_dih.append(coords)
                else:
                    #print('ind_neighbors: ', ind_neighbors)
                    icoord = ind_neighbors.index(j)
                    #print(icoord)
                    coord = pos_neighbors[icoord]
                    #print(coord)
                    pos_v_dih.append(coord)
            pos_cycle = [x for x in pos_v_dih]
            #print('pos_v_dih: ', pos_v_dih)
            size_cycle = len(ind_cycle)
            # look to generate a polygon
            while size_cycle < nindices:
                #if (icurr == 1):
                #    print('ind_cycle: ', ind_cycle)
                #print('size_cycle: ', size_cycle)
                #print('dummy_ind[-1]: ', dummy_ind[-1])
                #print('list_cycles: ', list_cycles)
                index_v1 = all_indices.index(dummy_ind[-1])
                dummy_sneighbors = [x for x in rx.findall(lfo1[index_v1])[1:]]
                ndneighbors = int(len(dummy_sneighbors)/4)
                dummy_ind_neighbors[:] = []
                dummy_pos_neighbors[:] = []
                for j in range(ndneighbors):
                    dummy_ind_neighbors.append(int(dummy_sneighbors[j*4]))
                    neighbor_pos[:] = []
                    for k in range(1,4,1):
                        jk = k + j*4
                        neighbor_pos.append(float(dummy_sneighbors[jk]))
                    tmp_neighbor_pos = [x for x in neighbor_pos]
                    dummy_pos_neighbors.append(tmp_neighbor_pos)
                #********need to be corrected for translation ************
                trans_dist = float(0.0)
                orig_coord = [float(x) for x in rx.findall(lfo2[dummy_ind[-1]+1])]
                curr_coord = [x for x in pos_v_dih[-1]]
                test_r = [x-y for x,y in zip(orig_coord,curr_coord)]
                test_norm = np.linalg.norm(test_r)
                #if (icurr == 1):
                #    print('dummy_ind[-1]: ', dummy_ind[-1])
                #    print('orig_coord: ', orig_coord)
                #    print('curr_coord: ', curr_coord)
                #    print('test_r: ', test_r)
                #    print('test_norm: ', test_norm)
                if (test_norm > eps):
                    #print('coords need translation')
                    #print('icurr, current set, dummy vertex: ', icurr, ',', item0, ',', dummy_ind[-1])
                    tmp_dummy_pos_neighbors = [x for x in dummy_pos_neighbors]
                    dummy_pos_neighbors[:] = []
                    for pos in tmp_dummy_pos_neighbors:
                        tmp_pos = [x-y for x,y in zip(pos,test_r)]
                        dummy_pos_neighbors.append(tmp_pos)
                #if (icurr == 1):
                   # print('dummy_ind[-1], dummy_ind_neighbors: ', dummy_ind[-1], ', ', dummy_ind_neighbors)
                i1dummy = dummy_ind_neighbors.index(dummy_ind[1])
                dummy_ind_neighbors.remove(dummy_ind[1])
                dummy_pos_neighbors.pop(i1dummy)
                #print('dummy_ind_neighbors after removing the immediate neighbor: ', dummy_ind_neighbors)
                # check whether there is a smaller cycle
                tmp_ind_cycle = [i for i in ind_cycle]
                tmp_pos_cycle = [x for x in pos_cycle]
                ind_ref = tmp_ind_cycle.index(ref_ind)
                #tmp_ind_cycle.remove(ref_ind)
                tmp_ind_cycle.pop(ind_ref)
                tmp_pos_cycle.pop(ind_ref)
                ind_member = -1
                for member in dummy_ind_neighbors:
                    if member in tmp_ind_cycle:
                        ind_member = tmp_ind_cycle.index(member)
                        break
                #print('ind_member = ', ind_member)
                if (ind_member > -1):
                    tmp3 = [i for i in tmp_ind_cycle[ind_member:]]
                    ind_cycle = [i for i in tmp3]
                    pos_tmp3 = [x for x in tmp_pos_cycle[ind_member:]]
                    pos_cycle = [x for x in pos_tmp3]
                    #if (icurr == 1):
                    #    print('ind_member = ', ind_member)
                    tmp3.sort()
                    #print('LIST_CYCLES: ', list_cycles)
                    if len(list_cycles) > 0:
                    # check for repitition
                        icollect = 1
                        for item2 in list_cycles:
                            tmp4 = [i for i in item2]
                            tmp4.sort()
                            #print('old_cycle: ', item2)
                            if (tmp3 == tmp4):
                                icollect = 0
                                #if (icurr == 1):
                                #    print('smaller cycle repeated')
                                break
                        if (icollect == 1):
                            #if (icurr == 1):
                            #    print('NEW CYCLE (smaller): ', ind_cycle)
                            #if (len(ind_cycle) > 6):
                            #    print('POLYGON LARGER THAN HEXAGON: ', ind_cycle)
                            #if (len(ind_cycle) == 5):
                            #    print('PENTAGON: ', ind_cycle)
                            tmp5 = [i for i in ind_cycle]
                            list_cycles.append(tmp5)
                            pos_tmp5 = [x for x in pos_cycle]
                            list_pos_cycles.append(pos_tmp5)
                    else:
                        #if (icurr == 1):
                        #    print('NEW CYCLE (smaller): ', ind_cycle)
                        #print('NEW CYCLE: ', ind_cycle)
                        #if (len(ind_cycle) > 6):
                        #    print('POLYGON LARGER THAN HEXAGON: ', ind_cycle)
                        #if (len(ind_cycle) == 5):
                        #    print('PENTAGON: ', ind_cycle)
                        tmp = [i for i in ind_cycle]
                        list_cycles.append(tmp)
                        pos_tmp = [x for x in pos_cycle]
                        list_pos_cycles.append(pos_tmp)
                    #icycle = icycle + 1
                    break
                else:
                    if (ref_ind in dummy_ind_neighbors):
                        #print('LIST_CYCLES: ', list_cycles)
                        # check for repitition
                        if len(list_cycles) > 0:
                            # check for repitition
                            icollect = 1
                            tmp1 = [i for i in ind_cycle]
                            tmp1.sort()
                            for item1 in list_cycles:
                                tmp2 = [i for i in item1]
                                tmp2.sort()
                                #print('tmp1: ', tmp1)
                                #print('tmp2: ', tmp2)
                                #print('old_cycle: ', item1)
                                if (tmp1 == tmp2):
                                    icollect = 0
                                    #if (icurr == 1):
                                    #    print('cycle repeated')
                                    break
                            if (icollect == 1):
                                #if (icurr == 1):
                                #    print('NEW CYCLE: ', ind_cycle)
                                #if (len(ind_cycle) > 6):
                                #    print('POLYGON LARGER THAN HEXAGON: ', ind_cycle)
                                #if (len(ind_cycle) == 5):
                                #    print('PENTAGON: ', ind_cycle)
                                tmp = [i for i in ind_cycle]
                                list_cycles.append(tmp)
                                pos_tmp = [x for x in pos_cycle]
                                list_pos_cycles.append(pos_tmp)
                        else:
                            #if (icurr == 1):
                            #    print('NEW CYCLE: ', ind_cycle)
                            #if (len(ind_cycle) > 6):
                            #    print('POLYGON LARGER THAN HEXAGON: ', ind_cycle)
                           #if (len(ind_cycle) == 5):
                           #    print('PENTAGON: ', ind_cycle)
                            tmp = [i for i in ind_cycle]
                            list_cycles.append(tmp)
                            pos_tmp = [x for x in pos_cycle]
                            list_pos_cycles.append(pos_tmp)
                        icycle = 3
                        tmp_set = [v1,v2]
                        tmp_set.sort()
                        icollect = 1
                        if len(tmp_complete_sets) > 0:
                            for item in tmp_complete_sets:
                                if (tmp_set == item):
                                    icollect = 0
                                    break
                            if (icollect == 1):
                                tmp_complete_sets.append(tmp_set)
                        else:
                            tmp_complete_sets.append(tmp_set)
                        break
                    else:
                        if (len(dummy_ind_neighbors) == 0):
                            break
                        elif (len(dummy_ind_neighbors) == 1):
                            kmin = dummy_ind_neighbors[0]
                        else:
                            max_cosdih = float(-2)
                            kmin = 0
                            # find neighbors of the neighbor and compute the min dihedral
                            indk = 0
                            while indk < len(dummy_ind_neighbors):
                                dummy_indk = dummy_ind_neighbors[indk]
                                pos_k = dummy_pos_neighbors[indk]
                                #print('len(pos_k) = ', len(pos_v_dih))
                                #print('pos_v_dih, pos_k: ', pos_v_dih, pos_k)
                                tmp6 = pos_v_dih + [pos_k]
                                #print('tmp6: ', tmp6)
                                if len(tmp6) > 4:
                                    print('icurr, k, dummy_ind_neighbors: ', icurr, k, dummy_ind_neighbors)
                                    print(tmp6)
                                    sys.exit('problem in the code')
                                tmp_cos = [x for x in calc_dih(tmp6, eps)]
                                cosdih = float(tmp_cos[0])
                                #if (len(ind_cycle) >=5):
                                #    print('dummy_indk = ', dummy_indk)
                                #    print('pos_4points: ', tmp6)
                                #    print('cosdih = ', cosdih)
                                #if (abs(cosdih) < eps):
                                #    sys.exit('Dihedral angle is less than cutoff')
                                #elif (cosdih < eps):
                                #    print('k = ', k)
                                dcosdih = cosdih - max_cosdih
                                if (dcosdih > eps):
                                    max_cosdih = cosdih
                                    kmin = dummy_indk
                                #print('kmin = ', kmin)
                                indk = indk + 1

                        #if (len(ind_cycle) >= 5):
                            #print('kmin = ', kmin)
                        ind_kmin = dummy_ind_neighbors.index(kmin)
                        pos_kmin = dummy_pos_neighbors[ind_kmin]
                        pos_v_dih.pop(0)
                        pos_v_dih.append(pos_kmin)
                        ind_cycle.append(kmin)
                        pos_cycle.append(pos_kmin)
                        dummy_ind.pop(0)
                        dummy_ind.append(kmin)
                        size_cycle = len(ind_cycle)

            icycle = icycle + 1
            #if (icurr == 1):
            #    print('icycle after exiting the cycle loop', icycle)
            # if the endpoints do not match even after traversing in either direction, save the set
            if (icycle == 2):
                #print('incomplete set: ', [v1,v2])
                tmp_set = [v1,v2]
                tmp_set.sort()
                icollect = 1
                if len(tmp_incomplete_sets) > 0:
                    for item in tmp_incomplete_sets:
                        if (tmp_set == item[1]):
                            icollect = 0
                            break
                    if (icollect == 1):
                        tmp_incomplete_sets.append([icurr,tmp_set])
                        #if (icurr == 1):
                        #    print('tmp_incomplete_sets: ', tmp_incomplete_sets)
                else:
                    tmp_incomplete_sets.append([icurr,tmp_set])

numcycles = len(list_cycles)
print('Total No. of cycles obtained by considering minimum dihedral = ', numcycles)

lenincomptmp = len(tmp_incomplete_sets)
print('Total No. of trial incomplete sets = ', lenincomptmp)
lencomp = len(tmp_complete_sets)
print('Total No. of complete sets = ', lencomp)

# check whether the incomplete sets are real? ***** Start here; not a trivial task **********
incomplete_sets = [x for x in tmp_incomplete_sets]
lenincomplete = len(incomplete_sets)
print('Total No. of unparsed incomplete sets = ', lenincomplete)
print(incomplete_sets)
# check for incomplete sets
list_cycles2 = []
list_pos_cycles2 = []
dummy_ind[:] = []
ind_cycle[:] = []
pos_cycle[:] = []
incomplete_indices = []
lindneighbors = []
dummy_vert_set = []
parse_set = []
lpostneighbors = []
lindtneighbors = []
lposneighbors = []
pos_vert_dih = []
complete_pairs = []
#incomplete_sets[:] = []
if lenincomplete > 0:
    for sets in incomplete_sets:
        incomplete_indices.append(sets[0])
    #print('incomplete_indices: ', incomplete_indices)
    #print('No. of indices with incomplete cycles = ', len(incomplete_indices))
    # retry to form a polygon with the elements of the incomplete set but do not use dihedral as the deciding criterion
    iset = 0
    for sets in incomplete_sets:
        iset = iset + 1
        lenincomplete = len(incomplete_sets)
        irepeat = 0
        print('iset = ', iset)
        #print('incomplete_sets')
        #print(incomplete_sets)
        print('Size of the incomplete_sets')
        print(lenincomplete)
        lindtneighbors[:] = []
        lpostneighbors[:] = []
        pos_vert_dih[:] = []
        icurr_ind = sets[0]
        #if (iset == 8):
        #    print('icurr_incomplete = ', icurr_ind)
        ind_sets = [x for x in sets[1]]
        #if (iset == 8):
        #    print('ind_sets: ', ind_sets)
        # compute neighbors of both the indices and check whether a 3 membered cycle is formed
        lnelem = []
        lnnelem = []
        ireal = 1
        for element in ind_sets:
            lnelem[:] = []
            lnnelem[:] = []
            # find neighbors of the set elements
            ind_elem = all_indices.index(element)
            #print('ind_elem = ', ind_elem)
            floats = [x for x in rx.findall(lfo1[ind_elem])[1:]]
            #print('floats: ', floats)
            lenf = len(floats)
            #print(lenf)
            nnelem = int(lenf/4)
            #print(nnelem)
            # check whether the number of neighbors is 2
            if (nnelem == 2):
                # eliminate the anchor
                for k in range(nnelem):
                    lnelem.append(int(floats[k*4]))
                lnelem.remove(icurr_ind)
                # check whether the remaining neighbor forms a cycle with the anchor
                ind_nelem = all_indices.index(lnelem[0])
                floats2 = [x for x in rx.findall(lfo1[ind_nelem])[1:]]
                lenf2 = len(floats2)
                nnelem2 = int(lenf2/4)
                for k in range(nnelem2):
                    lnnelem.append(int(floats2[k*4]))
                if icurr_ind in lnnelem:
                    ireal = 0
                    break
        #if (iset == 8):
        #    print('ireal = ', ireal)
        if (ireal == 0):
            print('Set ', sets, 'is NOT real')
            ifset = incomplete_indices.index(icurr_ind)
            incomplete_sets.pop(ifset)
            incomplete_indices.pop(ifset)
        #    if (iset == 8):
        #        print(incomplete_sets)
            continue
        dummy_ind = [ind_sets[0],icurr_ind,ind_sets[1]]
        ref_ind = dummy_ind[0]
        #if (iset == 8):
        #    print('ref_ind = ', ref_ind)
        ind_cycle = [i for i in dummy_ind]
        ind_icurr = all_indices.index(icurr_ind)
        sfloats = [x for x in rx.findall(lfo1[ind_icurr])[1:]]
        lensfloats = len(sfloats)
        ntneighbors = int(lensfloats/4)
        for k in range(ntneighbors):
            lindtneighbors.append(int(sfloats[k*4]))
            lpostneighbors.append([float(x) for x in sfloats[k*4+1:k*4+4]])
        for j in dummy_ind:
            #print(j)
            if (j == icurr_ind):
                coords = [float(x) for x in rx.findall(lfo2[j+1])]
                pos_vert_dih.append(coords)
            else:
                indj = lindtneighbors.index(j)
                #print(icoord)
                jcoord = lpostneighbors[indj]
                #print(coord)
                pos_vert_dih.append(jcoord)
        pos_cycle = [x for x in pos_vert_dih]
        size_cycle = len(ind_cycle)
        while size_cycle < (nindices+1):
            #if (iset == 8):
            #    print('ind_cycle: ', ind_cycle)
            lindneighbors[:] = []
            lposneighbors[:] = []
            #lindtneighbors[:] = []
            #lpostneighbors[:] = []
            dummy_vert = dummy_ind[-1]
            test_ind = dummy_ind[1]
            if (iset == 41):
                print('dummy_vert = ', dummy_vert)
                print('test_ind = ', test_ind)
            ind_dvert = all_indices.index(dummy_vert)
            sdfloats = [x for x in rx.findall(lfo1[ind_dvert])[1:]]
            lensdfloats = len(sdfloats)
            ndneighbors = int(lensdfloats/4)
            for j in range(ndneighbors):
                lindneighbors.append(int(sdfloats[j*4]))
                lposneighbors.append([float(x) for x in sdfloats[j*4+1:j*4+4]])
            if ref_ind in lindneighbors:
                print('trial cycle: ', ind_cycle)
                break
            #print('lindneighbors: ', lindneighbors)
            # ************** consider translation here ***************
            trans_dist = float(0.0)
            orig_coord = [float(x) for x in rx.findall(lfo2[dummy_vert+1])]
            curr_coord = [x for x in pos_vert_dih[-1]]
            test_r = [x-y for x,y in zip(orig_coord,curr_coord)]
            test_norm = np.linalg.norm(test_r)
            #if (icurr == 1):
            #    print('dummy_ind[-1]: ', dummy_ind[-1])
            #    print('orig_coord: ', orig_coord)
            #    print('curr_coord: ', curr_coord)
            #    print('test_r: ', test_r)
            #    print('test_norm: ', test_norm)
            if (test_norm > eps):
                #print('coords need translation')
                #print('icurr, current set, dummy vertex: ', icurr, ',', item0, ',', dummy_ind[-1])
                tmp_lposneighbors = [x for x in lposneighbors]
                lposneighbors[:] = []
                for pos in tmp_lposneighbors:
                    tmp_pos = [x-y for x,y in zip(pos,test_r)]
                    lposneighbors.append(tmp_pos)
            lind_dummy_vert = [index for index, element in enumerate(incomplete_indices) if element == dummy_vert]
            icheck_incomplete = 0
            set_incomplete = []
            if len(lind_dummy_vert) > 0:
                set_incomplete[:] = []
                for i in lind_dummy_vert:
                    set_incomplete = [j for j in incomplete_sets[i][1]]
                    if test_ind in set_incomplete:
                        icheck_incomplete = 1
                        break
            if (icheck_incomplete == 1):
                print('lind_dummy_vert: ', lind_dummy_vert)
                if len(lind_dummy_vert) == 1:
                    ind_dummy_vert = lind_dummy_vert[0]
                    if (iset == 41):
                        print(len(incomplete_sets))
                        print(incomplete_sets)
                    dummy_set = [i for i in incomplete_sets[ind_dummy_vert][1]]
                    if (iset == 41):
                        print('dummy_set: ', dummy_set)
                    dummy_set.remove(test_ind)
                    new_ind = dummy_set[0]
                else:
                    dummy_vert_set[:] = []
                    parse_set[:] = []
                    #if (iset == 8):
                    #    print('incomplete_sets: ', incomplete_sets)
                    for iv in lind_dummy_vert:
                        print('relevant index and incomplete set: ', iv, ',', incomplete_sets[iv])
                        dummy_vert_set.append([x for x in incomplete_sets[iv][1]])
                        print('dummy_vert_set: ', dummy_vert_set)
                    for ipv in dummy_vert_set:
                        if test_ind in ipv:
                            parse_set.append(ipv)
                    print('test_ind, parse_set: ', test_ind, ',', parse_set)
                    if len(parse_set) == 1:
                        new_trial = [i for i in parse_set[0]]
                        new_trial.remove(test_ind)
                        new_ind = new_trial[0]
                    elif len(parse_set) > 1:
                        # choose from multiple options using dihedral angle
                        # ********* CHECK CRITERION **********
                        kmin = 0
                        maxdih = float(-2)
                        for isetp in parse_set:
                            isetp.remove(test_ind)
                            ik = isetp[0]
                            #print('isetp[0] = ', ik)
                            indik = lindneighbors.index(ik)
                            posik = [x for x in lposneighbors[indik]]
                            tmp_pos = pos_vert_dih + [posik]
                            #print('tmp_pos: ', tmp_pos)
                            if len(tmp_pos) != 4:
                                sys.exit('length of tmp_pos should be 4')
                            tmp_cos = [x for x in calc_dih(tmp_pos, eps)]
                            cosdih = float(tmp_cos[0])
                            #print('cosdih = ', cosdih)
                            dcosdih = cosdih - maxdih
                            if (dcosdih > eps):
                                maxdih = cosdih
                                kmin = ik
                        new_ind = kmin
                    else:
                        print('lind_dummy_vert: ', lind_dummy_vert)
                        print('dummy_vert_set: ', dummy_vert_set)
                        print('parse_set: ', parse_set)
                        sys.exit('check code')
            if (icheck_incomplete == 0):
                # consider general case
                #print('dummy_vert is NOT an incomplete index')
                ind_test = lindneighbors.index(test_ind)
                lindneighbors.pop(ind_test)
                lposneighbors.pop(ind_test)
                ndneighbors = len(lindneighbors)
                if (ndneighbors == 1):
                    new_ind = lindneighbors[0]
                elif (ndneighbors > 1):
                    kmin = 0
                    maxdih = float(-2)
                    for ik in range(ndneighbors):
                        posik = [x for x in lposneighbors[ik]]
                        tmp_pos = pos_vert_dih + [posik]
                        #print('tmp_pos: ', tmp_pos)
                        if len(tmp_pos) != 4:
                            sys.exit('length of tmp_pos should be 4')
                        tmp_cos = [x for x in calc_dih(tmp_pos, eps)]
                        cosdih = float(tmp_cos[0])
                        #print('cosdih = ', cosdih)
                        dcosdih = cosdih - maxdih
                        if (dcosdih > eps):
                            maxdih = cosdih
                            kmin = ik
                    new_ind = lindneighbors[kmin]
                else:
                    sys.exit('Number of neighbors should be greater than or equal to 2')

            #print('new_ind = ', new_ind)
            if (new_ind in ind_cycle):
                ind_new_cycle = ind_cycle.index(new_ind)
                tmp_ind_cycle = [i for i in ind_cycle[ind_new_cycle:]]
                ind_cycle = [i for i in tmp_ind_cycle]
                tmp_pos_cycle = [i for i in pos_cycle[ind_new_cycle:]]
                pos_cycle = [i for i in tmp_pos_cycle]
                irepeat = 1
                break
            else:
                ind_cycle.append(new_ind)
                dummy_ind.pop(0)
                dummy_ind.append(new_ind)
                ind_new = lindneighbors.index(new_ind)
                pos_vert_dih.pop(0)
                pos_vert_dih.append(lposneighbors[ind_new])
                pos_cycle.append(lposneighbors[ind_new])

            size_cycle = len(ind_cycle)

        lenindc = len(ind_cycle)
        if lenindc < nindices:
            # find whether any set in the cycle is still part of incomplete set
            for i in range(lenindc):
                if i == (lenindc - 1):
                    ic = 0
                else:
                    ic = i + 1
                trial_pair = [ind_cycle[i-1],ind_cycle[ic]]
                trial_pair.sort()
                #if (iset == 8):
                #    print(incomplete_sets)
                for inset in incomplete_sets:
                    ind_inset = incomplete_sets.index(inset)
                    trial_set = [j for j in inset[1]]
                    trial_set.sort()
                    if (trial_set == trial_pair):
                        print('removed set and its index: ', inset, ',', ind_inset)
                        incomplete_sets.pop(ind_inset)
                        incomplete_indices.pop(ind_inset)
                        break
                #if (iset == 8):
                #    print(incomplete_sets)

            if len(list_cycles2) > 0:
                test_cycle = [i for i in ind_cycle]
                test_cycle.sort()
                icollect = 1
                for cycle in list_cycles2:
                    cycle_tmp = [i for i in cycle]
                    cycle_tmp.sort()
                    if (test_cycle == cycle_tmp):
                        icollect = 0
                        print('***** Trial cycle is a repeatition *****')
                        break
                if (icollect == 1):
                    print('NEW CYCLE: ', ind_cycle)
                    #if (len(ind_cycle) > 6):
                    #    print('POLYGON LARGER THAN HEXAGON: ', ind_cycle)
                    #if (len(ind_cycle) == 5):
                    #    print('PENTAGON: ', ind_cycle)
                    list_cycles2.append([i for i in ind_cycle])
                    pos_tmp = [x for x in pos_cycle]
                    list_pos_cycles2.append(pos_tmp)
            else:
                print('NEW CYCLE: ', ind_cycle)
                list_cycles2.append([i for i in ind_cycle])
                pos_tmp = [x for x in pos_cycle]
                list_pos_cycles2.append(pos_tmp)

        if (irepeat == 1):
            ind_icurr = incomplete_indices.index(icurr_ind)
            print('icurr_ind for repeatition = ', icurr_ind)
            incomplete_indices.pop(ind_icurr)
            incomplete_sets.pop(ind_icurr)
            incomplete_indices.append(icurr_ind)
            incomplete_sets.append(sets)
            #if (iset == 8):
            #    print(incomplete_sets)

numcycles2 = len(list_cycles2)
numtotalcycles = numcycles + numcycles2
print('Total No. of cycles before parsing = ', numtotalcycles)

# Make sure that the cycles in list_cycles2 are not part of list_cycles
if (numcycles2 > 0):
    tmp_list_cycles = [x for x in list_cycles]
    for cycle in list_cycles2:
        #if len(cycle) > 3:
        tmp_cycle = [i for i in cycle]
        tmp_cycle.sort()
        for tcycle in tmp_list_cycles:
            tmp_tcycle = [i for i in tcycle]
            tmp_tcycle.sort()
            if (tmp_tcycle == tmp_cycle):
                cycle_ind = list_cycles2.index(cycle)
                list_cycles2.pop(cycle_ind)
                list_pos_cycles2.pop(cycle_ind)
                break

numcycles2a = len(list_cycles2)
numtotalcycles = numcycles + numcycles2a
print('Total No. of cycles after parsing = ', numtotalcycles)

list_cycles_total = list_cycles + list_cycles2
list_pos_cycles_total = list_pos_cycles + list_pos_cycles2

list_len_cycle = []
for item in list_cycles_total:
    leni = len(item)
    list_len_cycle.append(leni)
list_freq_len = [0]*nindices
for item in list_len_cycle:
    list_freq_len[item-1] = list_freq_len[item-1] + 1

list_len_freq = []
print('len_freq')
ncycles = 0
i = 0
while i < len(list_freq_len):
    freq = list_freq_len[i]
    if freq > 0:
        len_freq = [i+1,freq]
        print(len_freq)
        list_len_freq.append(len_freq)
        ncycles = ncycles + freq
    if (ncycles == numcycles):
        break

    i = i + 1

with open(file3,'a') as fo3:
    for cycle in list_cycles_total:
        lcycle = ['{:6d}'.format(i) for i in cycle]
        string = ''.join(lcycle)
        fo3.write('%s\n' %string)

with open('coords_cycles.txt','a') as fo4:
    for cycle_pos in list_pos_cycles_total:
        #print(cycle_pos)
        string = ''
        for pos in cycle_pos:
            pos_form = ['{:13.6f}'.format(x) for x in pos]
            tmp_string = ''.join(pos_form)
            string = string + tmp_string
        fo4.write('%s\n' %string)
