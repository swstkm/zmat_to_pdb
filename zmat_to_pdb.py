'''
This code converts zmat to pdb
usage: zmat_to_pdb.py <zmat-file>
'''

# library for functions used for calculations
import numpy as np
# library for using argv function later, to read filepath of zmat file
import sys


# this functions reads the z-matrix file and
# returns the information in lists
def zmatr(filename):
    with open(filename, 'r') as zmatf:
        # read teh file object into a list with the same variable name
        zmatf = zmatf.readlines()
        # each of these is a list that will store
        # each of the columns in z matrix file
        anamelist = []
        rdto = []
        rdist = []
        ang_to = []
        angleslist = []
        dih_to = []
        dihlist = []

        for smline in zmatf:
            # for each line in zmatf check what the length of the list is
            # when the list is split based on space characters
            # then append the values to the correct list defined earlier,
            # based on the index of the element in the list
            splitline = smline.split()
            if len(splitline) > 0:
                anamelist.append(splitline[0])
            if len(splitline) > 1:
                rdto.append(int(splitline[1]))
            if len(splitline) > 2:
                rdist.append(float(splitline[2]))
            if len(splitline) > 3:
                ang_to.append(int(splitline[3]))
            if len(splitline) > 4:
                angleslist.append(float(splitline[4]))
            if len(splitline) > 5:
                dih_to.append(int(splitline[5]))
            if len(splitline) > 6:
                dihlist.append(float(splitline[6]))
    return (anamelist, rdto, rdist, ang_to,
            angleslist, dih_to, dihlist, filename)


# this functions computes the coordinates and writes out the pdb file
# it takes as input, the tuple that the zmatr function returns
def pdbr(anamelist, rdto, rdist, ang_to,
         angleslist, dih_to, dihlist, filename):
    # length of anamelist is the total number of atoms
    natoms = len(anamelist)
    # define output file name
    outfilename = filename.split('\\').pop().rsplit('.', 1)[0]

    # header line that says what compound this pdb represents
    compound = (filename.split('/')[-1]).split('.')[0]

    # with output file open in write mode, write the following
    with open(outfilename + '.pdb', 'w') as outf:
        # write compound and authors name
        outf.write("COMPND    " + compound.upper() +
                   "\nAUTHOR    SWASTIK MISHRA\n")

        # initialise coordinates array using numpy.zeros
        coor = np.zeros([natoms, 3])
        # this is the first atom placed at the origin
        if (natoms > 1):
            coor[1] = [rdist[0], 0.0, 0.0]

        # atom connected to previous other atom but no angle dependency
        if (natoms > 2):
            # init indices
            # index of atom to which distance is defined to
            i = rdto[1] - 1
            # index of atom to which angle is defined to
            j = ang_to[0] - 1
            r = rdist[1]             # distance from origin
            # find angle for this atom
            angq = angleslist[0] * np.pi / 180.0
            # find x and y coordinates based on projection along x and y axes
            x = r * np.cos(angq)
            y = r * np.sin(angq)
            a_i = coor[i]
            b_ij = coor[j] - coor[i]
            if (b_ij[0] < 0):
                x = a_i[0] - x
                y = a_i[1] - y
            else:
                x = a_i[0] + x
                y = a_i[1] + y
            coor[2] = [x, y, 0.0]
        # for all other atoms
        for n in range(3, natoms):
            # just like previous case,
            # find distance, angle and dihedral has been
            # defined with respect to indices of which atom
            r = rdist[n - 1]
            angq = angleslist[n - 2] * np.pi / 180.0
            dihangle = dihlist[n - 3] * np.pi / 180.0

            # find sin and cos for angles and dihedrals
            # and use those to find x,y,z, with respect
            # to this coordinate system
            sinPhi = np.sin(dihangle)
            cosPhi = np.cos(dihangle)
            sinTheta = np.sin(angq)
            cosTheta = np.cos(angq)

            x = r * cosTheta
            y = r * cosPhi * sinTheta
            z = r * sinPhi * sinTheta

            # transform coordinates based on previously defined coordinates
            i = rdto[n - 1] - 1
            j = ang_to[n - 2] - 1
            k = dih_to[n - 3] - 1

            a = coor[k]
            b = coor[j]
            c = coor[i]

            ab = b - a
            bc = c - b
            bc = bc / np.linalg.norm(bc)
            nv = np.cross(ab, bc)
            nv = nv / np.linalg.norm(nv)
            ncbc = np.cross(nv, bc)

            new_x = c[0] - bc[0] * x + ncbc[0] * y + nv[0] * z
            new_y = c[1] - bc[1] * x + ncbc[1] * y + nv[1] * z
            new_z = c[2] - bc[2] * x + ncbc[2] * y + nv[2] * z
            coor[n] = [new_x, new_y, new_z]

        # write out the coordinates in pdb format
        for i in range(natoms):
            outf.write('{:<6}{:>5}  {:<4}{:>3} {:^1}{:>4}    {:>8}{:>8}{:>8}{:>6}{:>6}          {:>2}'.format('HETATM',
              str(i + 1), anamelist[i].upper(), 'NIL', 'X', '1', str(
                  round(coor[i][0], 3)),
              str(round(coor[i][1], 3)), str(round(coor[i][2], 3)), '1.00', '0.00', anamelist[i].upper()))
            outf.write('\n')

        '''for CONECT statements
        the following works only for zmat files
        where the distance in third column
        is measured for bonds only,
        which should be mentioned in the second column'''

        # pad rdto with a 0 in the beginning
        rdto.insert(0, 0)
        # this is a dictionary that will contain index of atom as key
        # and atoms it is connected to in a list, which will be the value
        connect = {x: [] for x in range(1, natoms + 1)}
        # if there is an atom distance defined with respect to another
        # add both to each others' dictionary lists
        for i in connect.keys():
            for j in range(len(rdto)):
                if rdto[j] == i:
                    connect[i].append(str(j + 1))
                    connect[j + 1].append(str(i))
            if connect[i] != []:
                # if atom is connected to anything at all
                # write out a CONECT statement using the dictionary
                outf.write('{:<6}{:>5}'.format(
                    'CONECT', str(i)))
                for k in range(len(connect[i])):
                    outf.write('{:>5}'.format((connect[i][k])))
            outf.write('\n')


############################################

# call all the functions appropriately
# read the filepath
zmatfilename = sys.argv[1]
# use filepath to read out values for zmatrix file
anamelist, rdto, rdist, ang_to, angleslist, dih_to, dihlist, filename = zmatr(
    zmatfilename)
# write out pdb file after doing computation using the zmatrix info
pdbr(anamelist, rdto, rdist, ang_to,
     angleslist, dih_to, dihlist, filename)
