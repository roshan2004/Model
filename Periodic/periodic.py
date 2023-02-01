import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_angles, calc_bonds, calc_dihedrals
import argparse


# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", type=str, default='martini_graphene',
                    help='Name of the output, default = martini_graphene')

args = parser.parse_args()

filename = args.filename
rows = 6
columns = 6
positions = []
dist = 0
for i in range(rows):
    if i % 2 == 0:
        x = np.linspace(0, (0+(columns-1)*2.56), columns)
        for k in range(len(x)):
            positions.append([x[k], dist, 0])
        dist += 2.217

    else:
        x = np.linspace(-1.28, (-1.28+(columns-1)*2.56), columns)
        for k in range(len(x)):
            positions.append([x[k], dist, 0])
        dist += 2.217
w = mda.Universe.empty(n_atoms=len(positions), trajectory=True)
w.atoms.positions = positions
w.add_TopologyAttr('names')
w.atoms.names = [f'B{i}' for i in range(1, w.atoms.n_atoms + 1)]
w.add_TopologyAttr('resnames')
w.residues.resnames = ['GRA']
w.atoms.write(filename+'.gro')


u = mda.Universe(filename+'.gro')
c = np.unique(u.atoms.positions[:, 1])

u.atoms.masses = 36



# Setting mass of the virtual-site 0
for j in range(len(c)):
    if j == 0:
        group = u.atoms[u.atoms.positions[:, 1] == c[j]]
        idx = group.atoms.indices
        gr = np.arange(idx[2], idx[-1]+1, 3)
        for k in gr:
            u.atoms[k].mass = 0
    elif (j != len(c) - 1) and (j % 2 != 0) and (j != 0):
        group = u.atoms[u.atoms.positions[:, 1] == c[j]]
        gr = np.arange(1, len(group), 3)

        for k in gr:
            u.atoms[group[k].index].mass = 0
    elif (j != len(c) - 1) and (j % 2 == 0) and (j != 0):
        group = u.atoms[u.atoms.positions[:, 1] == c[j]]
        gr = np.arange(2, len(group), 3)

        for k in gr:
            u.atoms[group[k].index].mass = 0

    elif j == len(c) - 1:
        group = u.atoms[u.atoms.positions[:, 1] == c[j]]
        idx = group.atoms.indices
        gr = np.arange(idx[1], idx[-1]+1, 3)
        for k in gr:
            u.atoms[k].mass = 0


def regroup_lower(lst):
    final = []

    threes = [lst[i:i+3] for i in range(0, len(lst), 3)]
    for l, r in zip(threes, threes[1:]):
        final.append(l)
        final.append([l[-1], r[0]])
    if len(threes[-1]) > 1:
        final.append(threes[-1])

    data = []
    for i in final:
        if len(i) == 3:
            data.append(i)
        else:
            data.append(i)
    return data

def regroup_upper(lst):
    sublists = []
    data = []

    for idx in range(0, len(lst), 3):
        pair = lst[idx:idx+2]
        triplet = lst[idx+1:idx+4]

        if len(pair) == 2:
            sublists.append(pair)

        if len(triplet) == 3:
            sublists.append(triplet)
    for i in sublists:
        if len(i) == 3:
            data.append(i)
        else:
            data.append(i)
    return data

def hexagon(universe):
    idx_top = u.atoms[u.atoms.positions[:, 1] == c[0]]
    idx_bottom = u.atoms[u.atoms.positions[:, 1] == c[-1]]

    sel = u.atoms[(u.atoms.masses == 0)]
    sel1 = sel.atoms[sel.atoms.positions[:, 0]
                     == np.max(sel.atoms.positions[:, 0])]
    idx = u.atoms - idx_top - idx_bottom - sel1

    b = idx.atoms[idx.atoms.masses == 0].indices

    hexagon_indices = []
    for i in b:
        empty = []
        for j in u.atoms.indices:
            if (2.55 <= calc_bonds(u.atoms[j].position, u.atoms[i].position) <= 2.57):
                empty.append(j)
        empty.append(i)
        hexagon_indices.append(empty)

    '''
    Hexagons along the horizontal, or the x-axis
    '''
    upper = u.atoms[u.atoms.positions[:, 1] == c[0]]
    lower = u.atoms[u.atoms.positions[:, 1] == c[-1]]
    upper_down = u.atoms[u.atoms.positions[:, 1] == c[1]]
    lower_up = u.atoms[u.atoms.positions[:, 1] == c[-2]]

    groups_upper = regroup_upper(upper.atoms.indices)
    groups_lower = regroup_lower(lower.atoms.indices)
    groups_upper_down = regroup_lower(upper_down.atoms.indices)
    groups_lower_up = regroup_upper(lower_up.atoms.indices)
    for i in range(len(groups_upper)):
        if i % 2 == 0:
            hexagon_indices.append([groups_lower_up[i][0], groups_lower_up[i][1], groups_lower[i]
                                   [0], groups_lower[i][2], groups_upper[i][0], groups_upper[i][1], groups_lower[i][1]])

        else:
            hexagon_indices.append([groups_lower[i][0], groups_lower[i][1], groups_upper[i]
                                   [0], groups_upper[i][2], groups_upper_down[i][0], groups_upper_down[i][1], groups_upper[i][1]])
            
    return hexagon_indices
    


vs = hexagon(u)
hexagons = [hex[:-1] for hex in vs]
def virtual_sites(vs):
    virtual_indices = []
    for hex in vs:
        virtual_indices.append([hex[-1], hex[0], hex[1], hex[4], hex[5]])
    return virtual_indices

virtual_sites(vs)








def bonds(hexagons):
    bonds = []
    for hex in hexagons:
        bond1 = sorted([hex[0], hex[1]])
        bond2 = sorted([hex[0], hex[2]])
        bond3 = sorted([hex[2], hex[4]])
        bond4 = sorted([hex[4], hex[5]])
        bond5 = sorted([hex[5], hex[3]])
        bond6 = sorted([hex[1], hex[3]])

        if bond1 not in bonds:
            bonds.append(bond1)
        if bond2 not in bonds:
            bonds.append(bond2)
        if bond3 not in bonds:
            bonds.append(bond3)
        if bond4 not in bonds:
            bonds.append(bond4)
        if bond5 not in bonds:
            bonds.append(bond5)
        if bond6 not in bonds:
            bonds.append(bond6)

    return bonds


def find_angles(hexagons):
    angles = []
    for hex in hexagons:
        angle1 = [hex[0], hex[2], hex[4]]
        angle2 = [hex[2], hex[4], hex[5]]
        angle3 = [hex[4], hex[5], hex[3]]
        angle4 = [hex[5], hex[3], hex[1]]
        angle5 = [hex[3], hex[1], hex[0]]
        angle6 = [hex[1], hex[0], hex[2]]
        angles.append(angle1)
        angles.append(angle2)
        angles.append(angle3)
        angles.append(angle4)
        angles.append(angle5)
        angles.append(angle6)
    return angles


def get_angles(hexagons):
    angles = []
    for hex in hexagons:
        angle1 = [hex[2], hex[4], hex[0], hex[5]]
        angle2 = [hex[4], hex[0], hex[5], hex[1]]
        angle3 = [hex[0], hex[5], hex[1], hex[3]]
        angles.append(angle1)
        angles.append(angle2)
        angles.append(angle3)
    return angles


def get_between_hexagons_angles(hexagons, beads_per_row):
    def get_common_side(hexa1, hexa2):
        common = list(set(hexa1).intersection(set(hexa2)))
        common.sort()
        return common

    angles = []
    treated_side = []
    gap_for_after = int(beads_per_row / 3)
    gap_for_below = int(beads_per_row / 3) + int(beads_per_row/3 - 1)
    gap_for_before = int(beads_per_row / 3) - 1

    for i in range(len(hexagons)):
        hexa = hexagons[i]
        try:
            hexa_after = hexagons[i+gap_for_after]
        except IndexError:
            hexa_after = None
        try:
            hexa_below = hexagons[i+gap_for_below]
        except IndexError:
            hexa_below = None

        if i > 0:
            try:
                hexa_before = hexagons[i+gap_for_before]
            except IndexError:
                hexa_before = None
        else:
            hexa_before = None

        if hexa_below == hexa_before:
            hexa_before = None

        if hexa_after:
            side_after = get_common_side(hexa, hexa_after)
            if side_after:
                idx1 = hexa.index(side_after[0])
                idx2 = hexa_after.index(side_after[1])
                side = (side_after[0], side_after[1])
                if side not in treated_side:
                    treated_side.append(side)
                    angle = [hexa[idx1 - 2], side_after[0],
                             side_after[1], hexa_after[idx2 + 2]]
                    angles.append(angle)

        if hexa_below:
            side_below = get_common_side(hexa, hexa_below)
            if side_below:
                idx1 = hexa.index(side_below[0])
                idx2 = hexa_below.index(side_below[1])
                side = (side_below[0], side_below[1])
                if side not in treated_side:
                    treated_side.append(side)
                    angle = [hexa[idx1 - 4], side_below[0],
                             side_below[1], hexa_below[idx2 + 4]]
                    angles.append(angle)

        if hexa_before:
            side_before = get_common_side(hexa, hexa_before)
            if side_before:
                idx1 = hexa.index(side_before[0])
                idx2 = hexa_before.index(side_before[1])
                side = (side_before[0], side_before[1])
                if side not in treated_side:
                    treated_side.append(side)
                    angle = [hexa[idx1 - 2], side_before[0],
                             side_before[1], hexa_before[idx2 + 2]]
                    angles.append(angle)

    return angles


impropers = np.vstack(
    (get_angles(hexagons), get_between_hexagons_angles(hexagons, int(columns))))
impropers = impropers + 1


def exclusions(vs):
    return np.array(vs) + 1



#---------------#
# Topology File #
#---------------#


# Open the file for writing

topology_file = open(filename+".itp", 'w')

# Variables


# Header

topology_file.write(
    "; \n;  Graphene topology\n; for the Martini3 force field\n;\n; created by martini3-graphene-topology.py\n;\n")
topology_file.write("; Roshan Shrestha\n; CNRS\n;\n\n")

topology_file.write("[ moleculetype ]\n")
topology_file.write("; molname	 nrexcl\n")
topology_file.write("  GRA           1")


# Atoms

topology_file.write("\n[ atoms ]\n")
topology_file.write("; nr	 type	 resnr	 residue	 atom	 cgnr	 charge	 mass\n")
for i in range(1, u.atoms.n_atoms+1):
    topology_file.write(
        f"  {i:<5}     TC5     0     GRA     B{i:<5}     {i:<5}     0     {int(u.atoms[i-1].mass)}\n")

# Atoms

topology_file.write("\n[ atoms ]\n")
topology_file.write("; nr	 type	 resnr	 residue	 atom	 cgnr	 charge	 mass\n")
for i in range(1, u.atoms.n_atoms+1):
    topology_file.write(
        f"  {i:<5}     TC5     0     GRA     B{i:<5}     {i:<5}     0     {int(u.atoms[i-1].mass)}\n")

# Bonds

topology_file.write("\n[ bonds ]\n")
topology_file.write("; i	 j	  funct	 length	 kb\n")
for i in bonds(hexagons):
    topology_file.write(
        f"  {i[0]+1:<3}     {i[1]+1:<3}     1    0.24595     25000\n")
    




# Angles

topology_file.write("\n[ angles ]\n")
topology_file.write("; i	 j	 k	 funct	 angle	 force_k\n")
for i in find_angles(hexagons):
    topology_file.write(
        f"  {i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}     1    120     300\n")


# Improper Dihedrals

topology_file.write("\n[ dihedrals ]\n")
topology_file.write("; i	 j	 k	 l     funct	 ref.angle     force_k\n")
for i in impropers:
    topology_file.write(
        f"  {i[0]:<3}     {i[1]:<3}     {i[2]:<3}     {i[3]:<3}    2     180     200\n")


# Virtual sites

topology_file.write("\n[ virtual_sitesn ]\n")
topology_file.write("; site	 funct	 constructing atom indices\n")
for i in virtual_sites(u):
    topology_file.write(
        f"  {i[0]+1:<3}     1     {i[1]+1:<3}     {i[2]+1:<3}     {i[3]+1:<3}    {i[4]+1:<3}\n")

# Exclusions
topology_file.write("\n[ exclusions ]\n")
for i in exclusions(u):
    topology_file.write(
        f"  {i[0]+1:<3}     {i[1]+1:<3}     {i[2]+1:<3}     {i[3]+1:<3}    {i[4]+1:<3}     {i[5]+1:<3}     {i[6]+1:<3}\n")

topology_file.close()



