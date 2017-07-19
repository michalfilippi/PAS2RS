import sys

import numpy as np


def central_carbon_planes_deflection(coordinates):
    """Calculates torsion angles and torsion directions. For every residue at 
    position i, a angle between planes p and q is calculated, where p is defined
    by residues at position i - 2, i - 1 and i, q is defined by residues at 
    position i - 1, i and i + 1.
    
    :param coordinates: List of absolute positions. Every position is tuple.
    [x,y,z].
    :return: Tuple (angles, directions), where angles and directions are lists 
    of length len(coordinates) containing torsion angles and torsion directions.
    """

    angles = [0, 0]
    directions = [0, 0]

    for i in range(2, len(coordinates) - 1):

        # calculate vectors defining planes p and q
        u = np.subtract(coordinates[i - 1], coordinates[i - 2])
        v = np.subtract(coordinates[i], coordinates[i - 1])
        w = np.subtract(coordinates[i + 1], coordinates[i])

        # plane normal vectors
        n_p = np.cross(u, v)
        n_q = np.cross(v, w)

        # if v and w have the same direction, therefore torsion angle is 0
        if np.linalg.norm(n_q) == 0:
            angles.append(0)
            directions.append(0)
            continue

        # if u and v has the same direction, look at the previous plane
        j = i
        while np.linalg.norm(n_p) == 0:
            j -= 1
            if j == 1:
                # no previous plane
                break

            # redefine plane p
            new_u = np.subtract(coordinates[j - 1], coordinates[j - 2])
            new_v = np.subtract(coordinates[j], coordinates[j - 1])

            # new normal vector
            n_p = np.cross(new_u, new_v)

        # no previos plane found, therefore torsion angle is 0
        if np.linalg.norm(n_p) == 0:
            angles.append(0)
            directions.append(0)
            continue

        # calculate the torsion angle
        cos_a = np.dot(n_p, n_q) / (np.linalg.norm(n_p) * np.linalg.norm(n_q))

        #
        cos_a = np.round(cos_a, decimals=5)

        # scaled torsion angle (0, 180) ~ (0, 1)
        angle = np.arccos(cos_a) / np.pi

        # plane p: ax + by + cz + d = 0
        # n_p = [a, b, c]
        # calculate d
        s = coordinates[i]
        t = coordinates[i + 1]
        d = - n_p[0] * s[0] - n_p[1] * s[1] - n_p[2] * s[2]

        # calculate the torsion direction
        direc = np.sign(n_p[0] * t[0] + n_p[1] * t[1] + n_p[2] * t[2] + d)

        angles.append(angle)
        directions.append(direc)

    angles.append(0)
    directions.append(0)

    return angles, directions


def central_carbon_angles(coordinates):
    """Calculates angles between every two neighboring residues. 
    
    :param coordinates:  List of absolute positions. Every position is tuple 
    [x,y,z].
    :return: List of angles of length len(coordinates).
    """

    cc_angles = [0]

    for i in range(1, len(coordinates) - 1):
        u = np.subtract(coordinates[i], coordinates[i - 1])
        v = np.subtract(coordinates[i + 1], coordinates[i])
        cos_uv = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))

        # scale angle (0, 180) ~ (0, 1)
        angle = np.arccos(cos_uv) / np.pi

        cc_angles.append(angle)

    cc_angles.append(0)

    return cc_angles


def transformPAS2RS(coordinates):
    """Transforms protein absolute structure to protein relative structure. 
    Relative structure is represented by 5 features for every residue in 
    protein. These are central carbon angle, forward  torsion angle and 
    direction and backward torsion angle and direction. 
    
    :param coordinates: List of absolute positions. Every position is tuple 
    [x,y,z].
    :return: Relative structure as list of length len(coordinates).
    """

    cc_angles = central_carbon_angles(coordinates)

    # carbon planes deflection and direction of deflection in forward direction
    ccp_a_f, ccp_d_f = central_carbon_planes_deflection(coordinates)

    # carbon planes deflection and direction of deflection in backward direction
    ccp_a_b, ccp_d_b = central_carbon_planes_deflection(coordinates[::-1])
    ccp_a_b, ccp_d_b = ccp_a_b[::-1], ccp_d_b[::-1]

    # merge all residue features together
    relative_structure = zip(cc_angles, ccp_a_f, ccp_d_f, ccp_a_b, ccp_d_b)

    return list(relative_structure)


def main():
    """Main function to be run when module is not used as library. Protein 
    absolute structure, i.e. list of central carbon coordinates is read from
    input file passed as sys.argv[1]
    
    :return: None.
    """

    # get the input file name
    try:
        coordinates_file = sys.argv[1]
    except IndexError:
        print('Missing parameter.\n')
        print('Usage: python3 protstr_converter.py coordinates_file\n')
        return

    # read absolute structure from 'coordinates_file'
    coordinates = np.loadtxt(coordinates_file)

    # transform structure
    relative_structure = transformPAS2RS(coordinates)

    # print relative structure to stdout
    print('CCA CCPTA_F CCPTD_F CCPTA_B CCPTD_B')
    for rs in relative_structure:
        print(*rs)


if __name__ == "__main__":
    main()
