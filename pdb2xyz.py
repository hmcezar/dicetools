import argparse as arg

# fmt: off
ELEMENTS = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12,
            'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23,
            'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
            'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45,
            'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56,
            'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
            'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
            'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89,
            'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
            'Md': 101, 'No': 102, 'Lr': 103, 'Xx': 104}   # yapf: disable
# fmt: on


DESCRIPTION = """***********************************************************************

Convert GROMACS trajectory file in pdb format to a DICE xyz file format

Last modification: 06/Mar/2021

Developed by Emanuel Mancio                

Supervised by Prof Kaline Coutinho

Instituto de Fisica da Universidade de Sao Paulo (IF/USP)

***********************************************************************"""

STEP = 0


class IncorrectNumberOfAtomsOnTrajectory(Exception):
    """ Exception raised when the number of atoms provided differs from the trajectory """

    def __init__(
        self,
        message="The number of atoms on the trajectory is different than the provided",
    ):

        super().__init__(message)


class IncorrectNumberOfAtomsOnTopology(Exception):
    """ Exception raised when the number of atoms provided differs from the topology """

    def __init__(self, mol, expected, given):
        """Constructor of the exception IncorrectNumberOfAtomsOnTopology

        Args:
            mol (int): txt file molecule`s number
            expected (int): expected number of atoms of the molecule (mol)
            given ([type]): quantity of atoms of the molecule (mol) in the txt file 
        """
        message = f"""The number of atoms given ({given}) of the molecule {mol} on the topology 
        (txt) is different than the expected {expected}"""
        super().__init__(message)


def process_config(pdbfile, xyzfile, natoms, elements, reset_step=False):
    """Process a configuration of the trajectory pdb file

    Args:
        pdbfile (arg.FileType): pdb input file opened by the argparse module
        xyzfile (arg.FileType): xyz output file created by the argparse module
        natoms (int): number of atoms in the configuration
        elements (list): list with the elements in the same order from the txt file
        reset_step (bool, optional): if False the configuration number is the same from the pdb file.
                                     Otherwise the number will start from 0 and have a step of 1. Defaults to False.

    Raises:
        IncorrectNumberOfAtomsOnTrajectory: raised when the number of atoms provided differs from what there is in 
                                            the configuration

    Returns:
        int: 1 if the file ended or 0 otherwise
    """
    if pdbfile.readline() == "":
        return 1  # end of file

    global STEP

    xyzfile.write("{0: >12}\n".format(natoms))

    conf_info = " Configuration number :{0: >9}L = {1:>9.4f}{2:>9.4f}{3:>9.4f}\n"

    if not reset_step:
        step = (pdbfile.readline().split())[-1]

        pdbfile.readline()  # Pass the REMARK line

        STEP += 1
        try:
            dims = list(map(float, (pdbfile.readline().split())[1:4]))
        except:
            raise IncorrectNumberOfAtomsOnTrajectory()

        xyzfile.write(conf_info.format(step, dims[0], dims[1], dims[2]))

    else:
        pdbfile.readline()  # Pass the TITLE line
        pdbfile.readline()  # Pass the REMARK line

        STEP += 1
        try:
            dims = list(map(float, (pdbfile.readline().split())[1:4]))
        except:
            raise IncorrectNumberOfAtomsOnTrajectory()
        xyzfile.write(conf_info.format(STEP, dims[0], dims[1], dims[2]))

    pdbfile.readline()  # MODEL line

    dims[0], dims[1], dims[2] = dims[0] / 2, dims[1] / 2, dims[2] / 2

    pos_info = "  {0: >2}{1:>-15.5f}{2:>-15.5f}{3:>-15.5f}\n"

    for i in range(natoms):
        pos = list(map(float, (pdbfile.readline().split())[5:8]))

        try:
            xyzfile.write(
                pos_info.format(
                    elements[i], pos[0] - dims[0], pos[1] - dims[1], pos[2] - dims[2]
                )
            )
        except IndexError:
            raise IncorrectNumberOfAtomsOnTrajectory()

    pdbfile.readline()  # Pass ENDML line
    pdbfile.readline()  # Pass TER line

    return 0 # if the file didn`t end


def pass_config(pdbfile, natoms):
    """Used to ignore a configuration

    Args:
        pdbfile (arg.FileType): pdb input file opened by the argparse module
        natoms (int): number of atoms in the configuration
    """
    for _ in range(natoms + 7): # 7 is the number of lines with comments in the pdb file
        pdbfile.readline()


def get_key(val):
    """Gets the key from a dictionary from the value that the key holds

    Args:
        val (any): value held by the key

    Returns:
        any: key from the dictionary that hold the value (val)
    """
    global ELEMENTS

    for key, value in ELEMENTS.items():
        if val == value:
            return key


def get_elements(nmol, at_per_num):
    """Generate a list with the symbols of the atoms in the order of
       the txt file and in the qunatity provided

    Args:
        nmol (int): number of molecules
        at_num (list): list of list with the atomic number of the atoms of the molecules

    Returns:
        list: list with the atoms symbols
    """
    global ELEMENTS

    el = []

    for n, mol in zip(nmol, at_per_num):
        el_mol = []
        for num in mol:
            el_mol.append(get_key(num))

        el = el + el_mol * n

    return el


def read_txt(txtfile):
    """Function that read the txt file

    Args:
        txtfile (arg.FileType): txt input file opened by the argparse module

    Raises:
        IncorrectNumberOfAtomsOnTopology: when the number of atoms atoms that
        the molecule is supposed to be is different than the actual number

    Returns:
        list: list of lists of the atomic numbers of the atoms in the molecules
    """
    with txtfile as txt:
        txt.readline()
        nmol = int(txt.readline())

        at_per_mol = []

        for mol in range(nmol):
            natom = int((txt.readline().split())[0])
            at_num = []

            try:
                for at in range(natom):
                    num = int((txt.readline().split())[1])
                    at_num.append(num)

                at_per_mol.append(at_num)
            except IndexError:
                raise IncorrectNumberOfAtomsOnTopology(mol + 1, natom, at)

    return at_per_mol


def check_txt(txt, nmol):
    """Checks the txt file to have the same quantity of molecules as the providedd in the command line

    Args:
        txt (arg.FileType): txt input file opened by the argparse module
        nmol (int): quantity of molecules provided in the command line

    Raises:
        ValueError: when the quantity of molecules is different tahn the provided
    """
    txt.readline()
    n = int(txt.readline())
    txt.seek(0) # return to beginning of the file

    msg = " ".join(map(str, nmol))

    if n > len(nmol):
        raise ValueError(
            msg,
            "Not enough number of molecules were provided, expected {}, given {}".format(
                n, len(nmol)
            ),
        )
    elif n < len(nmol):
        raise ValueError(
            msg,
            "More than required number of molecules were provided, expected {}, given {}".format(
                n, len(nmol)
            ),
        )


def check_slice(start, final, intv, printinterval):
    """Check if the arguments of init, final or printinterval (ou intv) are valid

    Args:
        start (int): argument of start
        final (int): argument of final
        intv (int): argument of intv
        printinterval (int): argument of printinterval

    Raises:
        ValueError: when the arguments are invalid
    """
    if start < 1:
        raise ValueError(start, "The start argument should be more than or equal to 1")
    elif final < 1:
        raise ValueError(final, "The final argument should be more than or equal to 1")
    elif final < start:
        raise ValueError(
            final, "The final argument should be should be bigger than start"
        )
    elif intv < 1:
        raise ValueError(intv, "The intv argument should be more than or equal to 1")
    elif printinterval < 1:
        raise ValueError(
            printinterval,
            "The printinterval argument should be more than or equal to 1",
        )


if __name__ == "__main__":
    parser = arg.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=arg.RawDescriptionHelpFormatter,
        add_help=False,
    )

    grp_input = parser.add_argument_group("Input arguments")

    grp_output = parser.add_argument_group(
        "Output arguments",
        description="""The init, final, intv or printinterval should receive the frame number
        not the time step (only accepts integer positive numbers) """,
    )

    information = parser.add_argument_group("Informational arguments")

    grp_input.add_argument(
        "pdb",
        type=arg.FileType("r"),
        help="name of pdb trajectory file generated by GROMACS trjconv",
    )

    grp_input.add_argument(
        "txt", type=arg.FileType("r"), help="name of txt DICE topology file"
    )

    grp_input.add_argument(
        "nmol",
        type=int,
        nargs="+",
        help="""quantity of molecules of each type in the same order 
        as the txt file. Example: 1 1000 (if the system has 2 types)""",
    )

    grp_output.add_argument(
        "-o",
        type=arg.FileType("w+"),
        default="output.xyz",
        help="name of xyz trajectory file in DICE format (default: output.xyz)",
        metavar="name",
    )

    grp_output.add_argument(
        "-r", action="store_true", help="renumbers the frames of the trajectory from 0"
    )

    grp_output.add_argument(
        "-init",
        type=int,
        nargs=1,
        default=[1],
        help="first frame to be printed (default: 1)",
        metavar="#",
    )

    grp_output.add_argument(
        "-final",
        type=int,
        nargs=1,
        default=[-1],
        help="last frame to be considered (default: last config in pdb file)",
        metavar="#",
    )

    grp_output.add_argument(
        "-intv",
        type=int,
        nargs=1,
        default=[1],
        help="interval between the frames to be printed (default: 1)",
        metavar="#",
    )

    grp_output.add_argument(
        "-printinterval", type=int, nargs=1, default=[1], help="or -intv", metavar="#"
    )

    information.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )
    information.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0"
    )

    args = parser.parse_args()

    check_txt(args.txt, args.nmol)

    at_per_mol = read_txt(args.txt)


    el = get_elements(args.nmol, at_per_mol)

    natoms = len(el)

    start, final, intv, printinterval = (
        args.init[0],
        args.final[0],
        args.intv[0],
        args.printinterval[0],
    )

    if final == -1:
        final = 2147483647 # giant number to be considered as final if no argument is given

    check_slice(start, final, intv, printinterval)

    if intv < printinterval: # if both intv and printinterval, the higher will be used
        intv = printinterval

    for _ in range(start - 1): # when start is different than 1
        STEP += 1
        pass_config(args.pdb, natoms)

    if intv == 1:
        while True:
            end_of_file = process_config(args.pdb, args.o, natoms, el, args.r)
            if end_of_file or STEP >= final:
                break

    else:
        while True:
            end_of_file = process_config(args.pdb, args.o, natoms, el, args.r)

            for _ in range(intv - 1):
                STEP += 1
                pass_config(args.pdb, natoms)

            if end_of_file or STEP >= final:
                break

    args.o.close()
    args.pdb.close()
