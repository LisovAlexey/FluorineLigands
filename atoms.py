import re
import typing as tp
from pathlib import Path


match = re.compile("([A-z]*)(\d*)")


class Atom:
    def __init__(self, id_, str_id, x, y, z, atype, unk1, mol_name, unk2):
        self.type = None
        self.type_id = None

        self.id: int = id_
        self.str_id: str = str_id
        self.x: str = x
        self.y: str = y
        self.z: str = z
        self.atype: str = atype
        self.unk1: str = unk1
        self.mol_name: str = mol_name
        self.unk2: str = unk2

    def __post_init__(self):
        self.id = int(self.id)

    @property
    def str_id(self) -> str:
        return self.type + str(self.type_id)

    @str_id.setter
    def str_id(self, value: str) -> None:
        self.type, self.type_id = match.findall(value)[0]
        self.type_id = int(self.type_id)


class Mol2Atoms:
    def __init__(self):
        self.dict: tp.Dict[int, Atom] = OrderedDict()
        # self.atom_chains = {}

    def __get__(self, id_):
        pass

    def __set__(self, id_):
        pass

    def _replace_atom(self, atom_id: int, new_atom_type: str):
        """
        atom id: starts from 1
        if atom_form is not None
        """

        atom = self.dict[atom_id]
        replaced_atom_type = atom.type

        type_id = 1

        for atom_id_, atom in self.dict.items():
            if atom_id_ < atom_id:
                if atom.type == new_atom_type:
                    type_id += 1
            elif atom_id_ > atom_id:
                if atom.type == replaced_atom_type:
                    atom.type_id -= 1
                elif atom.type == new_atom_type:
                    atom.type_id += 1
            else:
                atom.type_id = type_id
                atom.type = new_atom_type
                atom.atype = new_atom_type


    @classmethod
    def replace_atom(cls, obj, atom_id, new_atom_type):
        obj_replaced = deepcopy(obj)
        obj_replaced._replace_atom(atom_id, new_atom_type)

        return obj_replaced

    @classmethod
    def from_strlist(cls, slist: tp.List[str]):
        obj = cls()

        for atom_str in slist[1:]:
            atom = Atom(*atom_str.split())
            obj.dict[int(atom.id)] = (atom)

        return obj

    def __str__(self):
        l = ["@<TRIPOS>ATOM"]
        for k, atom in self.dict.items():
            s = f"{str(k).rjust(8)} {atom.str_id.rjust(6)} {atom.x.rjust(12)} {atom.y.rjust(12)} {atom.z.rjust(12)} \
            {atom.atype.rjust(8)} {atom.unk1.rjust(3)} {atom.mol_name.rjust(16)} {atom.unk2.rjust(10)}"
            l.append(s)

        return "\n".join(l) + "\n"


class Mol2TriposMolecule:
    def __init__(self):
        self.header = None
        self.atoms: tp.Optional[Mol2Atoms] = None
        self.bonds = None

    def to_mol2(self, path: Path) -> None:
#         print("@@@@", str(self.bonds))
        with path.open("w") as wf:
            wf.write(str(self.header))
            wf.write(str(self.atoms))
            wf.write(str(self.bonds))

    def generate_atom_replacements(self, from_atom, to_atom, rep_lens = None) -> tp.List:
        """
        param: rep_lens tp.List[int] - lengths of replacements (number of atoms replaced in one replacement)
        """
            
        # Get ids of atoms that should be replaced
        atom_ids = []
        for _, atom in self.atoms.dict.items():
            if atom.type == from_atom:
#                 print("Detected atom with type: ", atom.type)
#                 print(atom)
                atom_ids.append(int(atom.id))
    
        if rep_lens is None:
            rep_lens = range(1, len(atom_ids) + 1)

#         print("Full list of atoms: ", atom_ids)

        # Get combinations of atoms
        combs = []
        for i in rep_lens:
            combs_i = list(combinations(atom_ids, i))
            combs.extend(combs_i)

#         print(f"Full list of combinations: ", combs)

        # replace atoms
        replacements = []
        for comb in combs:
            atom_copy = deepcopy(self.atoms)
            for atom_id in comb:
                atom_copy._replace_atom(atom_id, to_atom)
            replacements.append(atom_copy)

        return replacements

    @classmethod
    def from_mol2(cls, path: Path):
        obj = cls()

        with path.open("r") as rf:
            read_state = "begin"

            mol_raw = []
            atom_raw = []
            bond_raw = []

            for line in rf:
                if line.startswith("@<TRIPOS>MOLECULE"):
                    read_state = "mol"
                elif line.startswith("@<TRIPOS>ATOM"):
                    read_state = "atom"
                elif line.startswith("@<TRIPOS>BOND"):
                    read_state = "bond"

                if read_state == "mol":
                    mol_raw.append(line)
                elif read_state == "atom":
                    atom_raw.append(line)
                elif read_state == "bond":
                    bond_raw.append(line)

            # ats = Mol2Atoms.from_strlist(atom_raw)

            # print(str(ats))

            obj.atoms = Mol2Atoms.from_strlist(atom_raw)
            obj.header = "".join(mol_raw)
            obj.bonds = "".join(bond_raw)
            # print(mol_raw)
            # print(atom_raw)
            # print(bond_raw)

        return obj

    