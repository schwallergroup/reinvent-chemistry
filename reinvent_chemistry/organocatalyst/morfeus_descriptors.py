import os
import shutil

import numpy as np

from rdkit.Chem import Mol
from reinvent_chemistry.organocatalyst.geometry_optimizer import GeometryOptimizer

from morfeus import read_xyz, XTB


class MorfeusDescriptors:
    """Quantum mechanical descriptors from semi-empirical xTB. Wraps MORFEUS functionalities"""
    def __init__(self):
        self._geometry_optimizer = GeometryOptimizer()

    def ionization_potential(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        elements, coordinates = read_xyz(geometry_path)
        xtb = XTB(elements, coordinates)
        self._clean_up_temp_dir(path=temp_dir)

        return xtb.get_ip(corrected=True)

    def electron_affinity(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        elements, coordinates = read_xyz(geometry_path)
        xtb = XTB(elements, coordinates)
        self._clean_up_temp_dir(path=temp_dir)

        return xtb.get_ea(corrected=True)

    def homo(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        elements, coordinates = read_xyz(geometry_path)
        xtb = XTB(elements, coordinates)
        self._clean_up_temp_dir(path=temp_dir)

        return xtb.get_homo()

    def lumo(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        elements, coordinates = read_xyz(geometry_path)
        xtb = XTB(elements, coordinates)
        self._clean_up_temp_dir(path=temp_dir)

        return xtb.get_lumo()

    def dipole(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        elements, coordinates = read_xyz(geometry_path)
        xtb = XTB(elements, coordinates)
        dipoles = xtb.get_dipole()
        self._clean_up_temp_dir(path=temp_dir)

        return np.sqrt(np.array(dipoles).dot(np.array(dipoles).T))

    def global_electrophilicity(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        elements, coordinates = read_xyz(geometry_path)
        xtb = XTB(elements, coordinates)
        self._clean_up_temp_dir(path=temp_dir)

        return xtb.get_global_descriptor("electrophilicity", corrected=True)

    def global_nucleophilicity(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        elements, coordinates = read_xyz(geometry_path)
        xtb = XTB(elements, coordinates)
        self._clean_up_temp_dir(path=temp_dir)

        return xtb.get_global_descriptor("nucleophilicity", corrected=True)

    def _get_optimized_geometry_path(self, mol: Mol):
        temp_dir = self._geometry_optimizer.optimize_xtb_geometry(mol=mol)
        return temp_dir, os.path.join(temp_dir, "temp.xtbopt.xyz")

    @staticmethod
    def _clean_up_temp_dir(path: str):
        shutil.rmtree(path)

