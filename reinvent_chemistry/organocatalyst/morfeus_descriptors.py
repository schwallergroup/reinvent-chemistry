import os
import shutil

import numpy as np
import math

from rdkit.Chem import Mol
from reinvent_chemistry.organocatalyst.geometry_optimizer import GeometryOptimizer

from morfeus import read_xyz, XTB


class MorfeusDescriptors:
    """Quantum mechanical descriptors from semi-empirical xTB. Wraps MORFEUS functionalities"""
    def __init__(self):
        self._geometry_optimizer = GeometryOptimizer()

    def ionization_potential(self, mol: Mol) -> float:
        try:
            temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)

            return xtb.get_ip(corrected=True)

        except Exception:
            return 0.0

    def electron_affinity(self, mol: Mol) -> float:
        try:
            temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)

            return xtb.get_ea(corrected=True)

        except Exception:
            return 0.0

    def homo(self, mol: Mol) -> float:
        try:
            temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)

            return xtb.get_homo()

        except Exception:
            return 0.0

    def lumo(self, mol: Mol) -> float:
        try:
            temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)

            return xtb.get_lumo()

        except Exception:
            return 0.0

    def dipole(self, mol: Mol) -> float:
        try:
            temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            return xtb.get_ip(corrected=True)

        except Exception:
            return 0.0

    def electron_affinity(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            return xtb.get_ea(corrected=True)

        except Exception:
            return 0.0

    def homo(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            return xtb.get_homo()

        except Exception:
            return 0.0

    def lumo(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            return xtb.get_lumo()

        except Exception:
            return 0.0

    def dipole(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            dipoles = xtb.get_dipole()
            self._clean_up_temp_dir(path=temp_dir)

            return np.sqrt(np.array(dipoles).dot(np.array(dipoles).T))

        except Exception:
            return 0.0

    def global_electrophilicity(self, mol: Mol) -> float:
        try:
            temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)

            return xtb.get_global_descriptor("electrophilicity", corrected=True)

        except Exception:
            return 0.0

    def global_nucleophilicity(self, mol: Mol) -> float:
        try:
            temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)

            return xtb.get_global_descriptor("nucleophilicity", corrected=True)

        except Exception:
            return 0.0

    def global_electrophilicity(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            return xtb.get_global_descriptor("electrophilicity", corrected=True)

        except Exception:
            return 0.0

    def global_nucleophilicity(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            return xtb.get_global_descriptor("nucleophilicity", corrected=True)

        except Exception:
            return 0.0

    def chemical_potential(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            lumo_energy = xtb.get_lumo()
            homo_energy = xtb.get_homo()
            return (lumo_energy + homo_energy) / 2

        except Exception:
            return 0.0

    def hardness(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            lumo_energy = xtb.get_lumo()
            homo_energy = xtb.get_homo()
            return (lumo_energy - homo_energy) / 2

        except Exception:
            return 0.0

    def electrophilicity_index(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            lumo_energy = xtb.get_lumo()
            homo_energy = xtb.get_homo()
            chemical_potential = (lumo_energy + homo_energy) / 2
            hardness = (lumo_energy - homo_energy) / 2
            return (math.pow(chemical_potential, 2)) / (2 * hardness)

        except Exception:
            return 0.0

    def nucleophilicity_index(self, mol: Mol) -> float:
        temp_dir, geometry_path = self._get_optimized_geometry_path(mol=mol)
        try:
            elements, coordinates = read_xyz(geometry_path)
            xtb = XTB(elements, coordinates)
            self._clean_up_temp_dir(path=temp_dir)
            lumo_energy = xtb.get_lumo()
            homo_energy = xtb.get_homo()
            chemical_potential = (lumo_energy + homo_energy) / 2
            hardness = (lumo_energy - homo_energy) / 2
            E_index = (math.pow(chemical_potential, 2)) / (2 * hardness)
            return 1 / E_index

        except Exception:
            return 0.0

    def _get_optimized_geometry_path(self, mol: Mol):
        temp_dir = self._geometry_optimizer.optimize_xtb_geometry(mol=mol)
        return temp_dir, os.path.join(temp_dir, "temp.xtbopt.xyz")

    @staticmethod
    def _clean_up_temp_dir(path: str):
        shutil.rmtree(path)

