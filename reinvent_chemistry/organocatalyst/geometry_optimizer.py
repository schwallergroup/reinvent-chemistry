import os
import tempfile
import subprocess

from rdkit.Chem import Mol
from reinvent_chemistry import Conversions


class GeometryOptimizer:
    """optimize SMILES geometry"""
    def __init__(self):
        self._conversions = Conversions()

    def optimize_xtb_geometry(self, mol: Mol):
        """optimize geometry using OpenBabel then xTB (normal accuracy)"""
        smiles = self._conversions.mol_to_smiles(mol)
        # make temp folder to generate and store optimized geometries
        temp_dir = tempfile.mkdtemp()
        # generate un-optimized geometry with OpenBabel
        openbabel_path = os.path.join(temp_dir, "openbabel.xyz")
        subprocess.call(["obabel", f"-:{smiles}", "--gen3d", "-O", openbabel_path])
        # call xTB to optimize geometry
        subprocess.call(["xtb", openbabel_path, "--opt", "--namespace", f"{temp_dir}/temp"])

        return temp_dir



