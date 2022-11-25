import unittest
from reinvent_chemistry.conversions import Conversions
from reinvent_chemistry.organocatalyst.morfeus_descriptors import MorfeusDescriptors


class TestMorfeusDescriptors(unittest.TestCase):
    def setUp(self) -> None:
        self.ammonia_smiles = "N"
        self.water_smiles = "O"
        self._conversions = Conversions()
        self.morfeus_descriptors = MorfeusDescriptors()

    def test_ionization_potential(self):
        self.assertAlmostEqual(self.morfeus_descriptors.ionization_potential(self._conversions.smile_to_mol(self.ammonia_smiles)), 11.5859, 3)
        self.assertAlmostEqual(self.morfeus_descriptors.ionization_potential(self._conversions.smile_to_mol(self.water_smiles)), 13.4469, 3)

    def test_electron_affinity(self):
        self.assertAlmostEqual(self.morfeus_descriptors.electron_affinity(self._conversions.smile_to_mol(self.ammonia_smiles)), -10.6308, 3)
        self.assertAlmostEqual(self.morfeus_descriptors.electron_affinity(self._conversions.smile_to_mol(self.water_smiles)), -11.5145, 3)

    def test_homo(self):
        self.assertAlmostEqual(self.morfeus_descriptors.homo(self._conversions.smile_to_mol(self.ammonia_smiles)), -0.3859, 3)
        self.assertAlmostEqual(self.morfeus_descriptors.homo(self._conversions.smile_to_mol(self.water_smiles)), -0.4464, 3)

    def test_lumo(self):
        self.assertAlmostEqual(self.morfeus_descriptors.lumo(self._conversions.smile_to_mol(self.ammonia_smiles)), 0.06358, 3)
        self.assertAlmostEqual(self.morfeus_descriptors.lumo(self._conversions.smile_to_mol(self.water_smiles)), 0.08233, 3)

    def test_dipole(self):
        self.assertAlmostEqual(self.morfeus_descriptors.dipole(self._conversions.smile_to_mol(self.ammonia_smiles)), 0.7009, 3)
        self.assertAlmostEqual(self.morfeus_descriptors.dipole(self._conversions.smile_to_mol(self.water_smiles)), 0.8717, 3)

    def test_global_electrophilicity(self):
        self.assertAlmostEqual(self.morfeus_descriptors.global_electrophilicity(self._conversions.smile_to_mol(self.ammonia_smiles)), 0.00513, 3)
        self.assertAlmostEqual(self.morfeus_descriptors.global_electrophilicity(self._conversions.smile_to_mol(self.water_smiles)), 0.01870, 3)

    def test_global_nucleophilicity(self):
        self.assertAlmostEqual(self.morfeus_descriptors.global_nucleophilicity(self._conversions.smile_to_mol(self.ammonia_smiles)), -11.5859, 3)
        self.assertAlmostEqual(self.morfeus_descriptors.global_nucleophilicity(self._conversions.smile_to_mol(self.water_smiles)), -13.4469, 3)
