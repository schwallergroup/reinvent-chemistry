import numpy as np
from rdkit.Chem import DataStructs
import rdkit.Chem as Chem
from espsim import EmbedAlignScore


class Similarity:

    def calculate_tanimoto(self, query_fps, ref_fingerprints) -> np.array:
        return np.array([np.max(DataStructs.BulkTanimotoSimilarity(fp, ref_fingerprints)) for fp in query_fps])

    def calculate_jaccard_distance(self, query_fps, ref_fingerprints) -> np.array:
        tanimoto = self.calculate_tanimoto(query_fps, ref_fingerprints)
        jaccard = 1 - tanimoto
        return jaccard

    
    def smiles_to_mol(self, smiles_strings):
        if not isinstance(smiles_strings, list) or not all(isinstance(s, str) for s in smiles_strings):
            raise TypeError("smiles_strings must be a list of strings")

        mols = [Chem.rdmolops.AddHs(Chem.MolFromSmiles(smi)) for smi in smiles_strings]
        return mols

    def calc_espsim(self, ref_smis, query_smis):
        if not isinstance(ref_smis, list) or not all(isinstance(s, str) for s in ref_smis) or not isinstance(query_smis, list) or not all(isinstance(s, str) for s in query_smis):
            raise TypeError("ref_smis and query_smis must be lists of strings")

        similarity_scores = []
        ref_mols = self.smiles_to_mol(ref_smis)
        query_mols = self.smiles_to_mol(query_smis)

        for ref in ref_mols:
            for query in query_mols: 
                simShape, simEsp = EmbedAlignScore(ref, query)
                similarity_scores.append((simShape[0], simEsp[0]))

        for i, scores in enumerate(similarity_scores):
            print(f"Scores for molecule pair {i+1}:")
            print(f"simShape: {scores[0]}, simEsp: {scores[1]}")

        return similarity_scores