# -*- coding: utf-8 -*-
"""
Cheminformatics feature extraction for QCElemental dimers (monomer A/B) to study
correlations with SAPT errors.

Notes / rationale for features:
- Donor/acceptor counts, formal charge, heteroatom content, rotors: capture HBond/ionic propensities and flexibility.
- RDKit Crippen molar refractivity (MR) is a rough proxy for molecular polarizability
  (empirically correlated; see e.g., Ghose & Crippen parameterization).
- Gasteiger-charges dipole (from 3D coords) provides an approximate gas-phase dipole;
  while crude, it’s useful as a directional electrostatics feature (Debye).
- Aromatic ring flags and ring-plane normals enable π–π vs T-shaped heuristics.
- Principal moments/asphericity provide coarse shape/polarizability anisotropy cues.

Caveats:
- Bond perception from 3D is heuristic; where possible, supply connectivity (e.g., from SDF).
- Gasteiger charges can fail on odd cases; we guard and fall back accordingly.

Author: (you)
"""
import argparse
import os
import tempfile
from typing import Dict, Tuple, Optional, List

import numpy as np
import pandas as pd
import qcelemental as qcel

from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    rdDetermineBonds,
    Descriptors,
    rdMolDescriptors,
    Crippen,
    Lipinski,
    rdMolTransforms,
)

import matplotlib.pyplot as plt
import seaborn as sns


DEBYE_PER_E_ANG = 4.80320427  # 1 e·Å = 4.8032 Debye


# ------------------------- RDKit conversion helpers ------------------------- #

def qcel_to_rdkit(mol_obj: qcel.models.Molecule) -> Optional[Chem.Mol]:
    """
    Convert a qcelemental Molecule -> RDKit Mol with 3D coords and perceived bonds.
    Uses a temporary XYZ (in Å) + RDKit bond perception (connectivity + bond orders).
    """
    # qcel stores geometry in Bohr by default; export explicitly in Å for XYZ
    tmp_path = 'temp.xyz'
    mol_obj.to_file(tmp_path)
    charge = int(mol_obj.molecular_charge)
    mol = Chem.MolFromXYZFile(tmp_path)
    # os.remove(tmp_path)
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
        rdDetermineBonds.DetermineConnectivity(mol, charge=charge)
        os.remove(tmp_path)
        # Sanitize & perceive aromaticity
        Chem.SanitizeMol(mol)
        # Keep 3D conformer from XYZ
        if mol.GetNumConformers() == 0:
            # Shouldn't happen for XYZ; but ensure a 3D conformer exists
            AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)

        return mol
    except Exception as e:
        # Bond perception or sanitization failed
        print(f"RDKit conversion failed: {e}")
        return None

# ------------------------------ Feature blocks ------------------------------ #

def compute_gasteiger_dipole(mol: Chem.Mol) -> Tuple[float, np.ndarray]:
    """
    Compute an approximate dipole from Gasteiger charges and 3D coordinates.
    Returns (magnitude in Debye, vector in Debye).
    """
    m = Chem.Mol(mol)
    try:
        Chem.rdPartialCharges.ComputeGasteigerCharges(m)
    except Exception:
        return 0.0, np.zeros(3)

    conf = m.GetConformer()
    dip_vec = np.zeros(3, dtype=float)
    for atom in m.GetAtoms():
        qi = atom.GetDoubleProp("_GasteigerCharge") if atom.HasProp("_GasteigerCharge") else 0.0
        pos = conf.GetAtomPosition(atom.GetIdx())
        r = np.array([pos.x, pos.y, pos.z], dtype=float)  # Å
        dip_vec += qi * r  # e * Å
    dip_vec *= DEBYE_PER_E_ANG
    return float(np.linalg.norm(dip_vec)), dip_vec


def ring_plane_normals(mol: Chem.Mol) -> List[np.ndarray]:
    """
    For each aromatic ring, compute a plane normal using best-fit plane (PCA) of ring atoms.
    Returns a list of unit normals (np.ndarray, shape (3,)).
    """
    conf = mol.GetConformer()
    ri = mol.GetRingInfo()
    normals = []
    for ring in ri.AtomRings():
        ring_atoms = list(ring)
        if not any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring_atoms):
            continue
        coords = []
        for idx in ring_atoms:
            pos = conf.GetAtomPosition(idx)
            coords.append([pos.x, pos.y, pos.z])
        X = np.array(coords, dtype=float)
        # subtract centroid
        Xc = X - X.mean(axis=0, keepdims=True)
        # PCA: smallest singular vector is plane normal
        try:
            _, _, vh = np.linalg.svd(Xc, full_matrices=False)
            n = vh[-1, :]
            n /= np.linalg.norm(n) + 1e-15
            normals.append(n)
        except np.linalg.LinAlgError:
            continue
    return normals


def monomer_features(mol: Chem.Mol) -> Dict[str, float]:
    """
    Compute monomer descriptors plausibly relevant to noncovalent interactions.
    """
    props: Dict[str, float] = {}
    # Basic composition / counts
    props["mol_weight"] = Descriptors.MolWt(mol)
    props["num_rings"] = rdMolDescriptors.CalcNumRings(mol)
    props["num_aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
    props["is_aromatic_any"] = 1 if props["num_aromatic_rings"] > 0 else 0
    props["num_rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
    props["formal_charge"] = Chem.GetFormalCharge(mol)
    props["num_hetero"] = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6))
    props["num_heavy"] = mol.GetNumHeavyAtoms()
    props["frac_csp3"] = rdMolDescriptors.CalcFractionCSP3(mol)
    props["hba"] = Lipinski.NumHAcceptors(mol)
    props["hbd"] = Lipinski.NumHDonors(mol)
    props["halogen_count"] = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (9, 17, 35, 53, 85))

    # Physchem: logP, TPSA, MR (polarizability proxy)
    props["logp"] = Crippen.MolLogP(mol)
    props["tpsa"] = rdMolDescriptors.CalcTPSA(mol)
    props["molar_refractivity"] = Crippen.MolMR(mol)

    # Approximate dipole from Gasteiger charges
    dip_mag, dip_vec = compute_gasteiger_dipole(mol)
    props["dipole_debye"] = dip_mag
    props["dipole_x"] = float(dip_vec[0])
    props["dipole_y"] = float(dip_vec[1])
    props["dipole_z"] = float(dip_vec[2])

    # Ring plane normals (store first two angles later in dimer features)
    normals = ring_plane_normals(mol)
    props["num_ring_normals"] = len(normals)

    return props


def centroid(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    coords = np.array([[p[0], p[1], p[2]] for p in conf.GetPositions()], dtype=float)
    return coords.mean(axis=0)


def min_heavy_atom_distance(molA: Chem.Mol, molB: Chem.Mol) -> float:
    confA, confB = molA.GetConformer(), molB.GetConformer()
    coordsA = np.array([[p[0], p[1], p[2]] for p in confA.GetPositions()], dtype=float)
    coordsB = np.array([[p[0], p[1], p[2]] for p in confB.GetPositions()], dtype=float)
    # heavy only
    heavyA = [i for i, a in enumerate(molA.GetAtoms()) if a.GetAtomicNum() > 1]
    heavyB = [i for i, a in enumerate(molB.GetAtoms()) if a.GetAtomicNum() > 1]
    if not heavyA or not heavyB:
        useA = range(molA.GetNumAtoms())
        useB = range(molB.GetNumAtoms())
    else:
        useA, useB = heavyA, heavyB
    A = coordsA[np.array(useA)]
    B = coordsB[np.array(useB)]
    # pairwise min distance
    dmin = np.min(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1))
    return float(dmin)


def angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
    n1 = v1 / (np.linalg.norm(v1) + 1e-15)
    n2 = v2 / (np.linalg.norm(v2) + 1e-15)
    cosang = np.clip(np.dot(n1, n2), -1.0, 1.0)
    return float(np.degrees(np.arccos(cosang)))


def best_ring_normal(mol: Chem.Mol) -> Optional[np.ndarray]:
    normals = ring_plane_normals(mol)
    if len(normals) == 0:
        return None
    # Heuristic: pick the normal from the largest ring system by atoms
    ri = mol.GetRingInfo()
    sizes = []
    for ring in ri.AtomRings():
        sizes.append(len(ring))
    # map normals in order they were added (only aromatic rings were included)
    # choose the longest aromatic ring; approximate mapping
    idx = int(np.argmax(sizes)) if sizes else 0
    idx = min(idx, len(normals) - 1)
    return normals[idx]


def dimer_features(molA: Chem.Mol, molB: Chem.Mol,
                   dipA_vec: np.ndarray, dipB_vec: np.ndarray) -> Dict[str, float]:
    feats: Dict[str, float] = {}
    comA, comB = centroid(molA), centroid(molB)
    vecAB = comB - comA
    feats["com_distance"] = float(np.linalg.norm(vecAB))
    feats["min_heavy_dist"] = min_heavy_atom_distance(molA, molB)

    # Dipole-dipole angle & alignment with intermolecular axis
    if np.linalg.norm(dipA_vec) > 1e-8 and np.linalg.norm(dipB_vec) > 1e-8:
        feats["dipole_angle_AB_deg"] = angle_between(dipA_vec, dipB_vec)
        feats["dipoleA_axis_angle_deg"] = angle_between(dipA_vec, vecAB)
        feats["dipoleB_axis_angle_deg"] = angle_between(dipB_vec, vecAB)
    else:
        feats["dipole_angle_AB_deg"] = np.nan
        feats["dipoleA_axis_angle_deg"] = np.nan
        feats["dipoleB_axis_angle_deg"] = np.nan

    # π–π / T-shape heuristics from ring plane normals (if any)
    nA = best_ring_normal(molA)
    nB = best_ring_normal(molB)
    if nA is not None and nB is not None:
        feats["ring_normal_angle_deg"] = angle_between(nA, nB)
        feats["axis_to_planeA_deg"] = angle_between(vecAB, nA)
        feats["axis_to_planeB_deg"] = angle_between(vecAB, nB)
        # crude classification flags (non-exclusive)
        # parallel planes and axis roughly perpendicular => stacked
        feats["pi_pi_like"] = int(
            (feats["ring_normal_angle_deg"] <= 30.0 or feats["ring_normal_angle_deg"] >= 150.0)
            and (80.0 <= feats["axis_to_planeA_deg"] <= 100.0)
            and (80.0 <= feats["axis_to_planeB_deg"] <= 100.0)
        )
        # one plane ~perpendicular to axis (T-like), the other ~parallel
        feats["t_shape_like"] = int(
            (70.0 <= feats["axis_to_planeA_deg"] <= 110.0 and feats["axis_to_planeB_deg"] <= 30.0)
            or (70.0 <= feats["axis_to_planeB_deg"] <= 110.0 and feats["axis_to_planeA_deg"] <= 30.0)
        )
    else:
        feats["ring_normal_angle_deg"] = np.nan
        feats["axis_to_planeA_deg"] = np.nan
        feats["axis_to_planeB_deg"] = np.nan
        feats["pi_pi_like"] = 0
        feats["t_shape_like"] = 0

    return feats


def smiles_safe(mol: Optional[Chem.Mol]) -> str:
    return Chem.MolToSmiles(mol) if mol is not None else "N/A"


# ----------------------------------- Main ----------------------------------- #

def main(input_file: str, output_file: str):
    df = pd.read_pickle(input_file)

    # Error (kcal/mol); sign convention matches user's original line
    df["error"] = (df["SAPT0 TOTAL ENERGY atqz"] - df["SAPT2+3(CCD) TOTAL ENERGY atqz"]) * 627.5

    rows = []
    for _, row in df.iterrows():
        # Build QCElemental dimer (geometry likely provided in bohr; user's original divides by bohr2angstroms
        # when creating qcel Molecule. Here we build using their original logic to preserve behavior.)
        dimer = qcel.models.Molecule.from_data({
            "symbols": row["atomic_numbers"],
            "geometry": row["coordinates"] / qcel.constants.bohr2angstroms,  # preserve user’s convention
            "molecular_charge": row["dimer_charge"],
            "molecular_multiplicity": row["dimer_multiplicity"],
            "fragments": [row["monAs"], row["monBs"]],
            "fragment_charges": [row["monA_charge"], row["monB_charge"]],
            "fragment_multiplicities": [row["monA_multiplicity"], row["monB_multiplicity"]],
        })

        # Extract monomers
        monA_qc = dimer.get_fragment(0)
        monB_qc = dimer.get_fragment(1)

        # Convert to RDKit
        molA = qcel_to_rdkit(monA_qc)
        molB = qcel_to_rdkit(monB_qc)

        # If conversion failed, skip entry gracefully
        if molA is None or molB is None:
            continue

        # SMILES (canonical; no stereochem from 3D assumed)
        smilesA = smiles_safe(molA)
        smilesB = smiles_safe(molB)
        print(f"{smilesA} -- {smilesB}")

        # Monomer features
        fa = monomer_features(molA)
        fb = monomer_features(molB)

        # Dimer features (use dipole vectors we already computed)
        dipA_vec = np.array([fa["dipole_x"], fa["dipole_y"], fa["dipole_z"]], dtype=float)
        dipB_vec = np.array([fb["dipole_x"], fb["dipole_y"], fb["dipole_z"]], dtype=float)
        fd = dimer_features(molA, molB, dipA_vec, dipB_vec)

        # Aggregate row
        out = {
            "error": float(row["error"]),
            "smiles_a": smilesA,
            "smiles_b": smilesB,
        }
        # Copy selected monomer features with suffixes
        for k, v in fa.items():
            out[f"{k}_a"] = v
        for k, v in fb.items():
            out[f"{k}_b"] = v
        # Cross features
        out.update(fd)

        # Handy composite flags
        out["contains_aromatic"] = int(out["is_aromatic_any_a"] or out["is_aromatic_any_b"])
        out["any_halogen"] = int((fa["halogen_count"] > 0) or (fb["halogen_count"] > 0))
        out["any_donor_acceptor"] = int((fa["hbd"] > 0 and fb["hba"] > 0) or (fb["hbd"] > 0 and fa["hba"] > 0))

        rows.append(out)

    results_df = pd.DataFrame(rows)
    if results_df.empty:
        print("No RDKit conversions succeeded; nothing to plot.")
        return

    # Preview
    print(results_df.head(10))

    # Visualization: a small, informative subset
    pair_plot_vars = [
        "error",
        "molar_refractivity_a", "molar_refractivity_b",
        "dipole_debye_a", "dipole_debye_b",
        "logp_a", "logp_b",
        "tpsa_a", "tpsa_b",
        "com_distance", "min_heavy_dist",
    ]
    # Keep only columns that exist (robustness)
    pair_plot_vars = [c for c in pair_plot_vars if c in results_df.columns]

    # Drop NaNs for plotting
    plot_df = results_df[pair_plot_vars + ["contains_aromatic"]].dropna()

    if len(plot_df) >= 2 and len(pair_plot_vars) >= 2:
        sns.pairplot(plot_df, vars=pair_plot_vars, hue="contains_aromatic", diag_kind="kde", corner=True)
        plt.tight_layout()
        plt.savefig(output_file, dpi=200, bbox_inches="tight")
        print(f"Plot saved to {output_file}")
    else:
        print("Insufficient data for pairplot; skipping plot.")

    # Optionally write the enriched dataframe for downstream modeling
    enriched_path = os.path.splitext(output_file)[0] + "_features.parquet"
    results_df.to_parquet(enriched_path, index=False)
    print(f"Feature table written to {enriched_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze SAPT errors vs. monomer/dimer cheminformatics features.")
    parser.add_argument("--input_file", type=str, default="combined_df_4569.pkl", help="Pickle path with dimer & SAPT columns.")
    parser.add_argument("--output_file", type=str, default="cheminformatics_analysis.png", help="Where to save the pairplot.")
    args = parser.parse_args()
    main(args.input_file, args.output_file)
