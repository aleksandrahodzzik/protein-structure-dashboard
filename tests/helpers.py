"""Synthetic mmCIF builder shared by the structure and viewer tests."""

import itertools

ATOM_SITE_HEADER = """data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""

BACKBONE = ("N", "CA", "C")


def build_mmcif(
    records: list[tuple[str, str, str, tuple[str, ...]]],
    gaps_after: frozenset[int] = frozenset(),
    models: int = 1,
) -> str:
    """records: (group, residue name, chain, atom names). Residues are laid out along x,
    3.8 A apart, so consecutive backbones are 1.3 A C->N apart and PPBuilder links them.
    Residues after an index in ``gaps_after`` are shifted 10 A away, i.e. a chain break."""
    lines = [ATOM_SITE_HEADER.rstrip("\n")]
    serial = 1
    for model, (index, (group, resname, chain, atoms)) in itertools.product(
        range(1, models + 1), enumerate(records, start=1)
    ):
        for offset, atom in enumerate(atoms):
            shift = 10.0 * sum(gap < index for gap in gaps_after)
            x = 3.8 * index + shift + 1.2 * offset
            y = 0.0 if atom in BACKBONE else 2.0 + offset
            element = atom[0] if atom != "CA" or resname != "CA" else "CA"
            lines.append(
                f"{group} {serial} {element} {atom} . {resname} {chain} 1 {index} ? "
                f"{x:.3f} {y:.3f} 0.000 1.00 10.00 ? {index} {resname} {chain} {atom} {model}"
            )
            serial += 1
    return "\n".join(lines) + "\n#\n"


def residue(
    resname: str, chain: str = "A", group: str = "ATOM"
) -> tuple[str, str, str, tuple[str, ...]]:
    return (group, resname, chain, BACKBONE)
