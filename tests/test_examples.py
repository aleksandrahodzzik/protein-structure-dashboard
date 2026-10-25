"""The gallery is data, so its invariants are cheap to guard."""

from collections import Counter

from protein_dashboard import examples
from protein_dashboard.rcsb_client import is_valid_pdb_id


def test_there_are_64_examples_in_8_groups_of_8() -> None:
    assert len(examples.ALL) == 64
    assert [len(c.items) for c in examples.CATEGORIES] == [8] * 8


def test_ids_are_valid_and_unique() -> None:
    ids = [e.pdb_id for e in examples.ALL]

    assert all(is_valid_pdb_id(pdb_id) and pdb_id == pdb_id.upper() for pdb_id in ids)
    assert [pdb_id for pdb_id, n in Counter(ids).items() if n > 1] == []


def test_group_names_are_unique_and_every_card_has_text_that_fits() -> None:
    names = [c.name for c in examples.CATEGORIES]

    assert len(names) == len(set(names))
    for example in examples.ALL:
        assert 0 < len(example.name) <= 28, example
        assert 0 < len(example.blurb) <= 45, example
