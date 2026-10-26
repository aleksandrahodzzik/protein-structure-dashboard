"""The example gallery: 128 hand-picked entries in 16 themed groups of 8.

Pure data, so it is easy to extend and to test. Every ID was fetched from RCSB and parsed
into at least one protein chain before it was added; a few of them (e.g. 1AON, 1RYP, 3WU2)
are above the heavy-structure threshold on purpose, so the backbone-only mode has
something to show.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Example:
    pdb_id: str
    name: str
    blurb: str


@dataclass(frozen=True)
class Category:
    name: str
    items: tuple[Example, ...]


CATEGORIES: tuple[Category, ...] = (
    Category(
        "Classics",
        (
            Example("1UBQ", "Ubiquitin", "Small, 76 residues"),
            Example("1MBN", "Myoglobin", "Sperm whale, an early 3D protein model"),
            Example("4HHB", "Hemoglobin", "Multi-chain, heme ligands"),
            Example("1CRN", "Crambin", "Only 46 residues, from a plant seed"),
            Example("1LYZ", "Lysozyme", "Hen egg white, a textbook enzyme"),
            Example("5PTI", "BPTI", "Bovine trypsin inhibitor, 58 residues"),
            Example("2LZM", "T4 lysozyme", "Phage enzyme used in stability studies"),
            Example("4INS", "Insulin", "Hormone made of A and B chains"),
        ),
    ),
    Category(
        "Membrane",
        (
            Example("1BL8", "KcsA K+ channel", "Hydrophobic transmembrane region"),
            Example("1F88", "Rhodopsin", "Vision receptor, seven helices"),
            Example("2RH1", "Beta-2 adrenergic receptor", "Human GPCR with a T4 lysozyme fusion"),
            Example("1J4N", "Aquaporin-1", "A channel that passes water"),
            Example("1C3W", "Bacteriorhodopsin", "Light-driven proton pump"),
            Example("2POR", "Porin", "Beta-barrel channel in a bacterial membrane"),
            Example("3SN6", "Beta-2 receptor + Gs", "Receptor bound to its G protein"),
            Example("4DJH", "Kappa opioid receptor", "Human GPCR bound to JDTic"),
        ),
    ),
    Category(
        "Enzymes",
        (
            Example("1TIM", "Triosephosphate isomerase", "Chicken muscle, TIM-barrel fold"),
            Example("1AKE", "Adenylate kinase", "E. coli, bound to an inhibitor"),
            Example("3PTB", "Trypsin", "Bovine digestive protease"),
            Example("7RSA", "Ribonuclease A", "Bovine, a classic folding model"),
            Example("1RX2", "DHFR", "E. coli, with folate and NADP"),
            Example("1CA2", "Carbonic anhydrase II", "Human zinc enzyme"),
            Example("1ACJ", "Acetylcholinesterase", "Torpedo ray, deep active-site gorge"),
            Example("4CHA", "Chymotrypsin", "Bovine protease, refined at 1.68 Å"),
        ),
    ),
    Category(
        "Signalling",
        (
            Example("1EMA", "Green fluorescent protein", "A beta-barrel that glows green"),
            Example("2H5Q", "mCherry", "A red fluorescent protein"),
            Example("1CLL", "Calmodulin", "Calcium sensor with a dumbbell shape"),
            Example("5P21", "H-Ras", "Small GTPase, a frequent cancer driver"),
            Example("1ATP", "Protein kinase A", "Catalytic subunit with a peptide"),
            Example("1M17", "EGFR kinase", "Kinase domain with an inhibitor"),
            Example("1STP", "Streptavidin", "Binds biotin with extreme affinity"),
            Example("1HCK", "CDK2", "A cell-cycle kinase"),
        ),
    ),
    Category(
        "Immune",
        (
            Example("1IGT", "Mouse IgG antibody", "Intact immunoglobulin"),
            Example("1HZH", "Human IgG b12", "Antibody active against HIV"),
            Example("1AO6", "Serum albumin", "Human blood carrier protein"),
            Example("1HSA", "HLA-B27", "Antigen-presenting MHC protein"),
            Example("4ZQK", "PD-1 / PD-L1", "An immune checkpoint pair"),
            Example("1BJ1", "VEGF + antibody", "Growth factor bound by a Fab"),
            Example("1MLC", "Fab D44.1", "Antibody bound to lysozyme"),
            Example("3HFM", "HyHEL-10 Fab", "Another antibody-lysozyme complex"),
        ),
    ),
    Category(
        "DNA & RNA",
        (
            Example("1TUP", "p53 on DNA", "Tumour suppressor, three chains"),
            Example("1AOI", "Nucleosome", "Histones that package DNA"),
            Example("4OO8", "CRISPR-Cas9", "Genome editor, 1,300 residues"),
            Example("1CGP", "CAP on DNA", "Bacterial gene activator"),
            Example("1AAY", "Zif268", "Zinc fingers gripping DNA"),
            Example("1TAQ", "Taq polymerase", "The enzyme behind PCR"),
            Example("1CEZ", "T7 RNA polymerase", "Transcribing from a T7 promoter"),
            Example("1LMB", "Lambda repressor", "Phage switch bound to DNA"),
        ),
    ),
    Category(
        "Viruses",
        (
            Example("6VXX", "SARS-CoV-2 spike", "Large, three chains"),
            Example("6LU7", "SARS-CoV-2 protease", "Main protease with an inhibitor"),
            Example("6M0J", "Spike RBD + ACE2", "How the virus grips human cells"),
            Example("7BV2", "SARS-CoV-2 polymerase", "Cryo-EM, with remdesivir"),
            Example("1HVR", "HIV protease", "A drug target, two chains"),
            Example("1RU7", "Hemagglutinin", "Influenza surface protein"),
            Example("2HU4", "Neuraminidase", "Influenza enzyme with oseltamivir"),
            Example("1RTD", "HIV reverse transcriptase", "A catalytic complex"),
        ),
    ),
    Category(
        "Machines",
        (
            Example("1BMF", "F1-ATPase", "Bovine ATP-making rotor"),
            Example("1AON", "GroEL-GroES", "Chaperone cage, backbone-only view"),
            Example("1RYP", "20S proteasome", "Yeast protein shredder, backbone-only"),
            Example("1BG2", "Kinesin", "Motor domain that walks on microtubules"),
            Example("1ATN", "Actin + DNase I", "Muscle filament subunit"),
            Example("2FHA", "Ferritin", "Human iron-storage protein"),
            Example("1HRC", "Cytochrome c", "Horse heart electron carrier"),
            Example("1OCC", "Cytochrome c oxidase", "Bovine respiratory complex, 26 chains"),
        ),
    ),
    Category(
        "Plants & light",
        (
            Example("1AUS", "Rubisco", "Spinach enzyme that fixes CO2"),
            Example("1PRC", "Photosynthetic reaction centre", "Purple bacterium, 4 chains"),
            Example("1JB0", "Photosystem I", "Cyanobacterial complex, 12 chains"),
            Example("1KZU", "LH2 light harvester", "A ring-shaped antenna"),
            Example("1FNB", "Ferredoxin reductase", "Spinach, passes electrons to NADP+"),
            Example("3WU2", "Photosystem II", "Splits water, backbone-only view"),
            Example("3CNA", "Concanavalin A", "Jack bean lectin"),
            Example("1PLC", "Plastocyanin", "Poplar copper electron carrier"),
        ),
    ),
    Category(
        "Metalloproteins",
        (
            Example("1AZU", "Azurin", "A blue copper protein"),
            Example("2CCP", "Cytochrome c peroxidase", "Yeast heme enzyme"),
            Example("1DGF", "Catalase", "Human, breaks down hydrogen peroxide"),
            Example("2SOD", "Superoxide dismutase", "Bovine Cu,Zn enzyme"),
            Example("1RB9", "Rubredoxin", "Iron-sulfur protein at 0.92 Å"),
            Example("1N2C", "Nitrogenase complex", "Fixes atmospheric nitrogen"),
            Example("1YCC", "Yeast cytochrome c", "Iso-1, 108 residues"),
            Example("1CYO", "Cytochrome b5", "Bovine heme-binding protein"),
        ),
    ),
    Category(
        "Brain & disease",
        (
            Example("2BEG", "Amyloid-beta fibril", "Alzheimer's Abeta(1-42), NMR"),
            Example("1QLX", "Human prion protein", "Folded domain, NMR"),
            Example("2N0A", "Alpha-synuclein fibril", "Linked to Parkinson's disease"),
            Example("5O3L", "Tau filament", "Paired helical filament, Alzheimer's brain"),
            Example("1FKN", "BACE1", "Beta-secretase with an inhibitor"),
            Example("1I9B", "AChBP", "Snail acetylcholine-binding protein"),
            Example("3KG2", "AMPA receptor", "Rat glutamate receptor"),
            Example("4M48", "Dopamine transporter", "Fly transporter with an antidepressant"),
        ),
    ),
    Category(
        "Cancer targets",
        (
            Example("1YCR", "MDM2 + p53 peptide", "A target for restoring p53"),
            Example("1G5M", "Bcl-2", "Protein that blocks cell death"),
            Example("1IEP", "Abl kinase + imatinib", "The Gleevec target"),
            Example("1UWH", "B-Raf + sorafenib", "Kinase with a bound drug"),
            Example("1ERE", "Estrogen receptor", "Ligand-binding domain with estradiol"),
            Example("1YET", "HSP90", "Bound to geldanamycin"),
            Example("4OBE", "KRAS", "GDP-bound human KRas"),
            Example("1N8Z", "HER2 + Herceptin", "Antibody drug bound to HER2"),
        ),
    ),
    Category(
        "Hormones & cytokines",
        (
            Example("3HHR", "Growth hormone + receptor", "Hormone gripped by its receptor"),
            Example("1EER", "Erythropoietin + receptor", "Red blood cell signal"),
            Example("1M47", "Interleukin-2", "Human T-cell cytokine"),
            Example("1TNF", "TNF-alpha", "Inflammation signal, a trimer"),
            Example("1BET", "Nerve growth factor", "Mouse beta-NGF"),
            Example("1AX8", "Leptin", "Human obesity protein"),
            Example("4ZXB", "Insulin receptor", "Ectodomain in a complex"),
            Example("1ITF", "Interferon alpha-2a", "Antiviral cytokine, NMR"),
        ),
    ),
    Category(
        "Cytoskeleton",
        (
            Example("1JFF", "Tubulin", "Alpha-beta dimer with taxol"),
            Example("2MYS", "Myosin S1", "The muscle motor head"),
            Example("1TIT", "Titin Ig domain", "Part of the giant muscle spring"),
            Example("1CAG", "Collagen-like peptide", "A triple helix"),
            Example("1SA0", "Tubulin + colchicine", "A drug bound to tubulin"),
            Example("1C1G", "Tropomyosin", "Coiled coil along actin filaments"),
            Example("1AJ3", "Spectrin repeat", "A three-helix bundle, NMR"),
            Example("3VKG", "Dynein motor domain", "Huge, about 2,900 residues"),
        ),
    ),
    Category(
        "Folding & design",
        (
            Example("1L2Y", "Trp-cage", "A 20-residue mini-protein"),
            Example("1PGA", "Protein G B1 domain", "A folding model, 56 residues"),
            Example("1BDD", "Protein A domain", "Three-helix bundle, 60 residues"),
            Example("1A2P", "Barnase", "Bacterial ribonuclease"),
            Example("2CI2", "Chymotrypsin inhibitor 2", "Barley, a classic folding model"),
            Example("1QYS", "Top7", "A computationally designed protein"),
            Example("1SHG", "SH3 domain", "A small signalling module"),
            Example("1VII", "Villin headpiece", "36 residues, a fast folder"),
        ),
    ),
    Category(
        "Toxins & venoms",
        (
            Example("7AHL", "Alpha-hemolysin", "A pore built from seven chains"),
            Example("3BTA", "Botulinum neurotoxin A", "Neurotoxin, 1,277 residues"),
            Example("1XTC", "Cholera toxin", "Multi-chain bacterial toxin"),
            Example("1DDT", "Diphtheria toxin", "Shuts down protein synthesis"),
            Example("1ACC", "Anthrax protective antigen", "The toxin's delivery protein"),
            Example("1HC9", "alpha-Bungarotoxin", "Snake toxin with a peptide"),
            Example("2AAI", "Ricin", "Plant toxin, two chains"),
            Example("1OMN", "Conotoxin MVIIC", "Cone snail peptide, 26 residues"),
        ),
    ),
)

ALL: tuple[Example, ...] = tuple(item for category in CATEGORIES for item in category.items)
