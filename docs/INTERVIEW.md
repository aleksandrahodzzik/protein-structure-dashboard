# Interview notes: Protein Structure Dashboard

Written for my own preparation. Everything below describes what is actually in this repository; where something is a limitation, I say so.

## 60-second pitch

I built a Streamlit dashboard for exploring protein structures. You type a PDB ID or a keyword, it fetches the entry from the RCSB Protein Data Bank over their REST API, shows an interactive 3D model with Py3Dmol, and analyses each protein chain with Biopython: amino acid composition, molecular weight, isoelectric point, GRAVY, and a sliding-window Kyte-Doolittle hydrophobicity profile, all as interactive Plotly charts.

The part I care most about is the structure of the code. The analysis module is pure Python: it takes a string and returns numbers, with no Streamlit or network inside, so I could test it against known values like the molecular weight of ubiquitin. The RCSB client has timeouts, retries with backoff for 429 and 5xx, and is tested with mocked HTTP, so the whole test suite (132 tests, about 97% coverage) runs in CI without a network. It is containerised with a multi-stage, non-root Docker image, and GitHub Actions runs ruff, mypy, pytest and a Docker build on every push.

A lesson from using real data: the sequence of a PDB chain is not what you think it is. Only residues that were actually modelled have coordinates, so I analyse the observed sequence and say so in the UI instead of silently showing numbers for a different sequence than the one in the 3D view.

## Questions and answers

### 1. Why mmCIF and not the PDB format?
The legacy PDB format has fixed-width columns that cap out at 99,999 atoms, 62 chains and 3-character residue names. Large structures such as ribosomes and viral capsids cannot be written in it, and the PDB archive no longer guarantees PDB-format files for new entries. mmCIF is the archive's primary format; it is a self-describing key-value/table format, extensible, and Biopython's `MMCIFParser` reads it. The cost is that it is verbose and slower to parse, which is why I trim files before they go to the browser.

### 2. What is the difference between the observed sequence and SEQRES?
SEQRES is the full sequence of the molecule that was deposited. The observed sequence is the subset of residues that have atomic coordinates. Disordered loops, flexible termini and expression tags are often missing from the electron density, so the observed sequence can be shorter and have gaps. In 6VXX, chain A of the SARS-CoV-2 spike has 972 modelled residues in 12 segments, against 1,273 in the full chain. I analyse the observed sequence so the metrics, the plots and the 3D view all describe the same residues, and I show a caption explaining it. The trade-off is that the metrics (molecular weight especially) are not those of the full protein, which I list as a limitation.

### 3. What is GRAVY, and how does it differ from the hydrophobicity profile?
GRAVY (grand average of hydropathy) is the mean Kyte-Doolittle value over all residues of the chain: one number for the whole protein. Positive means more hydrophobic overall, negative means more hydrophilic (ubiquitin is about -0.49). The profile is the same scale averaged over a sliding window, so it keeps position information and shows where the hydrophobic stretches are. The single number can hide that a protein is, for instance, mostly polar with one transmembrane helix.

### 4. How does the isoelectric point work?
The pI is the pH at which the protein has zero net charge. Each ionisable group (N-terminus, C-terminus, and the side chains of Asp, Glu, His, Cys, Tyr, Lys, Arg) contributes a charge that depends on its pKa and the pH, via Henderson-Hasselbalch. Biopython sums these and searches for the pH where the total is zero. It is an estimate: pKa values shift with the local environment in a folded protein, which this simple model ignores, and different tools use different pKa tables, so values differ slightly between them.

### 5. Why a window of 9 by default, and 19 for the transmembrane preset?
These come from the original Kyte and Doolittle paper (1982). A short window (about 9) gives a profile that highlights local hydrophobic and hydrophilic stretches, such as buried versus exposed segments. A membrane-spanning alpha helix is around 20 residues long, so a window of 19 averages across roughly one helix, and a mean above about 1.6 is the paper's suggested cut-off for a likely transmembrane segment. The window must be odd so that it has a centre residue, and each point is plotted at that centre. For a chain of length L and window w the profile has L - w + 1 values, which is one of my tests. The UI caps the window at the chain length and at 39. For KcsA (1BL8) the preset shows peaks above the threshold around the helices, which I checked in the app; it is a heuristic, not a prediction tool.

### 6. Average versus monoisotopic mass: which one do you report?
I report the average mass, which weights each element by its natural isotope abundances; it is what you get from the molecular formula and is the usual "molecular weight" for a protein on a gel or in a table. Monoisotopic mass uses only the most abundant isotope of each element, and is what high-resolution mass spectrometry resolves for small peptides. For large proteins the monoisotopic peak is tiny and not practically observed, so average is the right default. The analysis function has a `monoisotopic` flag, but the UI does not expose it. Ubiquitin comes out at 8,564.7 Da, and a test checks that to within 1 Da.

### 7. Why did you structure the code as pure functions, and how does that help testing?
`analysis.py` has no I/O and no Streamlit imports, so the tests are just input and output: known sequence in, known number out, and no mocks, fixtures or browser. The same applies to `charts.py` (data in, Plotly figure out) and most of `structure.py` (text in, dataclasses out). The side effects are pushed to the edges: HTTP in `rcsb_client.py` and the UI in `app.py`. For the UI, the one place where I construct a client, `get_client()`, is the seam where the test swaps in a fake, and Streamlit's `AppTest` then runs the real app code against it.

### 8. How do you test an HTTP client without the network?
With the `responses` library, which patches the `requests` transport layer. The tests cover success, 404 (mapped to a specific "not found" error and not retried), 429 followed by success (checking it retries and that `Retry-After` is honoured), timeouts, and exhaustion of the retries. The sleep function is injectable, so a test that simulates three retries takes milliseconds instead of seconds. I also test ID validation, so bad input never reaches the network.

### 9. How do you handle caching and RCSB rate limits?
Metadata, file text, parsed structures and search results are cached with `st.cache_data` with a 24-hour TTL. PDB entries are practically immutable once released, so the cache is safe, and it means that moving a slider or switching a chart does not re-fetch anything. Streamlit reruns the script on every interaction, so without caching the app would hammer RCSB. The client sends an identifying `User-Agent`, uses timeouts, and backs off exponentially on 429/5xx. The cache is per-process and in memory, and is cleared on restart, which is fine here but is not a solution for many users (see question 13).

### 10. Why `components.html` instead of `stmol`?
`stmol` is a third-party wrapper that adds a dependency and a compatibility risk with Streamlit versions, and Py3Dmol already generates a complete HTML/JS page with `_make_html()`. Embedding that with Streamlit's own `components.html` needs no additional package. I put it inside a single function, `viewer.show_view`. The honest caveat is that newer Streamlit versions mark `components.html` as deprecated in favour of `st.iframe`, so there will be a migration, and it is a one-function change because of that wrapper.

### 11. What did you do about very large structures?
The browser is the bottleneck. A ribosome entry (4V6X) has about 237,000 atoms and took about 22 seconds to load. I first trim the mmCIF on the server: first model only (NMR entries contain many models that would otherwise be stacked on top of each other), and only the atoms needed for a cartoon. Above 50,000 atoms the viewer goes further: backbone only, with no surface, no ligands and no chain dimming, and the UI shows a warning. I picked the threshold from testing a few entries, not from a formal benchmark, so it is a judgement call.

### 12. How do you deal with modified residues like selenomethionine?
Structures contain non-standard residues (for example MSE, selenomethionine, common in X-ray phasing) which are written as HETATM records. A naive parser treats them as ligands and drops them, leaving a hole in the sequence. I treat a residue as part of the polymer if it is a standard amino acid by name or has a full N/CA/C backbone, and I map it to its parent amino acid using Biopython's extended table (MSE becomes M, SEP becomes S, and so on). Anything the table does not know, including `UNK`, becomes `X` rather than a guess, so positions in the sequence stay aligned with the structure; `X` residues are excluded from composition and mass and the UI reports how many there were. The limitation is that a modified residue is analysed as its parent, so for example the mass of selenomethionine is that of methionine, which is slightly off.

### 13. How would you scale this to many users or to the whole PDB?
Right now it is a single process with an in-memory cache, which is fine for a demo. To scale up I would, in order: (1) precompute the metrics (they are deterministic per entry) for popular or all entries in a batch job and store them in a database such as Postgres, so the page reads a row and never parses a file; (2) put the RCSB access behind a shared cache (Redis or an object store of the mmCIF files) so many app instances do not each hit RCSB; (3) move heavier work, like parsing large entries, to a job queue with workers (Celery or RQ) and show progress in the UI; (4) serve the app from several containers behind a load balancer. For the whole PDB, RCSB also publishes bulk downloads, which are a much better source than crawling the API.

### 14. Why uv, a lockfile, and an exported `requirements.txt`?
`uv.lock` gives reproducible installs locally, in CI and in Docker (`uv sync --frozen`). Streamlit Community Cloud reads `requirements.txt` and not `uv.lock`, so I generate it from the lockfile with `uv export`, adding `-e .` so the `src/` package is installed. CI regenerates it and fails if it differs from the committed file, so it cannot drift. Using the lockfile also means the app that CI tested is the app that gets deployed.

### 15. What would you change if you had another week?
First, analyse the SEQRES sequence by fetching it from the polymer entity API, and let the user switch between observed and full. Second, treat Kyte-Doolittle windows more carefully across gaps (right now they slide across missing residues as if they were adjacent). Third, add an end-to-end browser test in CI; today the unit and `AppTest` tests do not exercise the actual 3D rendering, which I verified manually with headless Chromium. Fourth, shrink the Docker image, which is about 850 MB.

## Read the code in this order

1. `src/protein_dashboard/analysis.py` - the core: pure functions for composition, MW, pI, GRAVY and Kyte-Doolittle windows. Start here, as it is the most self-contained.
2. `tests/test_analysis.py` - how the analysis is tested with known values (ubiquitin MW, profile length, invalid input).
3. `src/protein_dashboard/rcsb_client.py` - the HTTP layer: retries, backoff, `Retry-After`, ID validation, error types. Then `tests/test_rcsb_client.py` for the mocked scenarios.
4. `src/protein_dashboard/structure.py` - mmCIF parsing: observed sequence, segments and gaps, modified residues, ligands, and `trim_mmcif`.
5. `src/protein_dashboard/viewer.py` - the Py3Dmol scene, the three colour modes, the heavy-structure guard, and the single `show_view` wrapper.
6. `src/protein_dashboard/app.py` - how it all comes together: query params as state, caching, error handling. Note how little logic it contains.
7. `Dockerfile` and `.github/workflows/ci.yml` - multi-stage build, non-root user, healthcheck, and the lint → test → docker pipeline.
