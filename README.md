# Host-Infection × Residual-Burden Coupling as a Delayed-Risk Constraint Graph

**Thesis #14** (series label NP-08). Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Can infection, marrow-suppression-like host stress, and residual-burden relapse risk be represented as a qualitative delayed-risk graph with explicit non-parameters, such that host-context constraints change hypothesis rank without becoming PK/PD coefficients or care pathways?

Four simple temporal networks share a frozen triple Θ = (CL, V, ke). Delay bounds are hypothesis syntax. Evidence windows live in a separate store. With residual-burden and relapse-like windows only, the lexicographic rank is H0, H1, H3, H2. Attaching an infection window and a marrow-stress-like window, still outside Θ, moves the order to H2, H1, H0, H3. H3 becomes empty. A volume-only sort still leaves H0 first, so the rank change is the unhosted-evidence term. A map from those windows onto k_inf, k_host, CL and ke is refused. The SHA-256 digest of Θ does not change.

The windows are declared. They are not a cohort. This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/delay_graph.py` | Seeded constraint graph and refusal (seed 20260921) |
| `sim/results.json` | Ranks and the refused proposal cited in Chapter Four |
| `sim/figures/` | Topologies, windows, rank change, refusal |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/delay_graph.py
python3 build_pdf.py
```

NumPy and Matplotlib are required for the sketch. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Host-infection × residual-burden coupling as a delayed-risk constraint graph [Internet]. Thesis #14 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-14-infection-residual-burden-delay-graph

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
