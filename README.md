# PicoHFAnalysis

A lightweight heavy-flavor analysis package for STAR **PicoDst** data.

## Features
* **JPsiMaker** – di‐electron J/ψ reconstruction
* **DMesonMaker** – D⁰ → Kπ (and \bar D⁰) reconstruction
* **HFEMaker** – non‐photonic electron (NPE) spectra & QA
* Central steering via `PicoDstAnalysisManager`
* ROOT histograms written to a single output ROOT file (`hfOut.root` by default)

## Building
```bash
mkdir build && cd build
cmake ..
make -j4
```
Ensure `$ROOTSYS` is set and STAR libraries (`StRoot` and `StPicoDst`) are in your environment (`$STAR`).

## Running
Edit or create a text file `pico.list` with full paths to your PicoDst files.
Run from a ROOT session:
```C++
.L macros/runPicoHF.C
runPicoHF("pico.list", "hfOut.root", 1000000);
```
Histograms will be appended to `hfOut.root`.

## Output objects
| Maker        | Histograms                             |
|--------------|----------------------------------------|
| JPsiMaker    | `hJPsiInvMass`, `hJPsiPtY`             |
| DMesonMaker  | `hD0Mass`, `hD0PtY`                    |
| HFEMaker     | `hNPEPt`, `hEoverPvsP`                 |

Use standard ROOT to plot signals after background subtraction / fits.
