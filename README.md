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
cons
```

## Running

```C++
.L macros/runPicoHF.C
runPicoHF("pico.list", "hfOut.root", 1000000);
```
Histograms will be appended to `hfOut.root`.
