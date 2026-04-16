# PathwayViz Shiny App — v1.0.2

Interactive visualization of enriched pathways or biological functions from
normalized proteomics or transcriptomics data.

Copyright 2026 RGM  
Licensed under the MIT License — see [License](#license) section below.

---

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Input Files](#input-files)
5. [App Structure](#app-structure)
6. [Tabs and Features](#tabs-and-features)
   - [Tab 1 — Bar Plot](#tab-1--bar-plot)
   - [Tab 2 — Chord Diagram](#tab-2--chord-diagram)
   - [Tab 3 — Heatmap](#tab-3--heatmap)
   - [Tab 4 — Boxplot](#tab-4--boxplot)
7. [Sidebar Controls Reference](#sidebar-controls-reference)
8. [Download Outputs](#download-outputs)
9. [Tips and Troubleshooting](#tips-and-troubleshooting)
10. [License](#license)

---

## Overview

PathwayViz is a self-contained R Shiny application that takes three CSV files
as input and produces four interactive, publication-ready visualizations:

| Visualization | What it shows |
|---|---|
| **Bar Plot** | Top N enriched pathways ranked by significance |
| **Chord Diagram** | Pathway–molecule connectivity and shared molecules |
| **Heatmap** | Z-score expression patterns per pathway/function |
| **Boxplot** | Per-group normalized abundance for a single protein/gene |

All plots are adjustable in real time via the sidebar and can be downloaded
as high-resolution image files.

---

## Requirements

### R version
R ≥ 4.2.0 recommended.

### R packages

| Package | Version tested | Purpose |
|---|---|---|
| `shiny` | ≥ 1.8 | App framework |
| `tidyverse` | ≥ 2.0 | Data wrangling and ggplot2 plotting |
| `pheatmap` | ≥ 1.0.12 | Clustered heatmap rendering |
| `circlize` | ≥ 0.4.15 | Chord diagram rendering |

Install all dependencies at once:

```r
install.packages(c("shiny", "tidyverse", "pheatmap", "circlize"))
```

---

## Installation

1. Clone or download this repository.
2. Place `app.R` in a folder of your choice.
3. Open R or RStudio and run:

```r
shiny::runApp("path/to/your/folder")
```

Or open `app.R` in RStudio and click **Run App**.

> **Upload limit:** The app accepts files up to **500 MB** per upload.

---

## Input Files

The app requires **three CSV files**. All must have a header row.

---

### 1. Normalized Data file

Rows = proteins or genes. Columns = sample measurements plus an identifier column.

| Column type | Description |
|---|---|
| **Gene/Protein ID** | One column with gene symbols or protein IDs (e.g. `Gene.Symbol`). Selected in sidebar under *Normalized counts or protein abundance — Column*. |
| **Sample columns** | One column per sample. Column names must match the Sample ID values in the Annotation file. |

**Minimal example:**

```
Gene.Symbol, Sample_A1, Sample_A2, Sample_B1, Sample_B2
TP53,        12.4,      11.9,      8.1,        7.8
EGFR,        9.0,       9.3,       14.2,       13.9
```

---

### 2. Sample Annotation file

Rows = samples. Maps each sample to a group.

| Column type | Description |
|---|---|
| **Sample ID** | Must match the column names in the Normalized Data file |
| **Group** | Experimental group label (e.g. `Control`, `Treated`) |

**Minimal example:**

```
SampleID,  Group
Sample_A1, Control
Sample_A2, Control
Sample_B1, Treated
Sample_B2, Treated
```

---

### 3. Enriched Pathways or Functions file

Output from an enrichment tool such as IPA, g:Profiler, Enrichr, or similar.
Rows = pathways or functions.

| Column type | Description |
|---|---|
| **Category** | Pathway or function name. Defaults to column 1. |
| **p-value / -log10(p)** | Significance column. Auto-detected by keyword matching (`fdr`, `adj`, `pval`, `-log`, etc.). |
| **Molecules** | Comma-separated list of gene/protein IDs in that pathway. Defaults to the last column. |

**Minimal example:**

```
Pathway,             p-value, Molecules
Cell cycle,          0.001,   TP53,EGFR,CDK2
Apoptosis,           0.005,   TP53,BCL2
DNA repair,          0.020,   BRCA1,TP53
```

> **p-value auto-detection:** The app checks whether values are between 0 and 1
> (raw p-values) or larger (already -log10 transformed) and applies the
> appropriate transformation automatically. A confirmation message is shown in
> the Bar Plot status line.

---

## App Structure

```
app.R
│
├── UI
│   ├── Sidebar (width = 3)
│   │   ├── File uploads
│   │   ├── Column mapping selectors
│   │   ├── Bar Plot Settings
│   │   ├── Chord Diagram Settings
│   │   ├── Heatmap Settings
│   │   └── Boxplot Settings
│   │
│   └── Main panel (width = 9)
│       ├── Tab 1 — Bar Plot
│       ├── Tab 2 — Chord Diagram
│       ├── Tab 3 — Heatmap
│       └── Tab 4 — Boxplot
│
└── Server
    ├── Reactive data loaders
    ├── Dynamic column-selector widgets
    ├── p-value transformation detection
    ├── Group colour inputs
    ├── Bar plot logic
    ├── Chord diagram logic
    ├── Heatmap logic
    └── Boxplot logic
```

---

## Tabs and Features

### Tab 1 — Bar Plot

Displays the top N enriched pathways ranked by significance as a horizontal bar chart.

| Feature | Detail |
|---|---|
| Ranking | Sorted descending by -log10(p-value) |
| p-value handling | Raw p-values are -log10 transformed automatically; pre-transformed values are used as-is |
| Top N | Configurable from 5 to 200 (default 20) |
| Bar colour | Any R colour name or hex code |
| Status line | Shows p-value column used, transformation applied, total pathways, and number displayed |
| Download | 300 dpi PNG, 10 × 8 inches |

---

### Tab 2 — Chord Diagram

Visualizes the connections between pathways and the molecules they contain.

| Feature | Detail |
|---|---|
| Pathways shown | Top N by significance (default 10) |
| Molecule filter | Minimum number of times a molecule must appear across pathways to be included (default 1) |
| Sector arc width | Proportional to number of connections — wider arc = more molecules (for pathways) or more pathways (for molecules) |
| Ribbon | One ribbon per pathway–molecule connection |
| Sector colours | Pathway family = tomato/red gradient; Molecule family = steelblue gradient. Shade within each family is alphabetical order only and carries no biological meaning |
| Inner radius | Adjustable via slider (0.1–0.7); smaller = more label space |
| Label font sizes | Separately configurable for pathways and molecules |
| Status line | Shows counts of pathways, molecules, and total connections |
| Download | 300 dpi PNG, 10 × 10 inches |

> **Reading the chord diagram:**
> - **Arc width** encodes connectivity — a wide pathway arc means many molecules;
>   a wide molecule arc means that molecule is shared across many pathways (hub protein).
> - **Colour family** (red vs. blue) distinguishes pathways from molecules.
> - **Colour shade** within a family is alphabetical only — it does not encode
>   significance, fold change, or any other variable.

---

### Tab 3 — Heatmap

Shows the Z-score normalized expression pattern for all proteins in a selected
pathway, with samples colour-annotated by group.

| Feature | Detail |
|---|---|
| Category selector | Dropdown populated from all unique pathway/function names in the enrichment file |
| Scaling | Row-wise Z-score (mean = 0, SD = 1 per protein) |
| Clustering | Rows (proteins) clustered; columns (samples) not clustered |
| Colour scale | Dark blue → white → dark red |
| Row labels | Shown when ≤ 100 proteins are displayed |
| Group annotation bar | One colour per group, configurable |
| Status line | Shows category name, gene list size, matched proteins, and sample columns used |
| Download (single) | Selected heatmap — 300 dpi PNG, 8 × 10 inches |
| Download (all) | All categories with ≥ 2 matched proteins — ZIP archive of PNGs |

---

### Tab 4 — Boxplot

Shows the raw normalized abundance distribution for a single selected
protein/gene, one box per experimental group.

| Feature | Detail |
|---|---|
| Protein selector | Dropdown populated from all unique molecules across all pathways in the enrichment file |
| Data source | Normalized Data file — same file used by the heatmap |
| Y-axis | Raw normalized values (no Z-score scaling) |
| Points | Individual samples shown as jittered dots overlaid on each box |
| Outlier display | Outlier points suppressed from the box geometry to avoid double-plotting with jitter |
| Group colours | Shared with the Heatmap Settings colour inputs |
| Status line | Shows selected protein, number of samples found, and group names |
| Download | 300 dpi PNG, 8 × 6 inches |

---

## Sidebar Controls Reference

### File Uploads

| Input | Description |
|---|---|
| Normalized Data (.csv) | Protein/gene abundance matrix |
| Metadata or Sample Annotation (.csv) | Sample-to-group mapping |
| Enriched Pathway or Functions (.csv) | Enrichment results with molecule lists |

### Column Mapping

| Input | Description |
|---|---|
| Category column | Pathway/function name column in the enrichment file |
| Molecules column | Comma-separated molecule list column (default: last column) |
| Gene/Protein ID column | Identifier column in the normalized data file |
| Sample ID column | Sample name column in the annotation file |
| Group column | Group/condition column in the annotation file |

### Bar Plot Settings

| Input | Default | Description |
|---|---|---|
| P-value column | Auto-detected | Column used for ranking pathways |
| Show top N pathways | 20 | Number of bars to display |
| Bar fill colour | steelblue | Any R colour name or hex code |
| Bar plot height (px) | 600 | Display height in the browser |

### Chord Diagram Settings

| Input | Default | Description |
|---|---|---|
| Top N pathways | 10 | Pathways included, selected by significance |
| Minimum appearances | 1 | Minimum pathway count for a molecule to appear |
| Inner circle size | 0.4 | Radius of the central hole (0.1–0.7) |
| Pathway label font size | 0.55 | cex value for pathway sector labels |
| Molecule label font size | 0.70 | cex value for molecule sector labels |
| Pathway sector colour | tomato | Base colour for pathway arcs |
| Molecule sector colour | steelblue | Base colour for molecule arcs |
| Chord plot height (px) | 750 | Display height in the browser |

### Heatmap Settings

| Input | Default | Description |
|---|---|---|
| Select category | First entry | Pathway/function to display |
| Colour per group | purple, darkorange, … | One text input per group; accepts any R colour |
| Heatmap height (px) | 700 | Display height in the browser |

### Boxplot Settings

| Input | Default | Description |
|---|---|---|
| Select protein / gene | First molecule | Dropdown from all molecules in the enrichment file |
| Boxplot height (px) | 500 | Display height in the browser |

> Group colours for the boxplot are the same inputs defined under
> **Heatmap Settings**.

---

## Download Outputs

Every visualization has two download buttons — one in the sidebar and one
below the plot in the main panel. Both produce identical files.

| Plot | Format | Size | DPI |
|---|---|---|---|
| Bar Plot | PNG | 10 × 8 in | 300 |
| Chord Diagram | PNG | 10 × 10 in | 300 |
| Heatmap (selected) | PNG | 8 × 10 in | 300 |
| Heatmap (all) | ZIP of PNGs | 8 × 10 in each | 300 |
| Boxplot | PNG | 8 × 6 in | 300 |

---

## Tips and Troubleshooting

**No proteins matched in heatmap**
> Check that the Gene/Protein ID column in the normalized data file uses the
> same identifiers (e.g. gene symbols) as the Molecules column in the
> enrichment file. Both are case-sensitive.

**Chord diagram is too cluttered**
> Increase *Minimum appearances* in Chord Diagram Settings to show only
> molecules shared across multiple pathways. Also try reducing *Top N pathways*.

**p-value column not auto-detected**
> Manually select the correct column from the dropdown. The status line in the
> Bar Plot tab confirms which transformation is being applied.

**Sample columns not found**
> The Sample ID column in the annotation file must contain values that exactly
> match column names in the normalized data file.

**Heatmap shows fewer proteins than expected**
> The app requires at least 2 matched proteins to draw a heatmap. Proteins
> present in the pathway list but absent from the normalized data file are
> silently skipped.

**Boxplot shows a flat line or single point**
> Only one sample was matched for that protein. Verify that sample IDs in the
> annotation file match the normalized data column names.

**Download produces a blank or corrupt file**
> Ensure the plot renders correctly in the browser before downloading.
> The download uses the same rendering function as the display.

---

## License

MIT License

Copyright 2026 RGM

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
