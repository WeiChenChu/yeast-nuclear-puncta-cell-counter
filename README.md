# yeast-nuclear-puncta-cell-counter

### Fiji/ InageJ Macro: Cellpose + CLIJ2 pipeline for positive-cell detection

This repository contains an ImageJ/Fiji macro for batch-processing `.czi` images to:
1) segment cells using Cellpose (cyto3, using the cellpose wrapper),
2) detect puncta/green spots (C2),
3) detect red membrane/nuclear boundary signal (C1),
4) classify cells as "positive" based on a merged mask criterion,
5) export label images and (optionally) ROIs for positive cells.

#### Features
- Recursive batch processing of folders/subfolders
- GPU-accelerated filtering + masking using CLIJ2
- Cell segmentation via Cellpose (cyto3)
- Excludes edge-touching and overexposed cells prior to quantification
- Outputs:
  - `*_label.tif` (filtered cell labels)
  - `*_positive_label.tif` (visualization / overlay output)
  - `*_ROI.zip` (ROIs for positive cells; only saved if positives exist)
- Console log per image: total cells, positive cells, positive ratio

#### Dependencies
- Fiji (ImageJ) with:
  - Bio-Formats Importer
  - CLIJ + CLIJ2
  - BioVoxxel 3D Box
  - PTBIOP (Cellpose Wrapper, a Python-based Cellpose environment requireded)

> Important: You must set your Cellpose environment path in the macro
> (`env_path=...`). See **Configuration** below.

#### Input requirements / assumptions
- Images are imported as Hyperstacks via Bio-Formats.
- The macro expects **3 channels**:
  - C1: red signal (membrane / nuclear boundary; used for mask)
  - C2: green signal (spots; used for maxima detection)
  - C3: used as Cellpose input (cyto3)
- The macro currently uses a **fixed crop ROI**, to remove the uneven illumination part of the image:
  - `x=315, y=0, width=1300, height=1216`
  - Adjust this ROI to match your acquisition layout.

#### Installation
1. Install Fiji: https://fiji.sc/
2. Install CLIJ2, BioVoxxel, and PTBIOP (Cellpose wrapper):
   - Fiji > Help > Update...
   - Manage update sites: enable **CLIJ** / **CLIJ2** / **BioVoxxel 3D Box** /**PTBIOP**/

#### Configuration (must edit before running)
Open the macro and set:
- File suffix: default `.czi`
- Cell size filter:
  - `cell_minimum_size`
  - `cell_maximum_size`
- Denoise / background:
  - `median_radius`
  - `TopHat_radius`
- Signal detection:
  - `green_signal_prominence`
  - `red_segment_method` (e.g., Otsu)
  - `number_of_dilations_and_erotions`
- Crop ROI:
  - `makeRectangle(x, y, w, h);`
- Cellpose environment + model:
  - `env_path=...` (path to your conda env)
  - `model=cyto3`
  - `diameter=55` (tune for your cells)
  - `model_path=...` if using a custom model

#### Usage
1. Launch Fiji.
2. Run the macro (Plugins > Macros > Run...).
3. Select:
   - Input directory (can contain subfolders)
   - Output directory
   - File suffix (default `.czi`)
4. The macro will process all matching files recursively.

#### Output files
For each input image `NAME.czi`, the macro produces in the output directory:
- `NAME_label.tif`  
  Filtered cell label map after edge exclusion and size filtering.
- `NAME_positive_label.tif`  
  Visualization output of positive labels (and contrast/channel LUT changes).
- `NAME_ROI.zip` (optional)  
  Saved only when at least 1 positive cell is detected.

#### How "positive" is defined (high level)
- Green spots: maxima detection on (C2 - median filtered C2)
- Red mask: top-hat + auto-threshold + morphological closing on C1
- Merge: `mask_merge = green_spot + red_mask`
- Positive cells are those whose label region satisfies the merged-mask criterion
  (see macro section near `positive_map` / `positive_label`).

#### Exclusion criteria
- Edge-touching cells were excluded (labels on image borders are removed).
- Overexposed cells (likely dead cells) were excluded (labels containing saturated/overexposed pixels are removed).

#### Known limitations
- Hard-coded crop ROI (dataset-specific).
- Channel title assumptions after "Split Channels" may differ across Fiji versions.

### References
#### CLIJ / CLIJ2
- GitHub: https://github.com/clij/clij2  
- Documentation: https://clij.github.io/  
- Robert Haase, Loic Alain Royer, Peter Steinbach, Deborah Schmidt, Alexandr Dibrov, Uwe Schmidt, Martin Weigert, Nicola Maghelli, Pavel Tomancak, Florian Jug, Eugene W Myers. CLIJ: GPU-accelerated image processing for everyone. *Nat Methods* 17, 5–6 (2020). doi:10.1038/s41592-019-0650-1  
  Link: https://doi.org/10.1038/s41592-019-0650-1  
- Daniela Vorkel, Robert Haase. GPU-accelerating ImageJ Macro image processing workflows using CLIJ. *arXiv preprint*  
  Link: (add arXiv URL here if you have the exact identifier)  
  Code: https://github.com/clij/clij
- Robert Haase, Akanksha Jain, Stéphane Rigaud, Daniela Vorkel, Pradeep Rajasekhar, Theresa Suckert, Talley J. Lambert, Juan Nunez-Iglesias, Daniel P. Poole, Pavel Tomancak, Eugene W. Myers. Interactive design of GPU-accelerated Image Data Flow Graphs and cross-platform deployment using multi-lingual code generation. *bioRxiv preprint*  
  Link: (add bioRxiv URL/DOI here if you have it)  
  CLIJ2 code: https://github.com/clij/clij2

#### Cellpose
- GitHub: https://github.com/MouseLand/cellpose  
- Documentation: https://cellpose.readthedocs.io/  
- Stringer, C., Wang, T., Michaelos, M., & Pachitariu, M. (2021). Cellpose: a generalist algorithm for cellular segmentation. *Nature Methods*, 18(1), 100–106.  
  Link: https://www.nature.com/articles/s41592-020-01018-x  
- Stringer, C. & Pachitariu, M. (2025). Cellpose3: one-click image restoration for improved segmentation. *Nature Methods*.  
  Link: https: https://www.nature.com/articles/s41592-025-02595-5
  Code: https://github.com/MouseLand/cellpose

#### Fiji / ImageJ
- Schindelin, J. et al. Fiji: an open-source platform for biological-image analysis. *Nature Methods* (2012).  
  Link: https://www.nature.com/articles/nmeth.2019  

#### Bio-Formats (OME) — for reading CZI files
- Linkert, M. et al. Metadata matters: access to image data in the real world. *Journal of Cell Biology* (2010).  
  Project: https://www.openmicroscopy.org/bio-formats/  

#### BioVoxxel 3D Box (Fiji plugin)
- GitHub: https://github.com/biovoxxel/3D_Box
- BioVoxxel 3D Box (provides “Labels to 2D Roi Manager” in Fiji).  
