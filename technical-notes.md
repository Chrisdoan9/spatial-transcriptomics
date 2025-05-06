## DBiT-seq (Deterministic Barcoding in Tissue for spatial omics sequencing)

UMI = Unique Molecular Identifier  

•	It’s a short random sequence (like ACGT) added to each mRNA during library prep.  
•	It helps tag each original mRNA molecule uniquely before PCR amplification.  
•	After sequencing, we can collapse PCR duplicates to avoid overcounting.   

So:  

•	If 100 mRNA molecules from a gene are captured, and PCR amplifies them to 10,000 reads.  
•	The UMI helps ensure we count only 100, not 10,000.

  ## QC
  •	in_tissue == 1 → The spot does overlap the tissue (keep it)  
	•	in_tissue == 0 → The spot does not overlap the tissue (discard it)

## Niche
Cells belong to the same niche if they are surrounded by a similar composition of neighboring cell types — even if they are far apart in the tissue.
Seurat v5 https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2

## Coordinate Differences: Visium vs MERFISH

| **Feature**                   | **10x Visium**                                                                 | **MERFISH**                                                                 |
|------------------------------|--------------------------------------------------------------------------------|-----------------------------------------------------------------------------|
| **Coordinate system**        | Spot-based grid (≈ 55 µm diameter spots)                                       | Single-molecule / subcellular resolution                                   |
| **X–Y values**               | Each spot has a centroid X–Y coordinate (in microns or pixels)                 | Each mRNA molecule has an exact X–Y position (can also be in pixels/microns) |
| **Scale / resolution**       | Low resolution — each spot captures transcripts from multiple cells           | High resolution — each transcript is tied to a specific location within a cell |
| **Unit spacing**             | Spots are regularly spaced (fixed grid layout)                                | Coordinates are irregular / dense, based on actual molecular positions     |
| **Alignment to tissue image**| Coordinates are registered to histology image (H&E)                            | Coordinates are tied to high-resolution imaging (e.g., DAPI, protein stains) |
| **Typical usage**            | Map regional expression, cluster spots, infer cell types per spot             | Map molecular localization, segment individual cells, analyze subcellular expression |

## Pixel Space vs Microns
| **Term**         | **Meaning**                                                | **Example**                                                           |
|------------------|------------------------------------------------------------|------------------------------------------------------------------------|
| **Pixel space**  | Coordinates measured in pixels of an image (e.g., 512×512 image) | A gene signal at X = 200, Y = 300 pixels                              |
| **Microns (µm)** | Coordinates measured in real-world physical distance on the tissue | The same gene is located at X = 55.2 µm, Y = 82.7 µm from the origin of the slide |
