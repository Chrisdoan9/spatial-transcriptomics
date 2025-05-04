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
