# Clustkmer
A Python module for identifying and analyzing conserved k-mers in protein sequences to improve structural classification accuracy

**Background**: Protein sequence alignment is crucial for understanding structure, function and evolution. Tools like BLAST enable rapid database searches, but sacrifice accuracy. MMseq2 allows sensitive searching for massive datasets. Foldseek uses MMseq2 techniques on 3D structure sequences, enabling fast structure database searches. However, its accuracy suffers due to using a sequence-based prefilter. Better prefilters considering protein structure are needed.

**Methods**: The exonuclease cluster C3Z734 was analyzed. 685 structures were obtained and aligned into a multiple sequence alignment (MSA). K-mer scoring matrices were generated from the MSA profiles. K-mers conserved in the cluster but rare in unrelated sequences were systematically identified. K-mer combinations were evaluated to optimize sensitivity and specificity for classifying sequences. The workflow was applied to 100 diverse clusters to assess the model’s performance.

**Results**: Conserved k-mers clustered spatially and in domains. With the top 10 k-mers requiring at least 2 matches, 99.3% sensitivity and 99.7% specificity were achieved for the test cluster. Across 100 clusters, average sensitivity and specificity were 85.6% and 83.1%. Less predictive clusters contained intrinsically disordered regions or repetitive structures lacking conserved motifs.

**Conclusion**: Conserved k-mers effectively identify protein structural motifs and classify sequences. The proposed workflow provides a promising approach to developing improved prefilters for protein structure search tools like Foldseek.
[GitHub](https://github.com)
