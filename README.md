# Clustkmer
A Python module for identifying and analyzing conserved k-mers in protein sequences to improve structural classification accuracy

**Background**: Protein sequence alignment is crucial for understanding structure, function and evolution. Tools like BLAST enable rapid database searches, but sacrifice accuracy. MMseq2 allows sensitive searching for massive datasets. Foldseek uses MMseq2 techniques on 3D structure sequences, enabling fast structure database searches. However, its accuracy suffers due to using a sequence-based prefilter. Better prefilters considering protein structure are needed.

**Methods**: The exonuclease cluster C3Z734 was analyzed. 685 structures were obtained and aligned into a multiple sequence alignment (MSA). K-mer scoring matrices were generated from the MSA profiles. K-mers conserved in the cluster but rare in unrelated sequences were systematically identified. K-mer combinations were evaluated to optimize sensitivity and specificity for classifying sequences. The workflow was applied to 100 diverse clusters to assess the model’s performance.

**Results**: Conserved k-mers clustered spatially and in domains. With the top 10 k-mers requiring at least 2 matches, 99.3% sensitivity and 99.7% specificity were achieved for the test cluster. Across 100 clusters, average sensitivity and specificity were 85.6% and 83.1%. Less predictive clusters contained intrinsically disordered regions or repetitive structures lacking conserved motifs.

**Conclusion**: Conserved k-mers effectively identify protein structural motifs and classify sequences. The proposed workflow provides a promising approach to developing improved prefilters for protein structure search tools like Foldseek.

## Papers 
1. Fast and accurate protein structure search with Foldseek

    Citation: [van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. Fast and accurate protein structure search with Foldseek. Nat Biotechnol. 2024 Feb;42(2):243-246. doi: 10.1038/s41587-023-01773-0. Epub 2023 May 8. PMID: 37156916; PMCID: PMC10869269.] (https://pubmed.ncbi.nlm.nih.gov/37156916/)

    Summary: This paper introduces Foldseek, a tool for fast and accurate protein structure searches. The authors demonstrate its efficiency in large-scale protein database searches, highlighting significant improvements in speed and accuracy over existing methods.

2. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets

    Citation: [Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol. 2017 Nov;35(11):1026-1028. doi: 10.1038/nbt.3988. Epub 2017 Oct 16. PMID: 29035372.] (https://www.nature.com/articles/nbt.3988)

    Summary: This paper presents MMseqs2, a tool designed for sensitive and efficient protein sequence searching. The authors highlight its capability to handle massive datasets, providing a significant advantage for large-scale proteomics and genomics research.


