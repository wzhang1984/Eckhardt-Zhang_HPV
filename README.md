# Eckhardt-Zhang_HPV
Codes for the paper: Multiple Routes to Oncogenesis are Promoted by the Human Papillomavirus-Host Protein Network

## Identification of Mutated Genes in Each Tumor
Genes were classified as wild type (0) or altered (1) in each of the 799 tumors with alterations defined as follows. Most oncogenes (e.g., EGFR) were considered altered (activated) if impacted by a missense mutation, in-frame indel or copy number amplification. For the subset of oncogenes typically altered only by amplification (CCND1, LMO1, MDM2, MDM4, MYC, MYCL, MYCN, NCOA3, NKX2-1 and SKP2), only copy number amplifications were considered as alterations and not SNVs or indels. All other genes including tumor suppressors (e.g., CDKN2A) were considered altered (inactivated) if there was any type of non-silent mutation or a copy number deletion.

## Difference in Mutation Rate by HPV Status
For each gene, we fit a logistic regression model of its alteration state g (0 = wild type; 1 = altered) as a function of v, the HPV status (1 = HPV(+); 0 = HPV(−)), controlling for the impact of the cancer types t (HNSCC vs. cervical cancer) and cohorts c (TCGA vs. University of Chicago) as co-variates.
 
 To assess whether the HPV status is significantly associated with a gene’s alterations, the likelihood of the complete model was compared to that of a simple model under null hypothesis of no association. The deviance D (i.e., log likelihood ratio) between the two nested models was calculated. For the genes with increased mutation rate in HPV(–), D was used as the test statistic for the following network propagation.
  
## Network Propagation
Network propagation was used to identify the clusters of genes within the Reactome functional interaction network (ReactomeFI) that are enriched in both HPV interactors and genes with an increased mutation rate in HPV(−). MiST scores (without threshold), which  characterize the confidence of HPV-human PPIs, or the deviances D, which characterize the significance of the increased mutation rates in HPV(−) tumors, were separately propagated across the ReactomeFI network based on a random walk model (equivalent to heat diffusion) with a restart probability of 0.5. After convergence, the score of each gene (temperature) represents its network proximity to HPV protein binding or increased mutation rates in HPV(–). The product of two propagated scores for each gene was used as the test statistic. To estimate the expected background of this product, 20,000 permutations were performed by randomly reassigning MiST scores or deviances D, performing the network propagation and calculating the product of the two propagated scores. The goal was to produce a distribution of this product for each gene under the null hypothesis that the network neighborhood of this gene is not significantly enriched in both HPV binding and increased genetic alterations in HPV(−), while retaining the network structure. This null distribution was indexed with the observed product to obtain an empirical p-value characterizing the significance in both propagations. Using the Storey approach, q-values were then calculated and a threshold of FDR ≤ 25% was applied.
