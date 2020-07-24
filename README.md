# maturation_networks
the PCIT_based_GRN_analyses.awk script describe steps to construct gene networks:
genes were used as nodes and significant connections (edges) between them were identified using the Partial Correlation and Information Theory (PCIT) algorithm (Reverter and Chan, 2008), considering all 48 samples. PCIT determinates the significance of the correlation between two nodes after accounting for all the other nodes in the network. 
differential connectivity:
two networks were created; one using 12 samples at T1 (pre-maturation) and a second using 36 samples at T2, T3 and T4 (post-maturation). The number of connections of each gene in each network was computed, making it possible to compare the same gene in the two networks to identify differentially connected genes (DCGs). A series of subnetworks were explored based on the top trio genes and the top regulators (TFs) based on their differential connectivity between pre-and post-maturation
