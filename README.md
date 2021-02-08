# An R code to track functional extinction
R code to conduct extinction simulations and obtain patterns of functional rarity and functional traits as species goes extinct.

This is an R code to assess patterns of functional rarity and functional categories loss in biological communities at local and regional scales. This code runs three extinction scenarios at local and regional scales: scarce/restrict species first; abundant/widespread species first; and random extinction. The code assigns different weights to species based on their abundance and occurrence by extinguishing individuals at a local scale and populations at a regional scale, instead of species. The average functional rarity (distinctiveness and uniqueness) and the number of remaining functional categories are then calculated for the community after the extinction of each individual/population. We also provide null models of extinction locally, for each stream, and regionally, for all the streams, to compare with the pattern found in the other extinction scenarios. In the null models, the values of taxonomic rarity (and thus, the order of extinction) is shuffle 1000 times while maintaining functional rarity and functional trait values of each individual.
Locally, the code defines nine levels of extinction, 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80% and 90% of extinct individuals. Then, it computes the average functional rarity and the average number of functional categories for each stream in each extinction level. To identify differences between pairs of scenarios (scarcest extinct first; most abundant extinct first; and random extinction) at each of the nine extinction levels, a non-parametric Friedman-paired test and post-hoc pairwise comparisons are provided in the code. Regionally, it computes the confidence interval where 95% of data generated by the null model is confined.

For any inquires, feel free to contact me at lucasfcolares@gmail.com.
