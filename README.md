# varveModel
Rigorous varve count age modeling, developed partly for the manuscript "A multicore 2250-year varve chronology from Eklutna Lake, Alaska, with a multifaceted approach to assessing age uncertainty", by David Fortin, Nore Praet, Nicholas P. McKay, Darrell S. Kaufman, Britta J.L. Jensen, Peter J. Haeussler, Casey Buchanan, and Marc De Batist. 

To reproduce the results of the manuscript, run the "varveModel.R" script. Be warned, the code is computationally intensive, and is best suited for High-Performance Computing. A parallelized version is under development.

To run the model on a different dataset, you'll need to create input files that mirror the Eklutna Lake data in the "Eklutna" folder, and run the "importCoreData.R" script, and adjust references as needed.

Although the model is fully-functional and ready for scientific application, it is still under development, and there are multiple avenues by which it could be improved, both conceptually, and in terms of the technical implementation. Please feel from to fork, or get in touch if you're interested in collaboration.

[![DOI](https://zenodo.org/badge/136086485.svg)](https://zenodo.org/badge/latestdoi/136086485)

## Coneptual overview of algorithm
This model effectively integrates both a sedimentological understanding of the likelihood of the correct delineation of varves and marker layers with the constraint that all sites must have the same number of years (varves) between marker layers, even though that number is unknown. The model takes advantage of the replication and marker layers at the site to quantify the likelihood of over- and under-counting, as well as the ultimate uncertainty in varve year as a function of depth in each core. Consequently, the model is completely independent from 14C or other age control, which serve as appropriate validation of the model and its uncertainty quantification. 

There are three primary inputs to the model. The first is varve thickness. The second is prior estimates of over- and under-counting probabilities for each varve in each core, that is, the likelihood that a counted varve should actually be two varves (under-counted) or that an identified varve should actually be part of an adjacent varve (over-counted). The third input is the position of the marker layer layers found in two or more cores, and prior probability estimates that each marker layer is properly identified and isochronous.

Upon initialization, the algorithm heuristically updates the over- and under-counting priors based on observed differences in the counts between marker layers. The heuristic updates are only used to increase miscounting probabilities, and thus the originally assigned priors ultimately serve as minimum miscounting estimates. The model then proceeds to simulate 1000 iterations of the varve counts for each core. In each simulation, there several steps. First, we sample from the marker-layers proportionally to the prior estimates of each marker layer to select a subset of marker-layers that are treated as correctly identified and isochronous for this simulation. Because this set of marker layers varies with each iteration, the impact of misidentification of the marker-layers is included in this analysis. Then, for each segment (the interval of varves between a pair of marker-layers, the probability of miscounting the number of varves between marker layers is used to simulate a “true” number of varves in that segment in each core. This simulation serves as the likelihood function in the simulation, and this steps repeats until, by random, the number of simulated years between marker layers in all cores is identical. This step can be numerically challenging, especially as the number of cores increases. For example, consider a segment with 100 counted years: in five cores, each year could be over-counted, under-counted or correct for each core. The total number of permutations is therefore 5 * (3100) = 2.5 x 1048. Consequently, it is not possible to numerically search the full parameter space, and even finding a reasonable number of valid solutions requires some accommodations, including the heuristic adjustment to the prior described above. The second accommodation is, if the algorithm fails to find a solution for all cores, it iteratively removes one core, selected by random, from the simulation and attempts again. This continues until a solution is found, or until only two cores remain. Consequently, most iterations include gaps between marker layers for each core, but across the-1000 member ensemble there is robust coverage. 

The 1000-member ensemble is then used to quantify the uncertainty in depth as a function of varve year, and can be transposed to estimated uncertainty in varve year as a function of depth in each core. By design, and largely due to the heuristic update of the miscounting probabilities, the central tendency of the model output resembles the average counts across the cores, as simulating the same number of varves is currently the only constraint on likelihood. Additional constraints on likelihood, such as independent age control and correlation of varve thicknesses across cores, are planned improvements to the model, but beyond the scope of this work.
