# varveModel
Rigorous varve count age modeling, developed partly of the manuscript "A multicore 2250-year varve chronology from Eklutna Lake, Alaska, with a multifaceted approach to assessing age uncertainty", by David Fortin, Nore Praet, Nicholas P. McKay, Darrell S. Kaufman, Britt J.L. Jensen, Peter J. Haeussler, Casey Buchanan, and Marc De Batist. 

To reproduce the results of the manuscript, run the "varveModel.R" script. Be warned, the code is computationally intensive, and is best suited for High-Performance Computing. A parallelized version is under development.

To run the model on a different dataset, you'll need to create input files that mirror the Eklutna Lake data in the "Eklutna" folder, and run the "importCoreData.R" script, and adjust references as needed.

Although the model is fully-functional and ready for scientific application, it is still under development, and there are multiple avenues by which it could be improved, both conceptually, and in terms of the technical implementation. Please feel from to fork, or get in touch if you're interested in collaboration.

