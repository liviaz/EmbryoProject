# EmbryoProject
This repository contains all the code used to gather and analyze data and generate figures for the embryo project.

It contains the following:

1. Labview code (Automated Pressure Control v20.vi) used to operate the system hardware (perform micropipette aspiration and record video). Requires drivers for Firgelli linear actuator (included in repo), as well as Thorlabs DCC1545M USB camera.

2. Matlab code to analyze the videos of micropipette aspiration, perform (mostly) automated detection of aspiration depth over time and calculate embryo mechanical parameters (depthDetectionAuto.m), visualize mechanical parameters (plotAllMouseEmbryosSoFar.m, plotAllHumanSoFar.m, plotAllMouseOocytesSoFar.m, saveNewMorphology.m), analyze cortical granule images (CGprocessingPipeline.m, analyzeMicroinjectionImages.m), and generate figures / perform further analysis (analyzeMicroinjectionData.m, makeMicroinjectionFigs.m, plotEmbryoToOocyteChange.m, makePaperFigs.m).

3. Shell scripts to run various aspects of the RNA-seq analysis pipeline (FastX Clipper/Trimmer, STAR aligner, HTSeq-Count, SamTools, RNASeQC) and prepare data for processing in R.

4. R scripts to run edgeR as well as ComBat, and generate figures for RNA-seq data. 


