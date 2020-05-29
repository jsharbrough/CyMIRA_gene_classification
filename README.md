**CyMIRA_gene_classification**
==========================


Plant genome gene categorization pipeline incorporating:
<ul><li>CyMIRA (http://cymira.colostate.edu)</li>
<li>Homology to Arabidopsis inferred by Orthofinder v2+ (https://github.com/davidemms/OrthoFinder)</li>
<li>De novo targeting predictions from four different sub-cellular targeting prediction programs
    TargetP, iPSORT, Predotar, LOCALIZER</li></ul>


Required dependencies not included in this repository:
<ul><li>Orthofinder v≥2</li>
<li>TargetP</li>
<li>iPSORT</li>
<li>Predotar</li>
<li>LOCALIZER</li>
</ul>

Example job submission bash scripts have been provided for use (assuming SLURM job scheduler, but all programs should be able to be run locally on any linux-capable machine (run time might be affected though). TargetP appears to be particularly sensitive to computer architecture.


Please address queries/comments to Joel Sharbrough (jsharbro [at] gmail.com). If you can't get it to work, I'd be happy to assist and/or run it for you.
