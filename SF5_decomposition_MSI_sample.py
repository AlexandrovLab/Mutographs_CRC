# Author: Marcos Diaz-Gay
# Date: Aug 3, 2024
# Local conda env: spa_0.1.4 (SigProfilerAssignment v0.1.4)

from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(samples="../data_for_figures/Sample_PD48980a_profile.all",
                   output="../data_for_figures/Sample_PD48980a_reconstruction",
                   input_type="matrix",
                   genome_build="GRCh38",
                   cosmic_version=3.4,
                   sample_reconstruction_plots = 'both')
