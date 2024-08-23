# Author: Marcos Diaz-Gay
# Date: Jul 30, 2024
# Local conda env: spa_0.1.4 (SigProfilerAssignment v0.1.4)

import sigProfilerPlotting as sigPlt
import os

groups = ['MSS']
subcontexts = ['SBS288', 'ID83', 'CN48', 'SV32', 'DBS78']
plottypes = ['288', '83', '48', '32', '78']
comparisons = ['country', 'age_eo']

for comparison in comparisons:
    for s,p in zip(subcontexts, plottypes):
        if s == 'SBS288':
            sigPlt.plotSBS('../data_for_figures/Avg_mut_profiles/Matrix_' + s + '_' + comparison + '.all',
                            '../data_for_figures/Avg_mut_profiles/',
                            comparison,
                            p,
                            percentage=True)
            sigPlt.plotSBS('../data_for_figures/Avg_mut_profiles/Avg_MSS_' + s + '.all',
                            '../data_for_figures/Avg_mut_profiles/',
                            'MSS',
                            p,
                            percentage=True)
            sigPlt.plotSBS('../data_for_figures/Avg_mut_profiles/Avg_MSI_' + s + '.all',
                            '../data_for_figures/Avg_mut_profiles/',
                            'MSI',
                            p,
                            percentage=True)                            
        if s == 'DBS78':
            sigPlt.plotDBS('../data_for_figures/Avg_mut_profiles/Matrix_' + s + '_' + comparison + '.all',
                            '../data_for_figures/Avg_mut_profiles/',
                            comparison,
                            p,
                            percentage=True)
            sigPlt.plotDBS('../data_for_figures/Avg_mut_profiles/Avg_MSS_' + s + '.all',
                            '../data_for_figures/Avg_mut_profiles/',
                            'MSS',
                            p,
                            percentage=True)
            sigPlt.plotDBS('../data_for_figures/Avg_mut_profiles/Avg_MSI_' + s + '.all',
                            '../data_for_figures/Avg_mut_profiles/',
                            'MSI',
                            p,
                            percentage=True)
        if s == 'ID83':
            sigPlt.plotID('../data_for_figures/Avg_mut_profiles/Matrix_' + s + '_' + comparison + '.all',
                            '../data_for_figures/Avg_mut_profiles/',
                            comparison,
                            p,
                            percentage=True)
            sigPlt.plotID('../data_for_figures/Avg_mut_profiles/Avg_MSS_' + s + '.all',
                            '../data_for_figures/Avg_mut_profiles/',
                            'MSS',
                            p,
                            percentage=True)
            sigPlt.plotID('../data_for_figures/Avg_mut_profiles/Avg_MSI_' + s + '.all',
                            '../data_for_figures/Avg_mut_profiles/',
                            'MSI',
                            p,
                            percentage=True)                
        # if s == 'CN48':
        #     os.chdir('../data_for_figures/Avg_mut_profiles/')
        #     sigPlt.plotCNV('Matrix_' + s + '_' + comparison + '.all',
        #                     './',
        #                     comparison,
        #                     percentage=True)
        #     os.chdir('../../Mutographs_CRC')
        # if s == 'SV32':
        #     os.chdir('../data_for_figures/Avg_mut_profiles/')
        #     sigPlt.plotSV('Matrix_' + s + '_' + comparison + '.all',
        #                     './',
        #                     comparison,
        #                     percentage=True)
        #     os.chdir('../../Mutographs_CRC')
# Names changed manually in the matrices
