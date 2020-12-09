import os
import sys
import numpy
import scanpy.external as sce
import diffxpy.api as de
import matplotlib.pyplot as plt
import logging
import seaborn
import pickle
import pandas as pd
import scanpy as sc
from anndata._core.anndata import AnnData
import numpy as np
import csv
import operator
import gseapy
from gseapy.plot import gseaplot, heatmap
import lostdata as LSD
import leopard as lp
from apetype.tasks import SKIP, SKIPCACHE
from apetype.reports import TaskReport, TaskSection
from collections import OrderedDict

# Utility functions
def reformat_gs_table(table, maxlength=30, roundto=2): 
    table['Adjusted P-value'] = table['Adjusted P-value'].apply(lambda x: np.format_float_scientific(x, roundto))
    table.Term = table.Term.apply(lambda x: x if len(x)<=maxlength else x[:maxlength-3]+'...')
    return table

def reformat_prerank_table(table, maxlength=30, roundto=2): 
    table.fdr = table.fdr.apply(lambda x: np.format_float_scientific(x, roundto))
    table.index = [x if len(x)<=maxlength else x[:maxlength-3]+'...' for x in table.index]
    return table

# Reloaded gseapy enrichr results
class Reloaded_Enrichr_Results:
    def __init__(self, location):
        import os, glob
        self.results = pd.concat(
            [pd.read_table(f) for f in glob.glob(os.path.join(location, '*.txt'))]
        )

class Analysis(TaskReport):
    verbose = True
    cache = './cache'

    # Report
    title = 'RPEMYCN single cell analysis'
    author = 'Christophe Van Neste'
    @property
    def outfile(_):
        return os.path.join(_.savedir, 'report')
    
    # CLI
    all_treatments: bool = True
    data: str = "./Preprocessing and exploratory analysis/RPE_MYCN_10x_after_seurat2.h5" # after mail VO 20/07/20
    savedir: str = './analysis_results_cvn'
    output_report: bool = False
    filter_cells: int = 1000 # cells with minimum # genes reported
    filter_genes: int = 3    # genes reported in minumum # of cells
    filter_mt: float = 12    # max percentage of mitochondrial reads
    traject_nbr: int = 2     # number of neighbours for trajectory calculation (2 or 3) more leads to many cluster conn
    diffxpy: bool = True     # run diffxpy on ctrl, 22, 23, 24 as categorical and 22, 23 and 24 quantitatively
    bulk_de_thresh: float = .01 # threshold for diff expressed genes in bulk data
    gsea_perms: int = 100    # permutations for gsea prerank
    exclude_cell_cycle_genes: bool = False # in bulk gene sigs, exclude the genes used for determining cell phase
    
    def configure(_) -> SKIPCACHE(type(None)): #TODO type None always skip cache
        # Scanpy settings
        sc.settings.figdir = os.path.join(_.savedir, 'figs_scanpy')
        _.gseadir = os.path.join(_.savedir, 'gsea')
        if not os.path.exists(_.savedir):
            os.makedirs(_.savedir)
            os.mkdir(_.gseadir)
            os.mkdir(sc.settings.figdir)


    def quality_control(_, print) -> TaskSection:
        '''As to my knowledge this is the 3rd report in total on the RPEMYCN single cell data set. 
        Two previous reports having been made by Vlodomir. In the previous reports there was a
        narrow focus on finding the senescence groups/cells. In this report I tried to do a more
        general cell cycle analysis.
        
        In Vlodomir's first report there the data had been analysed per treatment duration
        separately, whereas in his second report only ctrl and treatment 22 were taking
        in consideration and contrasted with each other. In this report, I have again
        considered all the data, but separate in each treatment group.
        
        If you look back to the umap's of Vlodomir with all data together, you can see
        overlapping circular structures. One circle is more or less the cell cycle of the
        ctrl condition, whereas the other circle is the mixed cell cycles of 22, 23 and 24.
        When you cluster all this data together, you will get confused clusters containing
        cells from the different conditions, in different cell phases. My intuition at this
        point is that it is better to first separate cells in their phases separately, and
        then doing a comparison across treatments for the same cell phase (this however I have
        yet to do.)
        
        A first remark on the quality control of the data. If you look at Figure 3, the ctrl
        stands out, in not having a lot of cells with a smaller amount of total count of genes.
        Looking at "Fraction of Reads in Cells" in Vlodomir's first report, the ctrl stands out
        in only having 58% opposed to a minimum of 93% for the treatment conditions. According
        to 10x, there are two possible explanations: 1) a high amount of ambient RNA coming
        from dead or lysed cells 2) maybe more than expected cells with a low amount of RNA
        and a subsequent error in their calling heuristic. 1 does not seem plausible, unless
        something went wrong with that sample that would explaing a higher amount of dead and
        lysed cells in comparison with the treatment conditions. Could 2 be true and relevant?
        More control cells with a lower amount of RNA?
        
        Most likely the preprocessing was with all data together, which could have created
        a disadvantage for control cells to be called if it is due to 2). The calling heuristic
        allows for a 10fold variance in RNA expression. A condtion with a lower amount of RNA
        could get a negative bias because of this. It is possible to rerun the preprocessing
        with: cellranger reanalyze --force-cells. If reanalyzing this way, we do not see this
        strong of a difference, this would make the hypothesis that the control cells have
        less RNA more plausible.
        
        If you look at figure 4, again the ctr stands out, in not having as many high gene count and
        high total count cells. This could also be compatible with MYCN leading to a higher per cell
        amount of transcription. This would also be compatible with what we see in the timeseries data:
        initially increasing library sizes after MYCN translocation. Of course it could also be
        that something went wrong with the single cell ctrl sample in handling or sequencing.
        '''
        print(_)
        return {}

    def methods(_) -> SKIPCACHE(TaskSection):
        """The RPEMYCN single cell analysis was performed with the Python
        scanpy package (v1.5.1). Cells that reported less than 1000
        genes were filtered, so were genes that were reported in only
        1 or 2 cells. Cells were allowed to have up to 12%
        mitochondrial reads. For determining the differential
        expression of single cell genes, the Python package diffxpy
        (v0.7.4) was used. To calculate the bulk RPEMYCN signatures in
        the single cell umap spaces, only genes were considered with
        an adjusted p-value < 0.01. Single cell library sizes were
        adjusted to 10000 and transformed by natural logarithm after
        addition of 1. Umaps with 3 components were produced by taking
        10 neighbors in a 25 component PCA space. Cell phase was
        determined with the S and G2M gene sets used in (PMID:
        31294801). In figures with all treatments together, cell phase
        was determined across all data, in figures were only one
        treatment is displayed, cell phase was determined for each
        cell only within that treatment data. Gene senescence
        signature scores were calculated as the number of genes that
        had expression in the cells. For the bulk RPEMYCN signatures
        in the single cell, the average rank of the genes in the set
        was calculated. Full scripts can be found at
        https://github.com/dicaso/apetype/examples/RPEMYCN.

        """
        return {'clearpage':True}
    
    def msigdb(_) -> dict:
        return LSD.get_msigdb6()
        
    def genesets(_) -> dict:
        return {
    'S': ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1',
               'UHRF1',
               'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3',
               'MSH2',
               'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN',
               'POLA1',
               'CHAF1B', 'BRIP1', 'E2F8'],
    'G2M': ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67',
                'TMPO', 'CENPF',
                'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B',
                'GTSE1', 'KIF20B',
                'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2',
                'CDCA8', 'ECT2',
                'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3',
                'CBX5', 'CENPA'],
    'GO_senescence': ['HMGI-C', 'HMGA1', 'HMGA2'],
    'GO_senescence2': ['HMGI-C', 'HMGA1', 'HMGA2', 'KRAS', 'PAWR', 'YPEL3', 'SIRT1'],
    'Courtois_senescence_trigers': ['BMPR1A', 'ZPF277', 'HMGI-C', 'ARNTL', 'PNPT1', 'ZHX2', 'F7I934'],
        }

    def senescence_table(_) -> pd.DataFrame:
        # https://genomics.senescence.info/cells/ https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3531213/
        #senescence_table = pd.read_html('https://genomics.senescence.info/cells/query.php?show=5&sort=1&page=1')[0]
        return pd.read_html('genomics_senescence_info.html')[0]
        
    def gseasets(_, print) -> list:
        _.failed_gsea = 0
        gseasets = [
            'ChEA 2016','GO Biological Process 2018',
            #'TF Perturbations Followed by Expression',
            'TRRUST Transcription Factors 2019'
        ]
        gseasets = [gs.replace(' ','_') for gs in gseasets]
        print('The following Enrichr genesets were used for GSEA:', ', '.join(gseasets))
        return gseasets

    def diffgenesbulk(_, bulk_de_thresh, exclude_cell_cycle_genes, genesets) -> dict:
        import zipfile
        from io import TextIOWrapper 
        #zf = zipfile.ZipFile('../RPE_MYCN_timeseries/rereportrpemycneranalysiscvn.zip')
        zf = zipfile.ZipFile('../RPE_MYCN_timeseries/diffexprgenes_BDC.zip')
        de_genes = {
            os.path.basename(f.filename)[:-4]: pd.read_table(
                TextIOWrapper(zf.open(f), encoding='latin-1')
            )
            for f in zf.infolist()
            if f.filename.endswith('.xls')
            and 'top' not in f.filename
            and '__MACOSX' not in f.filename
        }
        if bulk_de_thresh:
            de_genes = {
                comp:de_genes[comp].loc[de_genes[comp]['adj.P.Val'].str.replace(',','.').astype(float)<=bulk_de_thresh]
                for comp in de_genes
            }
        if exclude_cell_cycle_genes:
            de_genes = {
                comp:de_genes[comp].loc[~de_genes[comp].external_gene_name.isin(genesets['G2M']+genesets['S'])]
                for comp in de_genes
            }
        de_genes_sets = {}
        
        # All
        de_genes_sets['all'] = {
                comp:de_genes[comp]['external_gene_name']
                for comp in de_genes
            }
        de_genes_sets['all']['intersection'] = set.intersection(*[
            set(de_genes_sets['all'][comp])
            for comp in de_genes
        ])
        
        # Up
        de_genes_sets['up'] = {
                comp:de_genes[comp]['external_gene_name'][de_genes[comp].logFC.str.replace(',','.').astype(float)>0]
                for comp in de_genes
            }
        de_genes_sets['up']['intersection'] = set.intersection(*[
            set(de_genes_sets['up'][comp])
            for comp in de_genes
        ])
        
        # Down
        de_genes_sets['down'] = {
                comp:de_genes[comp]['external_gene_name'][de_genes[comp].logFC.str.replace(',','.').astype(float)<0]
                for comp in de_genes
            }
        de_genes_sets['down']['intersection'] = set.intersection(*[
            set(de_genes_sets['down'][comp])
            for comp in de_genes
        ])
        
        return de_genes_sets

    def adata(_, report, genesets) -> AnnData:
        adata = sc.read(_.data)
        # Set time as indicated by cell barcode name TODO check with nxtgnt
        adata.obs.time = [i.split('-')[-1] for i in adata.obs.index]
        #pd.Series(
        #    [i.split('-')[-1] for i in adata.obs.index], dtype="category"
        #)
        # Reorder to put ctrl first as the logical first treatment duration time point
        #adata.obs.time = adata.obs.time.cat.reorder_categories(['ctrl', '22', '23', '24'])

        print('Number of genes duplicated in index:', adata.var.index.duplicated().sum())
        sc.pl.highest_expr_genes(adata, n_top=20)
        if _.output_report: report.lastSection.figs['Highest expressing genes'] = plt.gcf()
        
        # Quality control
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(9,3))
        adata.obs.n_genes_by_counts.hist(bins='auto', ax=axes[0])
        axes[0].set_title('Number of genes by counts')
        adata.obs.total_counts.hist(bins='auto', ax=axes[1])
        axes[1].set_title('Total counts')
        adata.obs.pct_counts_mt.hist(bins='auto', ax=axes[2])
        axes[2].set_title('% mt counts')
        #sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        #             jitter=0.4, multi_panel=True)
        if _.output_report: report.lastSection.figs['QC metrics'] = fig

        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color='time')#, save=True)
        if _.output_report: report.lastSection.figs['Total counts vs % mt'] = plt.gcf()
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='time')#, save=True)
        if _.output_report: report.lastSection.figs['Total counts vs number of genes'] = plt.gcf()
        
        # Filtering
        print(f'Filtering minimum gene count ({_.filter_cells}), cell count ({_.filter_genes})')
        print(f'Filtering cells with more than {_.filter_mt}% mitochondrial reads')
        print(adata) # before
        sc.pp.filter_cells(adata, min_genes=_.filter_cells)
        sc.pp.filter_genes(adata, min_cells=_.filter_genes)
        adata = adata[adata.obs.pct_counts_mt < _.filter_mt, :].copy()
        print(adata) # after filtering
        
        # Normalize/transform data
        print('Normalizing and transforming data')
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Select highly variable genes
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata.raw = adata
        #adata = adata[:, adata.var.highly_variable]
        
        #sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        #sc.pp.scale(adata, max_value=10)
        
        #sc.tl.pca(adata, svd_solver='arpack')
        #sc.pl.pca(adata, color='BRCA1')
        #sc.pl.pca_variance_ratio(adata, log=True)
        #sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
        #sc.tl.umap(adata)
        #sc.pl.umap(adata, color=['CST3', 'BRCA1', 'MALAT1'])
        #sc.tl.leiden(adata)
        #sc.pl.umap(adata, color=['leiden', 'MALAT1'])
        #sc.tl.rank_genes_groups(adata, 'leiden', method='t-test') # method='wilcoxon'
        #sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
        #pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
        return adata
        
    # Separate analysis per treatment duration
    def adata_tps(_, adata, genesets, report, print) -> dict:
        adata_tps = {}
        for tp in ('ctrl', '22', '23', '24'):
            adata_tp = adata[adata.obs.time == tp].copy()
            sc.tl.pca(adata_tp, n_comps=25)
            sc.pp.neighbors(adata_tp, n_neighbors=10)
            sc.tl.umap(adata_tp, n_components=3)
            sc.tl.leiden(adata_tp)
            sc.tl.score_genes_cell_cycle(adata_tp, genesets['S'], genesets['G2M'], use_raw=False)
            adata_tps[tp] = adata_tp
        return adata_tps

    def separate_analysis_per_treatmen_duration(_, adata_tps, print) -> TaskSection:
        figs = OrderedDict()
        tabs = OrderedDict()
        for tp in ('ctrl', '22', '23', '24'):
            adata_tp = adata_tps[tp]
            #sc.pl.umap(adata_tp, color=['CST3', 'BRCA1', 'MALAT1'])
            sc.pl.umap(adata_tp, color=['leiden','phase','CDKN1A'])
            fig = plt.gcf()
            fig.axes[1].set_ylabel('')
            fig.axes[2].set_ylabel('')
            fig.suptitle(f'Treatment {tp}')
            figs[f'Treatment {tp}'] = fig
            
        cellphases = pd.DataFrame(
            {tp:(adata_tps[tp].obs.phase.value_counts()*100/len(adata_tps[tp])).round(2)
                                   for tp in adata_tps}
        )
        cellphases.index = cellphases.index.astype(str)
        tabs['Cell phases per treatment duration'] = cellphases.copy()
        print('''There is a clear shift to G1 for any of the treatment durations,
        apparently increasing with every extra hour of treatment.
        ''')
        
        return {
            'figures': figs,
            'tables': tabs,
            'clearpage': True   
        }

    def adata_GO_senescence_analysis(_, adata, msigdb, genesets) -> SKIPCACHE(TaskSection):
        # TODO preparing umap in other section
        sc.tl.pca(adata, n_comps=25)
        sc.pp.neighbors(adata, n_neighbors=10)
        sc.tl.umap(adata, n_components=3)
        sc.tl.score_genes_cell_cycle(adata, genesets['S'], genesets['G2M'], use_raw=False)
        
        # Make ranksums for senescence gene sets
        figs = []
        for goset in msigdb['C5'].keys() | msigdb['C2'].keys():
            if 'SENESCENCE' in goset:
                gogenes = msigdb['C5' if goset in msigdb['C5'] else 'C2'][goset] 
                adata.obs[goset] = adata.to_df().T.filter(gogenes, axis=0).sum()
                # (adata.to_df().T.filter(gogenes, axis=0)>0).sum()
                # adata.to_df().T.rank().filter(gogenes, axis=0).sum()/len(gogenes)
                sc.pl.umap(adata, color=[goset,'time','phase'])
                figs.append((goset, plt.gcf()))
        sc.pl.umap(adata, color=[
            'FRIDMAN_SENESCENCE_UP', 'TANG_SENESCENCE_TP53_TARGETS_UP'
            , 'GO_CELLULAR_SENESCENCE', 'time',
            'FRIDMAN_SENESCENCE_DN', 'TANG_SENESCENCE_TP53_TARGETS_DN'
            , 'GO_REPLICATIVE_SENESCENCE', 'phase'
        ])
        figs.append(('Senescence signatures', plt.gcf()))
        return figs
    
    def joint_cell_cycle_analysis(_, adata, adata_tps, msigdb): #not activated by type
        '''
        Based on the cluster cell phase content I renamed clusters in the following way:
        
        - phase component is equal or more than 90%: phase_timely, e.g. G1_timely
        
        - 1 phase is more or equal to 50% and a second phase between 20 and 40%: phase_early
        or phase_late, depending on that second phase, e.g. G1_late
        
        - other possibilities are represented by the sorted sequence of phases, with the
        most abundant first, e.g. G1SG2M
        
        Based on this naming convention, we can make some observations. The control condition
        does not have a G1_timely cluster, whereas all treatments do. Considering that the G1
        phase increased for all the treatment conditions, this in itself could be a consequence
        of having more cells to cluster and thus an increase in the chance of getting a more
        pure G1 cluster that is then considered G1_timely, in other words it could be a consequence
        of not normalizing. On the other hand, if we consider the increase in G1 is because there is
        a possible subfraction of senescent-like or quiescent cells, the G1_timely cluster that we
        identify within the treatment conditions, would be the best clusters to look further into,
        although we then don't have a similar cluster to compare to in the control condition.
        
        The reverse is true for the S-phase: only the control condition has an S-phase dominated
        cluster 'S_late'. The treatment condition do not, although it is mainly the G2M phase
        that diminished in comparison with the control and the S-phase as far as proportion of
        cells remained more constant.
        '''
        
        # Joint analysis
        sc.tl.score_genes_cell_cycle(adata, genesets['S'], genesets['G2M'], use_raw=False)
        
        # Comparing phase assignment separate or joint
        cf = pd.concat([adata_tps[tp].obs.phase
                        for tp in adata_tps])
        cf.sort_index(inplace=True)
        cfg = adata.obs.phase.copy()
        cfg.sort_index(inplace=True)
        (cf != cfg).sum()/len(cf)
        print(pd.Series([i.split('-')[-1] for i in cf[(cf != cfg)].index]).value_counts())
        
        def clustername_decisiontree(phase_distro_pc):
            # To make sure values are sorted correctly
            phase_distro_pc = phase_distro_pc.sort_values(ascending=False)
            if phase_distro_pc[phase_distro_pc >= .9].any():
                return f'{phase_distro_pc[phase_distro_pc >= .9].index[0]}_timely'
            elif phase_distro_pc[phase_distro_pc >= .5].any() and (
                    (phase_distro_pc >= .2) & (phase_distro_pc <= .4)).sum() == 1:
                phase = phase_distro_pc.index[1]
                second_phase = phase_distro_pc.index[1]
                if phase == 'G1':
                    if second_phase == 'G2M': return 'G1_early'
                    else: return 'G1_late'
                elif phase == 'S':
                    if second_phase == 'G1': return 'S_early'
                    else: return 'S_late'
                else:
                    if second_phase == 'S': return 'G2M_early'
                    else: return 'G2M_late'
            elif phase_distro_pc[phase_distro_pc < .1].any():
                phases = phase_distro_pc[phase_distro_pc >= .1].index
                if 'S' in phases and 'G1': return 'G1_S'
                elif 'S' in phases and 'G2M': return 'S_G2M'
                else: return 'G2M_G1'
            else:
                return f'{"".join(phase_distro_pc.index)}'
        
        phaseorders = [
            'G1_early', 'G1_timely', 'G1_late', 'G1_S',
            'S_early', 'S_timely', 'S_late', 'S_G2M'
            'G2M_early', 'G2M_timely', 'G2M_late', 'G2M_G1'
            #'G2MSG1', 'G2MG1S'
        ]
        
        # Filter genesets
        Sgenes = [g for g in genesets['S'] if g in adata.var.index]
        G2Mgenes = [g for g in genesets['G2M'] if g in adata.var.index]
        # GO_senescence group only has 2 genes in our data
        GOSgenes = [g for g in genesets['GO_senescence2'] if g in adata.var.index]
        sen_inducers = [
            g for g in senescence_table['Gene Name'][senescence_table['Senescence Effect'] == 'Induces']
            if g in adata.var.index
        ]
        sen_inhibitors = [
            g for g in senescence_table['Gene Name'][senescence_table['Senescence Effect'] == 'Inhibits']
            if g in adata.var.index
        ]
        sen_oncind = [
            g for g in senescence_table['Gene Name'][senescence_table['Senescence Type'] == 'Oncogene-induced']
            if g in adata.var.index
        ]
        P53genes = [g for g in msigdb['H']['HALLMARK_P53_PATHWAY'] if g in adata.var.index]
        
        # Analysing cluster group phase content
        if _.output_report: report.append('Analyzing cluster group cell cycle phase distribution', '', clearpage=True)
        gsea_cl_results = {}
        rename_clusters = {}
        for tp in adata_tps: 
            adata_tp = adata_tps[tp]
            sc.tl.rank_genes_groups(adata_tp, 'leiden', method='t-test') # method='wilcoxon'
            adata_tp_ranks = pd.DataFrame({
                group + '_' + key[:1]: adata_tp.uns['rank_genes_groups'][key][group] 
                for group in adata_tp.uns['rank_genes_groups']['names'].dtype.names for key in ['names', 'pvals']
            })
            gsea_cl_results[tp] = {}
            rename_clusters[tp] = {}
            for cltr in adata_tp.obs.leiden.dtype.categories:
                adata_tp_cltr = adata_tp[adata_tp.obs.leiden == cltr]
                phase_distro = adata_tp_cltr.obs.phase.value_counts()
                phase_distro_pc = phase_distro/phase_distro.sum()
                rename_clusters[tp][cltr] = clustername_decisiontree(phase_distro_pc)
                if _.output_report:
                    report.lastSection.tabs[f'Phase distribution treatment {tp} cluster {cltr}'] = pd.DataFrame(
                        {'ratio':phase_distro_pc}, index=phase_distro_pc.index.astype(str)
                    )
                else:
                    print(phase_distro_pc.round(2))
                    print('P-val < 0.05', (adata_tp_ranks[f'{cltr}_p']<.05).sum())
                #print('\n'.join(adata_tp_ranks[f'{cltr}_n']))
                #if os.path.exists(f'{_.gseadir}/cluster_{tp}_{cltr}'):
                #    gsea_cl_results[tp][cltr] = True
                #    continue
                try:
                    gsea_cl_results[tp][cltr] = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/cluster_{tp}_{cltr}') if os.path.exists(
                            f'{_.gseadir}/cluster_{tp}_{cltr}') else gseapy.enrichr(
                        gene_list=adata_tp_ranks[f'{cltr}_n'],
                        description='single cell cluster',
                        gene_sets=','.join(gseasets),
                        outdir=f'{_.gseadir}/cluster_{tp}_{cltr}'
                    )
                except:
                    gsea_cl_results[tp][cltr] = None
                    print('GSEA cluster failed', tp, cltr)
                    failed_gsea += 1
                if _.output_report and gsea_cl_results[tp][cltr]:
                    report.lastSection.tabs[f'GSEA for treatment {tp} cluster {cltr}'] = reformat_gs_table(gsea_cl_results[
                            tp][cltr].results.groupby('Gene_set').head(3)[['Term', 'Adjusted P-value']]
                    )
        
            adata_tps[tp].obs['leiden_cpn'] = pd.Categorical(
                adata_tps[tp].obs.leiden.map(rename_clusters[tp]),
                categories=sorted(
                    set(rename_clusters[tp].values()), key=lambda x: phaseorders.index(x) if x in phaseorders else 100
                )
            )
            sc.pl.dotplot(adata_tps[tp], Sgenes, groupby='leiden_cpn');
            if _.output_report:
                report.lastSection.figs[f'S genes for treatment {tp}'] = plt.gcf()
            sc.pl.dotplot(adata_tps[tp], G2Mgenes, groupby='leiden_cpn');
            if _.output_report:
                report.lastSection.figs[f'G2M genes for treatment {tp}'] = plt.gcf()
            sc.pl.violin(adata_tps[tp], ['RRM2', 'CDKN1A', 'ODC1'], groupby='leiden_cpn', scale='area')
            fig = plt.gcf()
            fig.suptitle(f'Treatment {tp}')
            for ax in fig.axes:
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
            #fig.savefig(os.path.join(sc.settings.figdir,f'violin_cellphases_{tp}.pdf'))
            if _.output_report:
                report.lastSection.figs[f'violin_cellphases_{tp}'] = fig
        
        # Searching for PCs linked to cell phase
        if _.output_report: report.append('PCA component cell phase association', '', clearpage=True)
        for tp in adata_tps:
            sc.pl.pca(adata_tps[tp], color=['phase'], components=['1,2','1,3','1,4','1,5'])
            if _.output_report:
                report.lastSection.figs[f'PCA overview for treatment {tp}'] = plt.gcf()

    def PC_gsea(_, adata_tps, report):
        '''
        The cluster and cell phase comparison in the previous section is based on the leiden
        clustering performed, which takes in account all the genes. In this section, we look
        at the PCA components and how they associate visually with the cell phase. For all
        conditions PC1 divides between G1 and S+G2M. PC3 for control and treatment 23 allows
        further separating S and G2M, similar to PC4 for treatment 22 and 24. Based on a minimum
        PC contributing cutoff, I selected the relevant PC genes and did enrichr. Results in the
        following section's tables.
        
        
        Given that we can align the PC's in the following way as per their associations with cell cycle
        - PC1 for all conditions//
        - PC2 for all conditions//
        - PC3 for ctrl and trt 23 and PC4 for trt 22 and 24//
        - PC4 for ctrl and trt 23 and PC3 for trt 22 and 24//
        in the following tables we look at the intersect of genes across a principal component,
        once for all conditions and once substracting genes in common with the PC for control.
        The tables show the enrichr results.
        '''
        # Analysing PCA component geneset enrichment
        if _.output_report: report.append('PCA component geneset enrichment', '', clearpage=True)
        print('Analysing PCA component geneset enrichment in treatments separately')
        gene_min_contr = 0.01 # gene minimum contribution to PC
        gsea_pc_results = {}
        for tp in adata_tps:
            adata_tp = adata_tps[tp]
            gsea_pc_results[tp] = []
            for pc in range(5): # number of principal components to analyse
                pc_geneset = list(adata_tp.var.index[(adata_tp.varm['PCs'][:,pc] >= gene_min_contr)])
                try:
                    gsea_pc_results[tp].append(Reloaded_Enrichr_Results(
                        f'{_.gseadir}/cluster_{tp}_pc_{pc}') if os.path.exists(
                            f'{_.gseadir}/cluster_{tp}_pc_{pc}') else gseapy.enrichr(
                        gene_list=pc_geneset,
                        description='single cell pc',
                        gene_sets=','.join(gseasets),
                        outdir=f'{_.gseadir}/cluster_{tp}_pc_{pc}'
                    ))
                except:
                    gsea_pc_results[tp][pc] = None
                    print('GSEA PC failed', tp, pc)
                    failed_gsea += 1
                if _.output_report and gsea_pc_results[tp][pc]:
                        report.lastSection.tabs[f'Treatment {tp} PC {pc+1}'] = reformat_gs_table(gsea_pc_results[
                            tp][pc].results.groupby('Gene_set').head(3)[['Term','Adjusted P-value']]
                        )
        
        pc_genesets = {}
        for tp in adata_tps:
            pc_genesets[tp] = {}
            for pc in range(5):
                adata_tp = adata_tps[tp]
                pc_geneset = list(adata_tp.var.index[abs(adata_tp.varm['PCs'][:,pc]) >= gene_min_contr])
                if not _.output_report:
                    print(
                        tp, pc, len(pc_geneset), len(set(pc_geneset)&set(Sgenes)),
                        len(set(pc_geneset)&set(G2Mgenes)),
                        len(set(pc_geneset)&set(P53genes)),
                        'CDKN1A' in pc_geneset
                    )
                pc_genesets[tp][pc+1] = set(pc_geneset)
        
        # Based on cell phase resolution of PCs, analysis of their intersects
        # Overview of the PC resolving cell phase visually
        # ctrl -> 1 and 3
        # 22   -> 1 and 4
        # 23   -> 1 and 3
        # 24   -> 1 and 4
        if _.output_report: report.append('Analysis of PC gene intersects', '', clearpage=True)
        
        # PC 1
        pc_1_intersect = set.intersection(*[pc_genesets[tp][1] for tp in pc_genesets])
        pc_1_intersect_trt = set.intersection(*[pc_genesets[tp][1] for tp in pc_genesets if tp != 'ctrl'])
        gsea_results = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/pc_1_intersect') if os.path.exists(
                            f'{_.gseadir}/pc_1_intersect') else gseapy.enrichr(
            gene_list=list(pc_1_intersect), description='single cell pc',
            gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/pc_1_intersect'
                           )                
        if _.output_report:
            report.lastSection.tabs['PC 1 common genes'] = reformat_gs_table(gsea_results.results.groupby('Gene_set').head(3)[
                ['Term', 'Adjusted P-value']]
            )
        gsea_results = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/pc_1_trt_intersect_unique') if os.path.exists(
                            f'{_.gseadir}/pc_1_trt_intersect_unique') else gseapy.enrichr(
            gene_list=list(pc_1_intersect_trt-pc_1_intersect), description='single cell pc',
            gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/pc_1_trt_intersect_unique'
                           )                
        if _.output_report:
            report.lastSection.tabs['PC 1 treatment unique genes'] = reformat_gs_table(
                gsea_results.results.groupby('Gene_set').head(3)[
                ['Term', 'Adjusted P-value']]
            )
        
        # PC 2
        pc_2_intersect = set.intersection(*[pc_genesets[tp][2] for tp in pc_genesets])
        pc_2_intersect_trt = set.intersection(*[pc_genesets[tp][2] for tp in pc_genesets if tp != 'ctrl'])
        gsea_results = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/pc_2_intersect') if os.path.exists(
                            f'{_.gseadir}/pc_2_intersect') else gseapy.enrichr(
            gene_list=list(pc_2_intersect), description='single cell pc',
            gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/pc_2_intersect'
                           )                
        if _.output_report:
            report.lastSection.tabs['PC 2 common genes']: reformat_gs_table(
                gsea_results.results.groupby('Gene_set').head(3)[['Term', 'Adjusted P-value']]
            )
        gsea_results = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/pc_2_trt_intersect_unique') if os.path.exists(
                            f'{_.gseadir}/pc_2_trt_intersect_unique') else gseapy.enrichr(
            gene_list=list(pc_2_intersect_trt-pc_2_intersect), description='single cell pc',
            gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/pc_2_trt_intersect_unique'
                           )                
        if _.output_report:
            report.lastSection.tabs['PC 2 treatment unique genes'] = reformat_gs_table(
                gsea_results.results.groupby('Gene_set').head(3)[['Term', 'Adjusted P-value']]
            )
        
        # PC 3 or 4 depending on visual analysis (choosing the one that best visualizes cell cycle)
        pc_SvsG2M_intersect = set.intersection(*[pc_genesets[tp][3 if tp in ('ctrl','23') else 4] for tp in pc_genesets])
        pc_SvsG2M_intersect_trt = set.intersection(
            *[pc_genesets[tp][3 if tp in ('ctrl','23') else 4] for tp in pc_genesets if tp != 'ctrl']
        )
        gsea_results = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/pc_SvsG2M_intersect') if os.path.exists(
                            f'{_.gseadir}/pc_SvsG2M_intersect') else gseapy.enrichr(
            gene_list=list(pc_SvsG2M_intersect), description='single cell pc',
            gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/pc_SvsG2M_intersect'
                           )                
        if _.output_report:
            report.lastSection.tabs['S vs G2M PC common genes'] = reformat_gs_table(
                gsea_results.results.groupby('Gene_set').head(3)[['Term', 'Adjusted P-value']]
            )
        gsea_results = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/pc_SvsG2M_trt_intersect_unique') if os.path.exists(
                            f'{_.gseadir}/pc_SvsG2M_trt_intersect_unique') else gseapy.enrichr(
            gene_list=list(pc_SvsG2M_intersect_trt-pc_SvsG2M_intersect), description='single cell pc',
            gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/pc_SvsG2M_trt_intersect_unique'
                           )                
        if _.output_report:
            report.lastSection.tabs['S vs G2M PC treatment unique genes'] = reformat_gs_table(
                gsea_results.results.groupby('Gene_set').head(3)[['Term', 'Adjusted P-value']]
            )
        
        # pc_SvsG2M complement
        pc_compl_intersect = set.intersection(*[pc_genesets[tp][4 if tp in ('ctrl','23') else 3] for tp in pc_genesets])
        pc_compl_intersect_trt = set.intersection(
            *[pc_genesets[tp][4 if tp in ('ctrl','23') else 3] for tp in pc_genesets if tp != 'ctrl']
        )
        gsea_results = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/pc_compl_intersect') if os.path.exists(
                            f'{_.gseadir}/pc_compl_intersect') else gseapy.enrichr(
            gene_list=list(pc_compl_intersect), description='single cell pc',
            gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/pc_compl_intersect'
                           )
        if _.output_report:
            report.lastSection.tabs['4th or 3rd PC common genes'] = reformat_gs_table(
                gsea_results.results.groupby('Gene_set').head(3)[['Term', 'Adjusted P-value']]
            )
        gsea_results = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/pc_compl_trt_intersect_unique') if os.path.exists(
                            f'{_.gseadir}/pc_compl_trt_intersect_unique') else gseapy.enrichr(
            gene_list=list(pc_compl_intersect_trt-pc_compl_intersect), description='single cell pc',
            gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/pc_compl_trt_intersect_unique'
                           )
        if _.output_report:
            report.lastSection.tabs['4th or 3rd PC treatment unique genes'] = reformat_gs_table(
                gsea_results.results.groupby('Gene_set').head(3)[['Term', 'Adjusted P-value']]
            )
        #.first

        return pc_genesets

    def bulk_sigs_in_tp_clusters(_, adata_tps, diffgenesbulk, print) -> TaskSection:
        from scipy.stats import spearmanr#, kstest, normaltest
        figs = OrderedDict()
        for direction in diffgenesbulk:
            print('Bulk signature logFC direction', direction)
            ranksum_means = {}
            ranksum_vars = {}
            cellphases = {} # could inject from other subtask
            for tp in adata_tps:
                ranksum_means[tp] = {}
                ranksum_vars[tp] = {}
                cellphases[tp] = {}
                adata_tp = adata_tps[tp]
                for cltr in adata_tp.obs.leiden.dtype.categories:
                    adata_tp_cltr = adata_tp[adata_tp.obs.leiden == cltr]
                    phase_distro = adata_tp_cltr.obs.phase.value_counts()
                    phase_distro_pc = phase_distro/phase_distro.sum()
                    cellphases[tp][int(cltr)] = phase_distro_pc
                    ranksums = {
                        gs: adata_tp_cltr.to_df().T.rank().filter(
                            diffgenesbulk[direction][gs], axis=0
                        ).sum()
                        for gs in diffgenesbulk[direction]
                    }
                    #kstest(ranksums, 'norm', args=(ranksums.mean(),ranksums.std())).pvalue < .01
                    ranksum_means[tp][int(cltr)] = {gs: ranksums[gs].mean() for gs in ranksums}
                    ranksum_vars[tp][int(cltr)] = {gs: ranksums[gs].var() for gs in ranksums}
    
            figs[f'Bulk {direction} signature'], axes = plt.subplots(
                nrows=2, ncols=2, figsize=(14,7), sharex=True, sharey=True
            )
            axes[0,0].set_ylabel('higher rank ~ higher expression')
            axes[1,0].set_ylabel('Mean gene rank')
            axes[1,0].set_xlabel('G1 content (ratio)')
            axes[1,1].set_xlabel('G1 content (ratio)')
            for tp, ax in zip(ranksum_means, axes.flatten()):
                df = pd.DataFrame(ranksum_means[tp]).T
                cdf = pd.DataFrame(cellphases[tp]).fillna(0)
                correlations = df.apply(lambda x: spearmanr(cdf.loc['G1'], x)).T
                print(tp, correlations, '', sep='\n')
                for col in df:
                    #ax.scatter(df.index, df[col], s=cdf.T.G1 * 25 + 5, label=col)
                    ax.scatter(cdf.T.G1, df[col]/len(diffgenesbulk[direction][col]), label=col.replace('RPEMYCN_',''))
                ax.set_title(f'Treatment {tp}')
            #ax.legend(loc='lower right')
            axes[0,0].legend()
        return {
            'title': 'Bulk signatures in treatment duration clusters',
            'clearpage': True,
            'figures': figs
            }
                
    def normality_tests_of_bulk_signatures(_, adata_tps, diffgenesbulk, print) -> TaskSection:
        from scipy.stats import kstest
        from statsmodels.stats.multitest import multipletests
        ktest_results = []
        figs = OrderedDict()
        for tp in adata_tps:
            adata_tp = adata_tps[tp]
            for direction in diffgenesbulk:
                for gs in diffgenesbulk[direction]:
                    _gs = gs.replace('RPEMYCN_','')
                    ranksums = adata_tp.to_df().T.rank().filter(
                        diffgenesbulk[direction][gs], axis=0
                    ).sum()/len(diffgenesbulk[direction][gs])
                    adata_tp.obs[f'{_gs}_{direction}_score'] = ranksums
                    kstest_ranks = kstest(ranksums, 'norm', args=(ranksums.mean(),ranksums.std()))
                    ktest_results.append({
                        'tp':tp,
                        'direction':direction,
                        'sig':gs,
                        'pval': kstest_ranks.pvalue}
                    )
                    if kstest_ranks.pvalue < .01:
                        print(f'{tp} - {direction} - {gs}: {kstest_ranks}')
                        #fig, ax = plt.subplots()
                        #ranksums.hist(bins=20, ax=ax)
                        #ax.set_title(f'{tp} - {direction} - {gs}')
            print()
            for gs in diffgenesbulk[direction]:
                _gs = gs.replace('RPEMYCN_','')
                sc.pl.umap(adata_tp, color=['phase',f'{_gs}_up_score',f'{_gs}_down_score','CDKN1A'])
                fig = plt.gcf()
                figs[f'{_gs} signature in {tp}'] = fig

        ktest_results = pd.DataFrame(ktest_results)
        ktest_results['fdr'] = multipletests(
            pvals=ktest_results.pval, alpha=0.05, method="fdr_bh"
        )[1]
        ktest_results.sort_values('fdr', inplace=True)
        return {
            'figures': figs,
            'tables':{'ktest results':ktest_results},
            'clearpage': True
        }
            
    def cell_cycle_trajectories(_, adata_tps, report) -> SKIP(type(None)):
        # Cell cycle trajectories
        # Trajectory
        for tp in adata_tps:
            print(tp)
            adata_tps[tp].obs['G1_content'] = adata_tps[tp].obs.phase == 'G1'
            adata_tps[tp].obs['S_content'] = adata_tps[tp].obs.phase == 'S'
            adata_tps[tp].obs['G2M_content'] = adata_tps[tp].obs.phase == 'G2M'
            sc.tl.pca(adata_tps[tp], svd_solver='arpack')
            sc.pp.neighbors(adata_tps[tp], n_neighbors=4, n_pcs=10) # n_pcs=20
            # minimum of 4 neighbors seems necessary to do diffmap later
            sc.tl.draw_graph(adata_tps[tp])
            sc.pl.draw_graph(adata_tps[tp], color='leiden', legend_loc='on data')
            sc.tl.diffmap(adata_tps[tp])
            sc.pp.neighbors(adata_tps[tp], n_neighbors=_.traject_nbr, use_rep='X_diffmap')
            sc.tl.draw_graph(adata_tps[tp])
            sc.pl.draw_graph(adata_tps[tp], color='leiden', legend_loc='on data')
            sc.tl.paga(adata_tps[tp], groups='leiden')
            sc.pl.paga(adata_tps[tp], color=['leiden', 'G1_content', 'S_content','G2M_content'])
        
        # Combined analysis of G1 timely cells across conditions
        if _.output_report:
            report.append('Combined analysis of G1 timely cells', clearpage=True)
        # renaming G1_S in ctrl group to G1_timely -> its composition (G1 85% S 11% G2M 4%) was just outside
        # my 'timely' threshold
        adata_tps['ctrl'].obs.leiden_cpn.cat.rename_categories({'G1_S' : 'G1_timely'}, inplace=True)
        # combine annotation
        adata.obs['separate_leiden_cpn'] = pd.concat(
            [adata_tps[tp].obs.leiden_cpn for tp in adata_tps]
        )
        adata.obs['separate_phase'] = pd.concat(
            [adata_tps[tp].obs.phase for tp in adata_tps]
        )
        # filter for G1_timely and phase G1 cells
        adata_G1 = adata[
            (adata.obs.separate_leiden_cpn == 'G1_timely') & (adata.obs.separate_phase == 'G1'),:
        ].copy()
        print(adata_G1.obs.time.value_counts())
        sc.tl.pca(adata_G1, n_comps=25)
        sc.pp.neighbors(adata_G1, n_neighbors=10)
        sc.tl.umap(adata_G1, n_components=3)
        sc.tl.leiden(adata_G1)
        sc.tl.score_genes_cell_cycle(adata_G1, genesets['S'], genesets['G2M'], use_raw=False)
        sc.pl.umap(adata_G1, color=['leiden','time','CDKN1A']) # optionally also plot phase
        sc.tl.rank_genes_groups(adata_G1, 'leiden', method='wilcoxon')
        sc.pl.rank_genes_groups(adata_G1, n_genes=25, sharey=False)
        #ranked_G1_clustergenes = pd.DataFrame(adata_G1.uns['rank_genes_groups']['names'])
        #ranked_G1_clustergenes_gsea = {}
        #for cltr in ranked_G1_clustergenes:
        #    ranked_G1_clustergenes_gsea[cltr] = Reloaded_Enrichr_Results(
        #                f'{_.gseadir}/G1_cluster_{cltr}') if os.path.exists(
        #                    f'{_.gseadir}/G1_cluster_{cltr}') else gseapy.enrichr( 
        #        gene_list=ranked_G1_clustergenes[cltr], description='single cell G1 cluster', 
        #        gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/G1_cluster_{cltr}' 
        #    )
        
        # PC components
        sc.pl.pca(adata_G1, color=['time'], components=['1,2','1,3','1,4','1,5'])
        adata_G1.obs['total_counts_qcut'] = pd.qcut(adata_G1.obs.total_counts,5)
        sc.pl.pca(adata_G1, color=['total_counts_qcut'], components=['1,2','1,3','1,4','1,5'])
        pc_geneset = list(adata_G1.var.index[(adata_G1.varm['PCs'][:,0] >= gene_min_contr)])
        
        # Trajectory
        sc.tl.pca(adata_G1, svd_solver='arpack')
        sc.pp.neighbors(adata_G1, n_neighbors=4, n_pcs=20)
        # minimum of 4 neighbors seems necessary to do diffmap later
        sc.tl.draw_graph(adata_G1) # pip install fa2 failed
        sc.pl.draw_graph(adata_G1, color='leiden', legend_loc='on data')
        sc.tl.diffmap(adata_G1)
        sc.pp.neighbors(adata_G1, n_neighbors=_.traject_nbr, use_rep='X_diffmap')
        sc.tl.draw_graph(adata_G1)
        sc.pl.draw_graph(adata_G1, color='leiden', legend_loc='on data')
        sc.tl.paga(adata_G1, groups='leiden')
        sc.pl.paga(adata_G1, color=['leiden', 'YBX1', 'CDKN1A'])
        
        #sc.pl.violin(adata_G1, ['YBX1', 'CDKN1A','ODC1'], groupby='leiden', scale='area')
        sc.pl.heatmap(adata_G1, ['YBX1', 'CDKN1A','TP53','ODC1'], groupby='leiden')
        
        # Unexplored strategies from Vladie
        #print('DPT analysis')
        #sc.tl.diffmap(adata,n_comps=15)
        #adata.uns['iroot'] = np.flatnonzero(adata.obs['time'] == 'ctrl')[0]
        #sc.tl.dpt(adata, n_branchings=1)
        
        # sc.tl.dendrogram(adata2, groupby='louvain', n_pcs=15,cor_method='pearson')
        # sc.tl.rank_genes_groups(adata2, groupby='louvain', method='wilcoxon',n_genes=50)
        # sc.pl.rank_genes_groups_heatmap(adata2, n_genes=20, standard_scale='var',save='heatmap20.pdf')
        # sc.pl.rank_genes_groups_matrixplot(adata2, n_genes=50, standard_scale='var', cmap='Blues',save='matrix_plot30.pdf')
        # sc.pl.diffmap(adata, color=senescence_signatures, wspace=0.4, size=10, color_map='Reds',save='dpt_overview.pdf')
        # sc.pl.dpt_groups_pseudotime(adata)


    def q_cut_total_counts(_) -> SKIP(type(None)):
        # Exploring dividing up in total_counts qcut across phases
        totcountqcut = lambda x: pd.qcut(x.total_counts, 4, labels=False)
        adata_tps[tp].obs['total_counts_phaq'] = adata_tps[tp].obs[
               ['phase','total_counts']
        ].groupby('phase').apply(totcountqcut).reset_index().set_index('index').total_counts
        adata_tq = adata_tps[tp][adata_tps[tp].obs.total_counts_phaq == 1,:].copy()
        sc.tl.pca(adata_tq, n_comps=25)
        sc.pp.neighbors(adata_tq, n_neighbors=10)
        sc.tl.umap(adata_tq, n_components=3)
        sc.tl.leiden(adata_tq)
        sc.tl.score_genes_cell_cycle(adata_tq, genesets['S'], genesets['G2M'], use_raw=False)
        sc.pl.umap(adata_tq, color=['leiden','phase','CDKN1A']) # optionally also plot phase
        sc.tl.rank_genes_groups(adata_tq, 'leiden', method='wilcoxon')
        sc.pl.rank_genes_groups(adata_tq, n_genes=25, sharey=False)
        #ranked_tq_clustergenes = pd.DataFrame(adata_tq.uns['rank_genes_groups']['names'])
        #ranked_tq_clustergenes_gsea = {}
        #for cltr in ranked_tq_clustergenes:
        #    ranked_tq_clustergenes_gsea[cltr] = Reloaded_Enrichr_Results(
        #                f'{_.gseadir}/G1_cluster_{cltr}') if os.path.exists(
        #                    f'{_.gseadir}/G1_cluster_{cltr}') else gseapy.enrichr( 
        #        gene_list=ranked_tq_clustergenes[cltr], description='single cell G1 cluster', 
        #        gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/G1_cluster_{cltr}' 
        #    )
        
        # PC components
        sc.pl.pca(adata_tq, color=['total_counts','phase'], components=['1,2','1,3','1,4','1,5'])
        pc_geneset = list(adata_tq.var.index[(adata_tq.varm['PCs'][:,0] >= gene_min_contr)])
        
        # Trajectory
        sc.tl.pca(adata_tq, svd_solver='arpack')
        sc.pp.neighbors(adata_tq, n_neighbors=4, n_pcs=20)
        # minimum of 4 neighbors seems necessary to do diffmap later
        sc.tl.draw_graph(adata_tq) # pip install fa2 failed
        sc.pl.draw_graph(adata_tq, color='leiden', legend_loc='on data')
        sc.tl.diffmap(adata_tq)
        sc.pp.neighbors(adata_tq, n_neighbors=_.traject_nbr, use_rep='X_diffmap')
        sc.tl.draw_graph(adata_tq)
        sc.pl.draw_graph(adata_tq, color='leiden', legend_loc='on data')
        sc.tl.paga(adata_tq, groups='leiden')
        sc.pl.paga(adata_tq, color=['leiden', 'G1_content', 'S_content','G2M_content'])

    def overlap_analysis(_) -> SKIP(type(None)):
        # Comparison ctrl vs 22, 23 and 24 separately
        diff_with_ctrl = {}
        diff_for_ctrl = {}
        for tp in ('22', '23', '24'):
            adata_subset = adata[adata.obs.time.isin(('ctrl',tp)), :].copy()
            sc.tl.rank_genes_groups(adata_subset, groupby='time', method='wilcoxon', n_genes=1000)
            adata_subset_ranks = pd.DataFrame({ 
                group + '_' + key[:1]: adata_subset.uns['rank_genes_groups'][key][group]  
                for group in adata_subset.uns['rank_genes_groups']['names'].dtype.names for key in ['names', 'pvals'] 
            })
            diff_with_ctrl[tp] = adata_subset_ranks
            diff_for_ctrl[tp] = set(
                diff_with_ctrl[tp][diff_with_ctrl[tp].ctrl_p < .01].ctrl_n
            )
                
        print(
            'Overlap 22 vs 23:',
            len(set(diff_with_ctrl['22']['22_n'])&set(diff_with_ctrl['23']['23_n'])),
            '\nOverlap 22 vs 24:',
            len(set(diff_with_ctrl['22']['22_n'])&set(diff_with_ctrl['24']['24_n'])),
            '\nOverlap 23 vs 24:',
            len(set(diff_with_ctrl['23']['23_n'])&set(diff_with_ctrl['24']['24_n'])),
            '\nOverlap ctrl (22) vs ctrl (23):',
            len(set(diff_with_ctrl['22']['ctrl_n'])&set(diff_with_ctrl['23']['ctrl_n'])),
            '\nOverlap ctrl (22) vs ctrl (23):',
            len(set(diff_with_ctrl['22']['ctrl_n'])&set(diff_with_ctrl['24']['ctrl_n'])),
            '\nOverlap ctrl (22) vs ctrl (24):',
            len(set(diff_with_ctrl['23']['ctrl_n'])&set(diff_with_ctrl['24']['ctrl_n'])),
            '\nOverlap 22 vs ctrl (23):',
            len(set(diff_with_ctrl['22']['22_n'])&set(diff_with_ctrl['23']['ctrl_n']))
        )
            
    # Differential expression test with diffxpy
    def diffxpy_test(_, adata, gseasets, report, print) -> list:
        # Manually load cache -> in the future this should be automated, as should be adding it to
        # report section
        #if os.path.exists('cvn_analysis.pickle'):
        #    previous_output = pickle.load(open('cvn_analysis.pickle', 'rb'))['diffxpy_test']
        #else:
        previous_output = False
        
        
        if _.output_report: report.append('Diffxpy non-quantitative analysis', '', clearpage=True)
        print('starting de test all data')
        if previous_output:
            test_summary = previous_output[1]
        else:
            test = de.test.wald(data=adata, formula_loc='~1 + time', factor_loc_totest='time')
            #test.plot_volcano(corrected_pval=True, min_fc=1.05, alpha=0.05, size=20)
            test_summary = test.summary()
            test_summary['pi'] = (test_summary.pval.apply(np.log10)*-1*test_summary.log2fc).fillna(0)
            #test_summary.sort_values('pi', inplace=True) # volgens Jo+Pieter issues met pi
            test_summary.sort_values('log2fc', inplace=True)
            test_summary.set_index('gene', inplace=True)

        for gseaset in gseasets:
            pre_res = gseapy.prerank(
                rnk=test_summary.log2fc.dropna().rank(),#log2fc,#.pi,
                gene_sets=gseaset,
                processes=4,
                permutation_num=_.gsea_perms,
                outdir=f'/tmp/prerank_report_{gseaset}_fc_all', format='png', seed=6
            )
            if _.output_report:
                report.lastSection.tabs[f'Prerank all GSEA {gseaset}'] = reformat_prerank_table(
                    pre_res.res2d.head(5)[['es', 'pval', 'fdr']]
                )
            else:
                print(gseaset)
                print(pre_res.res2d.sort_values('pval').head())
                print()
                
            for term in pre_res.res2d.index[:3]:
                gseaplot(pre_res.ranking, term=term, **pre_res.results[term])
                if _.output_report:
                    report.lastSection.figs[f'Prerank all GSEA {gseaset}: {term}'] = plt.gcf()
        
        # Regressing Timepoint
        if _.output_report: report.append('Diffxpy time-quantitative analysis', '', clearpage=True)
        print('starting de test regression')
        adata_tps = adata[adata.obs.time.isin({'22','23','24'}),:].copy()
        adata_tps.obs.time = adata_tps.obs.time.astype(int)
        if previous_output:
            test_regressed_treatdur_summary = previous_output[3]
        else:
            test_regressed_treatdur = de.test.wald(
                data=adata_tps,
                formula_loc="~ 1 + time",
                factor_loc_totest="time",
                as_numeric=["time"]
            )
            #test_regressed_treatdur.plot_volcano(
            #    corrected_pval=True, min_fc=1.05, alpha=0.05, size=20
            #)
            test_regressed_treatdur_summary = test_regressed_treatdur.summary()
            test_regressed_treatdur_summary['pi'] = (
                test_regressed_treatdur_summary.pval.apply(np.log10)*-1* test_regressed_treatdur_summary.log2fc
            ).fillna(0).rank()
            #test_regressed_treatdur_summary.sort_values('pi', inplace=True)
            test_regressed_treatdur_summary.sort_values('log2fc', inplace=True)
            test_regressed_treatdur_summary.set_index('gene', inplace=True)
        for gseaset in gseasets:
            pre_res_regr = gseapy.prerank(
                rnk=test_regressed_treatdur_summary.log2fc.dropna().rank(),#log2fc,#.pi,
                gene_sets=gseaset,
                processes=4,
                permutation_num=_.gsea_perms,
                outdir=f'/tmp/prerank_report_{gseaset}_fc_quantitative', format='png', seed=6
            )
            if _.output_report:
                report.lastSection.tabs[f'Prerank GSEA {gseaset}'] = reformat_prerank_table(
                    pre_res_regr.res2d.head(5)[['es', 'pval', 'fdr']]
                )
            else:
                print(gseaset)
                print(pre_res_regr.res2d.sort_values('pval').head())
                print()
                
            for term in pre_res_regr.res2d.index[:3]:
                gseaplot(pre_res_regr.ranking, term=term, **pre_res_regr.results[term])
                if _.output_report:
                    report.lastSection.figs[f'Prerank GSEA {gseaset}: {term}'] = plt.gcf()
            # genes = pre_res_regr.res2d.genes[0].split(";")
            # need to use scanpy heatmap for preranked
            #heatmap(df = pre_res_regr.heatmat.loc[genes], z_score=0, title=pre_res_regr.res2d.index[0], figsize=(18,6))
            
        trt_sig = test_regressed_treatdur_summary.loc[test_regressed_treatdur_summary.qval<.05]
        sc.pl.violin(adata, ['SNHG25', 'TSPAN3','AP005233.2'], groupby='time', scale='area')
        sc.pl.violin(adata, ['AC010883.1', 'NATD1','AP005233.2'], groupby='time', scale='area')
        adata_sig = adata[:,adata.var.index.isin(trt_sig.index)].copy()
        
        dendrogram_sig = sc.tl.dendrogram(adata, groupby='time', var_names=trt_sig.index)
        sc.pl.heatmap(adata, trt_sig[trt_sig.qval <= .001].sort_values('log2fc').index, groupby='time', dendrogram=False)
        sc.pl.dotplot(adata, trt_sig[trt_sig.qval <= .001].sort_values('log2fc').index, groupby='time')
        fig = plt.gcf()
        for ax in fig.axes:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

        # Geneset enrichment with enrichr
        gsea_results = Reloaded_Enrichr_Results(
                        f'{_.gseadir}/diffxpy_quant') if os.path.exists(
                            f'{_.gseadir}/diffxpy_quant') else gseapy.enrichr( 
                gene_list=list(trt_sig.index), description='single cell diffxpy', 
                gene_sets=','.join(gseasets), outdir=f'{_.gseadir}/diffxpy_quant' 
        )
        if _.output_report:
            report.lastSection.tabs['Diffxpy quantitative analysis'] = reformat_gs_table(
                gsea_results.results.groupby('Gene_set').head(5)#[['Term', 'Adjusted P-value']]
            )
            
        return [pre_res, test_summary, pre_res_regr, test_regressed_treatdur_summary]

    def bulk_signature_analysis(_, adata, diffgenesbulk, diffxpy_test, bulk_de_thresh, print) -> type(None):
        from scipy.stats import fisher_exact
        #from scipy.stats import hypergeom
        diffgenesbulk = diffgenesbulk['all']
        intersection = diffgenesbulk['intersection']
        intersection_ix = adata.var.index[adata.var.index.isin(intersection)]
        sc.pl.heatmap(adata, intersection_ix, groupby='time')
        sc.pl.dotplot(adata, intersection_ix, groupby='time')
        utbulk_vs24h_ix = adata.var.index[
            adata.var.index.isin(
                diffgenesbulk['RPEMYCN_24h_Vs_RPEMYCN_UT'].head(50))
        ]
        sc.pl.dotplot(adata, utbulk_vs24h_ix, groupby='time')

        qsc=diffxpy_test[3] # quantitative regression results
        bts=diffgenesbulk
        # TODO manually set selecting logfc direction
        sigset = set(qsc.loc[(qsc.qval <= bulk_de_thresh)].index)#&(qsc.log2fc < 0)].index)
        # previously set to .05 but should be same as bulk threshold for fair comparison

        # Alternative way to do it would be to use fisher.test in R :
        #list1=3000
        #list2=400
        #overlap=100
        #popsize=15000
        #a=popsize-list1-list2+overlap
        #b=list2-overlap
        #c=list1-overlap
        #d=overlap
        #fisher.test(matrix(c(a,b,c,d),nrow=2,byrow=T),alternative="greater")
        #reprfact=overlap/((list1*list2)/popsize)
        #reprfact # representation factor
        #OR=(a*d)/(b*c)
        #OR # odds ratio

        for contrast in bts:
            overlap = len(sigset&set(bts[contrast]))
            set1 = len(sigset)# - set(bts[contrast]))
            set2 = len((set(qsc.index)&set(bts[contrast])))#-sigset)
            popsize = len(qsc)
            a = popsize-set1-set2+overlap
            b = set2-overlap
            c = set1-overlap
            d = overlap
            print(
                contrast,
                overlap, set1, popsize, set2,
                fisher_exact([[a, b], [c, d]]) # alternative='greater'
                #hypergeom.cdf(overlap-1, set1, popsize-set1, set2)
            )
            
    def bulk_cluster_signatures(_, adata, print) -> SKIPCACHE(TaskSection):
        import zipfile, io
        figs = OrderedDict()
        zipfolder = zipfile.ZipFile('../RPE_MYCN_timeseries/bulk_analysis_report_cvn.zip')
        for i in range(10):
            clustername = f'Cluster_{i}'
            clusterset = pd.read_csv(
                io.TextIOWrapper(
                    zipfolder.open(
                        f's2_Cluster_center_correlations/table{1+i}_{clustername}.csv'
                    )
                ), sep=';', decimal=','
            ).dropna()
            print(clustername, (clusterset['0'] >= .8).sum())
            clusterset = set(clusterset.gene_id[clusterset['0'] >= .8])
            adata.obs[clustername] = adata.to_df().T.filter(clusterset, axis=0).sum()
            sc.pl.umap(adata, color=[clustername,'phase','time'])
            figs[clustername] = plt.gcf()
        return {'figures':figs,'clearpage':True}
    
if __name__ == '__main__':
    plt.ion() # important other wise can hang because of plt.show of scanpy
    analysis = Analysis(parse=True)
    analysis.run()#load_cache=True)
    analysis.report.outputZip('svg')
    plt.close('all')


# %run -i RPE_MYCN_all_timepoints.py
# plt.close('all')
# analysis._output.keys()
# adata = analysis._output['adata']
# adata.obs.columns
# adata.obs.phase
# adata.obs.time
# sc.pl.scatter(adata, x='CDKN1A', y='P53', color='time')
# adata.var.index
# 'CDKN1A' in adata.var.index
# 'P53' in adata.var.index
# 'TP53' in adata.var.index
# sc.pl.scatter(adata, x='CDKN1A', y='TP53', color='time')
# sc.pl.scatter(adata, x='CDKN1A', y='TP53', color=['time','phase'])
# sc.pl.scatter(adata, x='CDKN1A', y='TP53', color='phase')
# sc.pl.scatter?
# sc.pl.scatter(adata, x='CDKN1A', y='TP53', color='phase', use_raw=True)
# adata.raw
# sc.pl.scatter(adata.raw, x='CDKN1A', y='TP53', color='phase')
# adata_raw = sc.read(analysis.data)
# adata_raw.raw
# adata_raw.obs['phase'] = adata.obs.phase
# adata_raw.obs['time'] = adata.obs.time
# sc.pl.scatter(adata_raw, x='CDKN1A', y='TP53', color='phase')
# sc.pl.scatter(adata_raw, x='CDKN1A', y='TP53', color='time')
# sc.pl.scatter(adata_raw, x='CDKN1A', y='TP53', color='time', s=20)
# sc.pl.scatter(adata_raw, x='CDKN1A', y='TP53', color='time', size=20)
# sc.pl.scatter(adata_raw, x='CDKN1A', y='TP53', color='time', size=50)
# sc.pl.scatter(adata_raw, x='CDKN1A', y='TP53', color='time', size=100)
# sc.pl.scatter(adata_raw, x='CDKN1A', y='TP53', color='time', size=100, alpha=.1)
# sc.pl.scatter(adata_raw, x='CDKN1A', y='TP53', color='time', size=100, alpha=.3)
# sc.pl.scatter(adata_raw, x='CDKN1A', y='TP53', color='phase', size=100, alpha=.3)
# sc.pl.scatter?
# sc.pl.scatter(adata_raw, ['CDKN1A', 'TP53'], groups='phase')
# sc.pl.heatmap(adata_raw, ['CDKN1A', 'TP53'], groups='phase')
# sc.pl.heatmap(adata_raw, ['CDKN1A', 'TP53'], groupby='phase')
# sc.pl.heatmap(adata_raw, ['CDKN1A', 'TP53'], groupby='time')
# sc.pl.heatmap(adata, ['CDKN1A', 'TP53'], groupby='time')
# adata.obs['tc'] = adata.obs.time+adata.obs.phase
# adata.obs['tc'] = adata.obs.time.astype(str)+adata.obs.phase.astype(str)
# sc.pl.heatmap(adata, ['CDKN1A', 'TP53'], groupby='tc')
# sc.pl.dotplot(adata, ['CDKN1A', 'TP53'], groupby='tc')
# adata_raw.obs['tc'] = adata_raw.obs.time.astype(str)+' '+adata_raw.obs.phase.astype(str)
# sc.pl.dotplot(adata_raw, ['CDKN1A', 'TP53'], groupby='tc')
# sc.pl.dotplot(adata_raw, ['CDKN1A', 'TP53'], groupby=['phase','time'])
# adata_raw.to_df().index
# adata_raw.to_df().columns
# p21_53 = adata_raw.to_df()[['CDN1KA','TP53']].copy()
# p21_53 = adata_raw.to_df()[['CDKN1A','TP53']].copy()
# p21_53['time'] = adata_raw.obs.time
# p21_53['time'] = adata_raw.obs.tc
# p21_53.head()
# p21_53 = p21_53.groupby('time').sum()
# p21_53
# p21_53 = p21_53.loc[~p21_53.index.isna()].copy()
# p21_53
# p21_53 = p21_53.loc[p21_53.index!='nan nan'].copy()
# p21_53
# p21_53.index.str.split
# p21_53.index.str.split()
# p21_53.index = p21_53.index.str.split()
# p21_53
# p21_53.CDKN1A.unstack(level=1)
# pd.MultiIndex(p21_53.index)
# pd.MultiIndex?
# pd.MultiIndex.from_product(p21_53.index, names=['t','p'])
# p21_53.index
# pd.MultiIndex.from_arrays(list(p21_53.index))
# p21_53.index = pd.MultiIndex.from_arrays(list(p21_53.index))
# pd.MultiIndex.from_tuples([tuple(i) for i in p21_53.index])
# p21_53.index = pd.MultiIndex.from_tuples([tuple(i) for i in p21_53.index])
# p21_53
# p21_53.CDKN1A.unstack(level=1)
# import statsmodels.api as sm
# p21_ct = sm.stats.Table(p21_53.CDKN1A.unstack(level=1))
# p21_ct
# p21_ct.table_orig
# p21_ct.fittedvalues
# rslt = p21_ct.test_nominal_association()
# rslt.pvalue
# table.chi2_contribs
# p21_ct.chi2_contribs
# sc.pl.dotplot(adata_raw, ['CDKN1A', 'TP53', 'RRM2'], groupby='tc')
# sc.pl.dotplot(adata, ['CDKN1A', 'TP53', 'RRM2'], groupby='tc')
# adata.obs['tc'] = adata.obs.time.astype(str)+' '+adata.obs.phase.astype(str)
# sc.pl.dotplot(adata, ['CDKN1A', 'TP53', 'RRM2'], groupby='tc')
# sc.pl.stacked_violin(adata, ['CDKN1A', 'TP53', 'RRM2'], groupby='tc')
