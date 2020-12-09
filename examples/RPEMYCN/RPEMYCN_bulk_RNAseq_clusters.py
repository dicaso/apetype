import apetype
from apetype.reports import TaskReport, TaskSection
from apetype.tasks import SKIP
import leopard as lp
import zipfile
import os
import pandas as pd
from io import TextIOWrapper
import numpy as np
#from scipy.cluster.vq import vq, kmeans, whiten
from sklearn.preprocessing import quantile_transform
from sklearn import mixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from scipy.stats import pearsonr
from lostdata.dealer.genenames import get_genenames
import gseapy
import matplotlib.pyplot as plt
plt.ion()

# How to do time series analysis?
# - filter low expressing
# - normalization
# - filter non-differential
# - log transform
# - split in groups of similar expression
# - scale within each group to min/max of group
# - take mean of replicates
# - clustering method

class BulkAnalysis(TaskReport):
    """Methods: to analyze the RPEMYCN bulk data, duplicated gene ids
    were removed.  All genes were clustered according to their average
    counts across replicates.  Counts were quantile normalized across
    genes and across samples. Clusters were made with k-means method
    from the Python sklearn package (v0.22.2). To determine the ideal
    number of clusters, the elbow method was used.

    """
    cache = './cache'
    verbose = True
    
    # Report settings
    title = 'Timecourse bulk analysis for RPEMYCN experiment'
    outfile = './bulk_analysis_report_cvn'
    clearpage = True
    
    # CLI settings
    alpha: float = .2
    qnorm_samples: bool = True
    qnorm_genes: bool = True
    only_DE_genes: bool = False
    libsize_analysis: bool = False
    k: int = 10 # number of clusters
    plot_all_genes: bool = False
    gsea: bool = False
    outdir: str = 'results'

    def setup(_) -> type(None):
        if not os.path.exists(_.outdir):
            print('Creating outdir', _.outdir, '...')
            os.mkdir(_.outdir)

    def genenames(_) -> pd.DataFrame:
        from lostdata.dealer.ensembl import get_biomart
        genenames = get_biomart().rename(
            {'Gene stable ID':'ensembl_gene_id', 'Gene name':'symbol'},
            axis=1
        ).set_index('ensembl_gene_id')
        # #Using BDC DE results to extract mapping genenames
        # zf = zipfile.ZipFile('rereportrpemycneranalysiscvn.zip')
        # de_genes = {
        #     f.filename[:-4]: pd.read_table(
        #         TextIOWrapper(zf.open(f)),
        #         usecols = ['Row.names','external_gene_name']
        #     )
        #     for f in zf.infolist()
        #     if f.filename.endswith('.xls')
        # }
        # de_genenames = pd.concat(
        #     de_genes.values()
        # ).drop_duplicates().rename(
        #     {'Row.names':'ensembl_gene_id', 'external_gene_name':'symbol'},
        #     axis=1).set_index('ensembl_gene_id')
        # #de_genenames = de_genenames.loc[~de_genenames.index.duplicated()]

        # # Genenames for downstream analysis
        # # This gave incomplete mapping for unofficial gene symbols
        # genenames = get_genenames()[['ensembl_gene_id', 'symbol']]
        # genenames.set_index('ensembl_gene_id', inplace=True)
        # genenames = pd.concat([genenames, de_genenames])
        genenames = genenames.loc[~genenames.symbol.duplicated()]
        genenames = genenames.loc[~genenames.index.duplicated()]
        print('Removed duplicated gene id, genes with more than one symbol will be biased')
        return genenames

    def counts(_) -> pd.DataFrame:
        zf = zipfile.ZipFile('datasofia.zip')
        with zf.open('CountTable.csv') as f:
            counts = pd.read_csv(TextIOWrapper(f), index_col=0)
            counts_raw = counts
            counts.to_csv('/tmp/counts.csv')
        return counts

    def metadata(_, counts) -> pd.DataFrame:
        zf = zipfile.ZipFile('datasofia.zip')
        with zf.open('metadata.csv') as f:
            metadata = pd.read_csv(TextIOWrapper(f), sep=';')
            metadata.index = metadata['sample # RNA Seq'].str.upper()
            # Sort according to counts
            metadata = metadata.loc[counts.columns].copy()
            metadata['sampletype'] = metadata['Cell line'
                ]+'_'+metadata.Treatment+'_'+metadata.Timepoint
            metadata['libsize'] = counts.sum()
            metadata.to_csv('/tmp/metadata.csv')
        return metadata

    def de_genes_analysis(_) -> str:
        '''R
        library(limma)
        library(edgeR)
        counts = read.csv('/tmp/counts.csv', row.names=1)
        metadata = read.csv('/tmp/metadata.csv', row.names=1)
        counts <- counts[rowSums(counts) > dim(counts)[2], ]
        design = model.matrix(~ 0 + sampletype, data=metadata)
        contrasts = makeContrasts(
          "sampletypeRPEMYCN_4OHT_24h-sampletypeRPEMYCN_UT_0h",
          "sampletypeRPEMYCN_4OHT_48h-sampletypeRPEMYCN_UT_0h",
          "sampletypeRPEMYCN_4OHT_72h-sampletypeRPEMYCN_UT_0h",
          "sampletypeRPEMYCN_4OHT_48h-sampletypeRPEMYCN_4OHT_24h",
          "sampletypeRPEMYCN_4OHT_72h-sampletypeRPEMYCN_4OHT_24h",
          "sampletypeRPEMYCN_4OHT_72h-sampletypeRPEMYCN_4OHT_48h",
          "sampletypeRPEPAR_4OHT_24h-sampletypeRPEPAR_UT_0h",
          "sampletypeRPEPAR_4OHT_48h-sampletypeRPEPAR_UT_0h",
          "sampletypeRPEPAR_4OHT_72h-sampletypeRPEPAR_UT_0h",
          "sampletypeRPEPAR_4OHT_48h-sampletypeRPEPAR_4OHT_24h",
          "sampletypeRPEPAR_4OHT_72h-sampletypeRPEPAR_4OHT_24h",
          "sampletypeRPEPAR_4OHT_72h-sampletypeRPEPAR_4OHT_48h",
          levels=design
        )
        v <- voom(counts, design, plot=TRUE, normalize="quantile")
        fit <- lmFit(v, design)
        fit <- eBayes(fit)
        fit_contrasts <- contrasts.fit(fit, contrasts)
        fit_contrasts <- eBayes(fit_contrasts)

        for(i in 1:ncol(design)){
          print(colnames(design)[i])
          print(sum(topTable(fit, coef=i, adjust="BH", number=10000)$adj.P.Val < .05))
        }
        
        for(i in 1:ncol(contrasts)){
          print(colnames(contrasts)[i])
          print(sum(topTable(fit_contrasts, coef=i, adjust="fdr", number=10000)$adj.P.Val < .05))
        }
        '''
        
    def de_genes(_) -> set:
        # Load DE genes to withhold for analysis
        zf = zipfile.ZipFile('rereportrpemycneranalysiscvn.zip')
        de_genes = {
            f.filename[:-4]: pd.read_table(
                TextIOWrapper(zf.open(f))
            )
            for f in zf.infolist()
            if f.filename.endswith('.xls')
        }
        de_genes = set.union(*[
            set(de_genes[comp]['Row.names'])
            for comp in de_genes
        ])
        return de_genes
                        
        
    def libsize_regression(_, print) -> SKIP(TaskSection):
        # Analyze raw total counts for conditional bias
        # https://songhuiming.github.io/pages/2017/01/21/linear-regression-in-python-chapter-3-regression-with-categorical-predictors/
        import statsmodels.formula.api as smf
        size_pred_data = metadata[
            ['Cell line', 'Treatment', 'Timepoint', 'libsize']
        ].copy() # replicates could be included to check for batch effect
        size_pred_data.Timepoint = size_pred_data.Timepoint.apply(
            lambda x: int(x[:-1])
        )
        size_pred_data.rename(
            {'Cell line': 'Cell_line'}, axis=1, inplace=True
        )
        size_pred_data_norm_ls = size_pred_data.copy()
        size_pred_data_norm_ls.libsize = (
            size_pred_data_norm_ls.libsize -
            size_pred_data_norm_ls[
                size_pred_data_norm_ls.Timepoint == 0
            ].libsize.mean()
        )
        
        reg = smf.ols(
            #formula = "libsize ~ 0+Cell_line+Cell_line:Timepoint", # F-stat 0.994
            formula = "libsize ~ 0+Cell_line:Timepoint", # F-stat 0.00614, R2 0.371 on size_pred_data_norm_ls[size_pred_data.Timepoint < 72]
            #formula = "libsize ~ 0+Cell_line",
            #formula = "libsize ~ Cell_line:Timepoint",
            #formula = "libsize ~ Cell_line",
            #formula = "libsize ~ Cell_line + Cell_line:Timepoint",
            #formula = "libsize ~ Cell_line*Timepoint",
            #formula = "libsize ~ Cell_line + Timepoint",
            #formula = "libsize ~ Cell_line + Treatment + Timepoint",
            data = size_pred_data_norm_ls[size_pred_data.Timepoint < 72]
            #data = size_pred_data[size_pred_data.Timepoint < 72]
            #data = size_pred_data
        ).fit()
        print(reg.summary())
        
        ## Based on results visualize according to cell line
        fig, ax = plt.subplots()
        for gn, g in size_pred_data.groupby('Cell_line'):
            ax.scatter(g.Timepoint, g.libsize, label=gn)
            g_m = g.groupby('Timepoint').mean()
            ax.plot(g_m.index, g_m.libsize)#, c='k')
        ax.legend()
        ax.set_xlabel('Timepoint (h)')
        ax.set_ylabel('Library size (#)')
        ax.set_title('Bias in library sizes')
        return {'figures': {'Library sizes regression': fig}}

    def count_qq(_, counts, metadata) -> pd.DataFrame:
        # filter low expressing
        counts = counts.loc[counts.T.sum() > counts.shape[1]].copy()
        counts_raw = counts.copy()
        mean_gene_counts = counts.T.mean()
        
        # filter non differentially expressed
        if _.only_DE_genes:
            counts = counts.loc[counts.index.isin(de_genes)].copy()
            mean_gene_counts = counts.T.mean()
        
        # normalization
        if _.qnorm_samples:
            counts_norm = pd.DataFrame(
                quantile_transform(counts, n_quantiles=10, axis=0, copy=True),
                columns=counts.columns, index=counts.index
            )
            counts = counts_norm
        if _.qnorm_genes:
            counts_norm = pd.DataFrame(
                quantile_transform(counts, n_quantiles=10, axis=1, copy=True),
                columns=counts.columns, index=counts.index
            )
            counts = counts_norm
            
        # mean of replicates
        count_repmeans = counts.T.groupby(metadata.sampletype).mean().T
        mean_gene_counts = count_repmeans.T.mean()
        
        # log transform
        if not _.qnorm_samples:
            counts_log = (count_repmeans+1).apply(np.log2)
            mean_gene_counts = counts_log.T.mean()
            counts = counts_log
        
        # split in groups of similar expression
        # if not normalizing across genes, in that case this has no sense
        if not _.qnorm_genes:
            mean_gene_counts_qntls = pd.qcut(mean_gene_counts, q=5)
            count_groups = counts.groupby(mean_gene_counts_qntls)
            
            scatter = False
            for name, count_grp in count_groups:
                count_grp = count_grp.rename(mapper = lambda x: x.replace('UT','00'),
                                 axis=1)
                count_grp = count_grp.sort_index(axis=1)
                #whitened = whiten(features) # per feature normalization
                #codebook, distortion = kmeans(count_grp, 10)
                kmeans = KMeans(n_clusters=_.k, random_state=0).fit(count_grp)
                kmeans.cluster_centers_
                
                fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(14,12))
                for clusname, clustergrp in count_grp.groupby(kmeans.labels_):
                    ax = axes.flatten()[clusname-1]
                    if scatter: ax.scatter(
                        list(range(4))*2*clustergrp.shape[0],
                        clustergrp.values.flatten(),
                        c = (['b']*4 + ['g']*4)*clustergrp.shape[0],
                        alpha = _.alpha
                    )
                    else:
                        for l in clustergrp.T:
                            ax.plot(list(range(4)), clustergrp.loc[l][:4], c = 'b', linestyle='dashed', alpha = _.alpha)
                            ax.plot(list(range(4)), clustergrp.loc[l][4:], c = 'g', linestyle='dashed', alpha = _.alpha)
                        ax.plot(
                            list(range(4)),
                            kmeans.cluster_centers_[clusname][:4],
                            c = 'b', linewidth = 5, label='RPEMYCN'
                        )
                        ax.plot(
                            list(range(4)),
                            kmeans.cluster_centers_[clusname][4:],
                            c = 'g', linewidth = 5, label='RPEPAR'
                        )
                    ax.legend()
                fig.suptitle(f'Expression range {name}')
                fig.savefig(os.path.join(_.outdir,f'{name}_clusters.png'))
            
        if _.qnorm_genes:
            # Analysis with only 1 group when counts are also normalized per gene
            
            count_qq = count_repmeans.rename(mapper = lambda x: x.replace('UT','00'),
                             axis=1)
            count_qq = count_qq.sort_index(axis=1)
            return count_qq

    def ideal_k(_, count_qq) -> SKIP(TaskSection):
        # Find ideal number of clusters
        max_k_test = 20
        inertias = []
        silhouettes = []
        print('Calculating optimal k:')
        for k in range(1, max_k_test+1):
            print(k, end='\r')
            fit = KMeans(n_clusters=k, random_state=0).fit(count_qq)
            inertias.append(fit.inertia_)
            if k > 1:
                silhouettes.append(silhouette_score(count_qq, fit.labels_))
                # https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html
                
        fig1, ax = plt.subplots(figsize=(14,12))
        ax.plot(range(1, max_k_test+1), inertias, 'bx-')
        ax.set_xlabel('k')
        ax.set_ylabel('squared distances sum')
        fig2, ax = plt.subplots(figsize=(14,12))
        ax.plot(range(2, max_k_test+1), silhouettes, 'bx-')
        ax.set_xlabel('k')
        ax.set_ylabel('silhouette average')
        return {'figures':{
            'k_elbow.png': fig1,
            'k_silhouette.png': fig2
        }}
            
    def kmeans(_, count_qq) -> KMeans:
        # After visual examination of both plots set ideal k
        ideal_k = _.k
        print('ideal k was set to', ideal_k)
        
        # Get clusters for ideal k
        kmeans = KMeans(n_clusters=ideal_k, random_state=0).fit(count_qq)
        print(kmeans.cluster_centers_)
        return kmeans
    
    def analyzing_the_clusters(_, kmeans, count_qq) -> TaskSection:
        ideal_plot = (2,5) if _.k <= 10 else (4,5)
        fig, axes = plt.subplots(
            nrows=ideal_plot[0], ncols=ideal_plot[1], figsize=(22,12)
        )
        for clusname, clustergrp in count_qq.groupby(kmeans.labels_):
            ax = axes.flatten()[clusname]
            if _.plot_all_genes:
                for l in clustergrp.T:
                    ax.plot(list(range(4)), clustergrp.loc[l][:4], c = 'b', linestyle='dashed', alpha = _.alpha)
                    ax.plot(list(range(4)), clustergrp.loc[l][4:], c = 'g', linestyle='dashed', alpha = _.alpha)
            ax.plot(
                list(range(4)),
                kmeans.cluster_centers_[clusname][:4],
                c = 'b', linewidth = 5, label='RPEMYCN'
            )
            ax.plot(
                list(range(4)),
                kmeans.cluster_centers_[clusname][4:],
                c = 'g', linewidth = 5, label='RPEPAR'
            )
            ax.set_ylim((-.05,1.05))
            ax.set_title(f'{clusname} ({len(clustergrp)})')
        ax.legend()
        return {'figures':{'Quantile normalized samples and genes':fig}}
    
    def cluster_center_correlations(_, kmeans, count_qq, genenames) -> TaskSection:
        tables = []
        # Correlation of genes with cluster center
        for clusname, clustergrp in count_qq.groupby(kmeans.labels_):
            center_corrs = clustergrp.T.apply(lambda x: pearsonr(kmeans.cluster_centers_[clusname], x)[0])
            center_corrs.sort_values(inplace=True, ascending=False)
            center_corrs.index = center_corrs.index.map(genenames.symbol)
            #print('Cluster', clusname, center_corrs.head(10), sep='\n', end='\n\n')
            tables.append((f'Cluster {clusname}', center_corrs))

        return {'tablehead':10, 'tables':tables}
                
        # Determining optimal k with gmm and bic
        # lowest_bic = np.infty
        # bic = []
        # import numpy as np
        # lowest_bic = np.infty
        # cv_types = ['spherical', 'tied', 'diag', 'full']
        # for cv_type in cv_types:
        #     for k in range(1, max_k_test+1):
        #         gmm = mixture.GaussianMixture(n_components=k,
        #             covariance_type=cv_type)
        #         gmm.fit(count_qq)
        #         bic.append(gmm.bic(count_qq))
        #         if bic[-1] < lowest_bic:
        #             lowest_bic = bic[-1]
        #             best_gmm = gmm
        # best_gmm
        # bic = np.array(bic)
        # it
        # import itertools as it
        # color_iter = it.cycle(['navy', 'turquoise', 'cornflowerblue',
        #                               'darkorange'])
        # clf = best_gmm
        # bars = []
        # fig, ax = plt.subplots(figsize=(22,12))
        # for i, (cv_type, color) in enumerate(zip(cv_types, color_iter)):
        #     xpos = np.array(range(1, max_k_test+1))+ .2 * (i - 2)
        #     bars.append(ax.bar(xpos, bic[i * len(range(1, max_k_test+1)):(i + 1) *
        #  len(range(1, max_k_test+1))],width=.2, color=color))                    
        # ax.set_xticks((range(1, max_k_test+1)))
        # ax.set_ylim([bic.min() * 1.01 - .01 * bic.max(), bic.max()])
        # ax.set_title('BIC score per model')
        # n_components_range = range(1, max_k_test+1)
        # xpos = np.mod(bic.argmin(), len(n_components_range)) + .65 +\
        #     .2 * np.floor(bic.argmin() / len(n_components_range))
        # ax.text(xpos, bic.min() * 0.97 + .03 * bic.max(), '*', fontsize=14)
        # ax.set_xlabel('Number of components')
        # ax.legend([b[0] for b in bars], cv_types)
        # Plot clusters
        # labels = best_gmm.predict(count_qq)
        # for clusname, clustergrp in count_qq.groupby(labels):
        #     ax = axes.flatten()[clusname]
        #     ax.plot(
        #        list(range(4)),
        #        best_gmm.means_[clusname][:4],
        #        c = 'b', linewidth = 5, label='RPEMYCN'
        #     )
        #     ax.plot(
        #        list(range(4)),
        #        best_gmm.means_[clusname][4:],
        #        c = 'g', linewidth = 5, label='RPEPAR'
        #     )
        #     ax.set_ylim((-.05,1.05))
        #     ax.set_title(f'{clusname} ({len(clustergrp)})')
        #     ax.legend()
        
    def geneset_enrichment_analysis(_, genenames, counts, metadata) -> TaskSection:
        #avail_genesets = gseapy.get_library_name()
        genesets = [
            'ChEA 2016',
            #'ENCODE and ChEA Consensus TFs from ChIP-X', # too many overl sig gs
            'ARCHS4 TFs Coexp',
            'BioPlanet 2019',
            'GO Biological Process 2018',
            'TF Perturbations Followed by Expression',
            #'TRRUST Transcription Factors 2019',
            'lncHUB lncRNA Co-Expression',
            #'Enrichr Submissions TF-Gene Coocurrence',
            #'TRANSFAC and JASPAR PWMs',
            #'Epigenomics Roadmap HM ChIP-seq',
            #'TargetScan microRNA 2017',
            #'miRTarBase 2017',
            #'ENCODE TF ChIP-seq 2015',
            #'TF-LOF Expression from GEO',
            #'ENCODE Histone Modifications 2015',
            #'Transcription Factor PPIs',
            #'Genome Browser PWMs',
            #'Cancer_Cell_Line_Encyclopedia', # only 1 cluster with sig gs
            #'Chromosome_Location',           # no sig gs
            #'Jensen_DISEASES',               # 4 clusters with sig gs
            'KEGG_2019_Human',
            'MSigDB_Oncogenic_Signatures',
            'TargetScan_microRNA',
            'TRANSFAC_and_JASPAR_PWMs'
        ]
        genesets = [gs.replace(' ','_') for gs in genesets]

        section = {'figures':{}}
        if _.gsea:
            gseas = {}
            fig, ax = plt.subplots(figsize=(14,12))
            for clusname, clustergrp in count_qq.groupby(kmeans.labels_):
                grp_gnames = [
                    genenames.loc[i].symbol
                    for i in clustergrp.index
                    if i in genenames.index
                ]
                try:
                    gsea_results = gseapy.enrichr(
                        gene_list = grp_gnames,
                        description = 'pathway',
                        gene_sets = ','.join(genesets),
                        outdir=f'results/cluster_{clusname}'
                    )
                    gseas[clusname] = gsea_results
                    print(
                        clusname, 'enriched gene sets',
                        (gsea_results.results['Adjusted P-value']<.05).sum()
                    )
                except: print('Something went wrong GSEA for', clusname)
            
                with open(
                        os.path.join(
                            _.outdir,
                            f'cluster_{clusname}_genes.txt'),'wt'
                ) as f:
                    f.write('\n'.join(grp_gnames))
                
                ax.plot(
                        list(range(4))[::-1]+list(range(4)),
                        list(kmeans.cluster_centers_[clusname][:4])[::-1]+
                        list(kmeans.cluster_centers_[clusname][4:]),
                        linewidth = 5, label=clusname
                    )
            ax.legend()
            section['figures']['clusters_together'] = fig
            
            ## Analyze overlap in genesets
            ### retrieve signficant genesets
            gsea_sig = {
                c:gseas[c].results.loc[
                    gseas[c].results['Adjusted P-value']<.05
                ].copy() for c in gseas
            }
            overlap = np.array(
                [
                    [
                        gsea_sig[c].Term.isin(gsea_sig[o].Term).sum()
                        for o in sorted(gsea_sig)
                    ] for c in gsea_sig
                ]
            )
            print(overlap)
            
            ### overlap per geneset library
            fig, axes = plt.subplots(nrows=2, ncols=5, figsize=(22,12))
            for i,gs in enumerate(genesets):
                ax = axes.flatten()[i]
                gsea_lib = {
                    c:gsea_sig[c].loc[gsea_sig[c].Gene_set == gs].copy()
                    for c in gsea_sig
                }
                overlap_lib = np.array(
                    [
                        [
                            gsea_lib[c].Term.isin(gsea_lib[o].Term).sum()
                        for o in sorted(gsea_lib)
                        ] for c in gsea_lib
                    ]
                )
                #fig, ax = plt.subplots(figsize=(14,12))
                ax.matshow(np.nan_to_num(
                    overlap_lib/overlap_lib.max(1),
                    0
                ))
                ax.set_title(gs.replace('_', ' '))
                #print(gs)
                #print(overlap_lib)
                #input()
            section['figures']['overlap_genesets'] = fig
            
        # enrichr api to check instead of using gseapy
        # https://github.com/MaayanLab/enrichr_api/blob/master/main.py
        # import json
        # import requests
        # from time import sleep
        # def get_enrichr_results(gene_set_library, genelist, description):
        #     ADDLIST_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
        #     payload = {
        #         'list': (None, genelist),
        #         'description': (None, description)
        #     }
        
        #     response = requests.post(ADDLIST_URL, files=payload)
        #     if not response.ok:
        #         raise Exception('Error analyzing gene list')
        #     sleep(1)
        #     data = json.loads(response.text)
        
        #     RESULTS_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
        #     query_string = '?userListId=%s&backgroundType=%s'
        #     user_list_id = data['userListId']
        #     response = requests.get(RESULTS_URL + query_string % (user_list_id, gene_set_library))
        #     sleep(1)
        #     return [data['shortId'], json.loads(response.text)]
        
        
        # def get_libraries():
        #     libs_json = json.loads(requests.get('http://amp.pharm.mssm.edu/Enrichr/datasetStatistics').text)
        #     libs = [lib['libraryName']for lib in libs_json['statistics']]
        #     return libs
        
        # Plots for individual genes
        genesOfInterest = {
            'preP21':[
                'PRDX1', 'EXOSC8', 'PSMC4', 'EIF4A3', 'HSP90AA1', 'PCNA',
                'SLC3A2', 'CCT8', 'EMC2', 'PCMT1', 'ESD', 'FDPS', 'HPF1',
                'PSMA5', 'NPM1', 'FH', 'EIF6', 'FRG1'#, 'OSLA3'
            ],
            'P21':[
                'MT-CO2', 'MT-CO1', 'MT-ND4', 'MT-CO3', 'MT-CYB', 'RPS12',
                'RPLP1', 'RPS28', 'RPL37A', 'TPT1', 'RPS2', 'RPS27', 'EEF1A1',
                'MT-ATP6', 'RPL12', 'RPL34', 'RPL38', 'RPL32', 'RPL29', 'RPS6'
            ]
        }

        genenames_symbol_index = genenames.reset_index().set_index('symbol')
        
        for gset in genesOfInterest:
            fig, axes = plt.subplots(nrows=4,ncols=5,figsize=(14,12), sharex=True, sharey=True)
            for (g,ax) in zip(genesOfInterest[gset],axes.flatten()):
                egid = genenames_symbol_index.loc[g].ensembl_gene_id
                genecounts = pd.DataFrame({
                    'counts': counts.loc[egid],
                    'time': metadata.Timepoint.apply(lambda x: int(x[:-1])),
                    'cellline': metadata['Cell line']
                })
                for n,grp in genecounts.groupby('cellline'):
                    ax.scatter(grp.time, grp.counts, label=n)
                ax.set_title(g)
            ax.legend()
            fig.suptitle(gset)
            section['figures'][gset] = fig
        return section
        
    def integrating_with_single_cell_analysis(_, counts, metadata, genenames) -> TaskSection:
        """From the single cell data analysis, I took the diffxpy quantitative signatures
        to see how they respond in the bulk data.
        """
        genenames_symbol_index = genenames.reset_index().set_index('symbol')
        sc_sig_genes = pd.read_csv('../RPE_MYCN_10x/analysis_results_cvn/222324quantmodel.csv', index_col=0)
        print('Significant genes from single cell quant analysis across 22, 23, and 24 hours:', len(sc_sig_genes))
        sc_sig_genes['up'] = sc_sig_genes.log2fc > 0
        sc_sig_up_genes = sc_sig_genes.gene[sc_sig_genes.up]
        sc_sig_down_genes = sc_sig_genes.gene[~sc_sig_genes.up]
        sc_sig_up_genes_egid = [
            genenames_symbol_index.loc[g].ensembl_gene_id
            for g in sc_sig_up_genes
            if g in genenames_symbol_index.index
        ]
        sc_sig_down_genes_egid = [
            genenames_symbol_index.loc[g].ensembl_gene_id
            for g in sc_sig_down_genes
            if g in genenames_symbol_index.index
        ]
        
        for gset in (sc_sig_down_genes[:20],sc_sig_up_genes[:20]):
            fig, axes = plt.subplots(nrows=4,ncols=5,figsize=(14,12), sharex=True, sharey=True)
            for (g,ax) in zip(gset,axes.flatten()):
                try: egid = genenames_symbol_index.loc[g].ensembl_gene_id
                except KeyError:
                    print(g)
                    continue
                genecounts = pd.DataFrame({
                    'counts': counts.loc[egid],
                    'time': metadata.Timepoint.apply(lambda x: int(x[:-1])),
                    'cellline': metadata['Cell line']
                })
                for n,grp in genecounts.groupby('cellline'):
                    ax.scatter(grp.time, grp.counts, label=n)
                ax.set_title(g)
            ax.legend()
        
        count_ranks = counts.rank()
        count_ranks_up = count_ranks.loc[count_ranks.index.isin(sc_sig_up_genes_egid)]
        ranksums = pd.DataFrame({
                    'rsums': count_ranks_up.sum(),
                    'time': metadata.Timepoint.apply(lambda x: int(x[:-1])),
                    'cellline': metadata['Cell line']
                })
        fig1, ax = plt.subplots()
        for n,grp in ranksums.groupby('cellline'):
            ax.scatter(grp.time, grp.rsums, label=n)
        ax.set_title('Ranksums up')
        ax.legend()
        
        count_ranks_down = count_ranks.loc[count_ranks.index.isin(sc_sig_down_genes_egid)]
        ranksums = pd.DataFrame({
                    'rsums': count_ranks_down.sum(),
                    'time': metadata.Timepoint.apply(lambda x: int(x[:-1])),
                    'cellline': metadata['Cell line']
                })
        fig2, ax = plt.subplots()
        for n,grp in ranksums.groupby('cellline'):
            ax.scatter(grp.time, grp.rsums, label=n)
        ax.set_title('Ranksums down')
        ax.legend()
        return [
            ('Ranksums up', fig1),
            ('Ranksums down', fig2),
            ('clearpage', True)
        ]

if __name__ == '__main__':
    analysis = BulkAnalysis(parse=True)
    analysis.run()#load_cache=True)
    analysis.report.outputZip('svg')
    plt.close('all')
