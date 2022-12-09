library('vroom')
library('foreach')
library('doMC')
library('survminer')
library('grid')
library('survival')

tcga_key = list('AML'='aml', 'bladder'='blca', 'breast'='brca', 'cervical'='cesc', 
                'colon'='coad', 'esophageal'='esca', 'GBM'='gbm', 'glioma'='glm', 
                'head_neck'='hnsc', 'kidney'='kicc', 'liver'='lihc', 'lung_ad'='luad',
                'lung_sq'='lusq', 'melanoma'='skcm','ovarian'='ovsc', 'pancreatic'='paad', 
                'paraganglioma'='pcpg', 'prostate'='prad', 'rectal'='read', 'sarcoma'='sarc', 
                'stomach'='stad', 'testicular'='tgct', 'thymoma'='thym', 'thyroid'='thca', 
                'uterine_endometrial'='ucec')
gene_map = read.csv('data/gene_syn_map',sep='\t',header=F)

col_surv = rev(c(rgb(0.86, 0.2, 0.3, 0.75), rgb(0, 0.5, 1, 0.75)))
cmd_mkdir = "mkdir -p DIR"


calculateMI_v2 <- function(x, y) {
        #x = cMatrix[,1] ; y = cMatrix[,2] ; # = rowSums(cMatrix)
        #x[is.na(x)] = 0 ; y[is.na(y)] = 0
        rs_cMatrix = x + y ; num_genes = sum(rs_cMatrix)
        log_t1 = log((num_genes*x)/(sum(x) * rs_cMatrix)) ; log_t1[is.infinite(log_t1)] = 0
        log_t2 = log((num_genes*y)/(sum(y) * rs_cMatrix)) ; log_t2[is.infinite(log_t2)] = 0
        val_ret = sum((1/num_genes)* ((x)*log_t1 + (y)*log_t2))
}

getSurv <- function(all_clin_df, surv_type = 'OS', rand_iter = 1000, samp_name = 'MPS+', ctrl_name = 'MPS-') {
    library('survival')
    all_clin = all_clin_df
    s_name = samp_name #; if(s_name == 'positive') { s_name = 'MPS+' }
    c_name = ctrl_name #; if(c_name == 'negative') { c_name = 'MPS–' }
    if (surv_type == 'OS') { all_clin$mod_censor_time = as.numeric(all_clin$OS.time) ; all_clin$death_event_binary = as.numeric(all_clin$OS) }
    if (surv_type =='PFI') { all_clin$mod_censor_time =as.numeric(all_clin$PFI.time) ; all_clin$death_event_binary = as.numeric(all_clin$PFI) }
    len_s = sum(all_clin$MPS_groups == samp_name) ; len_c = sum(all_clin$MPS_groups == ctrl_name)
    s_all = survfit(Surv(as.numeric(all_clin$mod_censor_time), all_clin$death_event_binary) ~ 1)
    ss = Surv(as.numeric(as.character(all_clin$mod_censor_time)),all_clin$death_event_binary)
    val_ret = NA
    pv = s_fit = pv_cox = hz_ratio = s_coxph = z = c_ix = NA
    if (! (sum(is.na(ss)) == length(is.na(ss)))) {
        s_fit = survfit(Surv(mod_censor_time, death_event_binary) ~ MPS_groups, data = all_clin)
        s_diff = tryCatch({ 
            survdiff(Surv(mod_censor_time, death_event_binary) ~ MPS_groups, data = all_clin)
            }, error = function(err) { list(chisq=NA, n=dim(all_clin)[1]) } )
        pv <- ifelse(is.na(s_diff),1,(round(1 - pchisq(s_diff$chisq, length(s_diff$n) - 1),50)))[[1]]
        
        s_coxph = coxph(Surv(mod_censor_time, death_event_binary) ~ MPS_groups, data = all_clin)
        if (! is.na(s_coxph[['coefficients']])) { 
            s_c = summary(s_coxph)
            hz_ratio = as.numeric(data.frame(s_c$coefficients)['exp.coef.'])
            pv_cox = as.numeric(data.frame(s_c$coefficients)['Pr...z..'])
            z = as.numeric(data.frame(s_c$coefficients)['z'])
            c_ix = as.numeric(summary(s_coxph)$concordance['C'])
        }

        pv_r = c()
        #all_clin_r = all_clin ; prob_v = as.numeric(prop.table(table(all_clin$group)))
        #grp_r = rand_iter*sample(c(s_name, c_name), dim(all_clin)[1], replace=T, prob=prob_v)
        
        for(r_i in seq(rand_iter)) {
            all_clin_r = all_clin ; prob_v = as.numeric(prop.table(table(all_clin$MPS_groups)))
            all_clin_r$MPS_groups = sample(c(s_name, c_name), dim(all_clin)[1], replace=T, prob=prob_v)
            s1_r = tryCatch({ 
            survdiff(Surv(mod_censor_time, death_event_binary) ~ MPS_groups, data = all_clin_r)
            }, error = function(err) { list(chisq=NA, n=dim(all_clin_r)[1]) } )
            pv_r = c(pv_r, round(1 - pchisq(s1_r$chisq, length(s1_r$n) - 1),50))
        }
        fp = (sum(pv_r < pv)/rand_iter)
        median_surv = summary(s_fit)$table[,'median']
#        med_ab = round(as.numeric(median_surv[grep(samp_name,names(median_surv))]),1)
#        med_re = round(as.numeric(median_surv[grep(ctrl_name,names(median_surv))]),1)
        med_ab = round(as.numeric(median_surv[paste('MPS_groups=',s_name,sep='')]),1)
        med_re = round(as.numeric(median_surv[paste('MPS_groups=',c_name,sep='')]),1)
        m_sub = paste(s_name,' ', len_s, '(', med_ab,' mo)',
                                    '| ', c_name,' ', len_c, '(', med_re, ' mo)', sep='')
    
        out_df = data.frame(pval_surv = pv, pval_cox = pv_cox, C_index = c_ix, z_cox = z,HZR_cox = hz_ratio, info_survival = m_sub)
   
#        s_dt = survfit(Surv(mod_censor_time, death_event_binary) ~ (group+disease_type), data = all_clin)
#        s1_dt = tryCatch({ 
#            survdiff(Surv(as.numeric(as.character(all_clin$mod_censor_time)),all_clin$death_event_binary) ~ all_clin$group + all_clin$disease_type)
#            }, error = function(err) { list(chisq=NA, n=dim(all_clin)[1]) } )
#        pv_distype <- ifelse(is.na(s1),1,(round(1 - pchisq(s1$chisq, length(s1$n) - 1),50)))[[1]]
        val_ret = list(KM_pv = pv, KM_fit = s_fit, clin_data = all_clin, cox_pv = pv_cox, hzr = hz_ratio, 
                                cox_fit = s_coxph, avg_fit = s_all, rand_pv = pv_r, z=z,
                                concordance = c_ix, out_df = out_df, fdr = fp)
    }
}

plot_Surv <- function(clin_data) {
    all_clin = clin_data
    dis_ = as.character(unique(all_clin$disease_type))
    path_id = as.character(unique(all_clin$module_id))
    name_module = unique(as.character(all_clin$module_name))
    name_module_trim = strtrim(gsub('at','-',gsub('OUTPUT','',gsub('-','',gsub('_', '', unique(as.character(all_clin$module_name)))))),17)
    genes_in_module = sort(unique(unlist(strsplit(as.character(all_clin$genes_in_module), '\\|'))))
    num_genes_in_module = unique(as.numeric(all_clin$num_genes_in_module))
    pl_o = pl_p = nullGrob()
    if (sum(table(all_clin$MPS_groups) > 10) == length(unique(all_clin$MPS_groups))) {
        surv_ovs = getSurv(all_clin, 'OS',1) ; surv_pfs = getSurv(all_clin, 'PFI',1)
        fdr_ovs = as.numeric(surv_ovs[[12]]) ; fdr_pfs = as.numeric(surv_pfs[[12]])
        if (fdr_ovs == 0) { f_ovs = '1e-3' } ; if (fdr_ovs != 0) { f_ovs = formatC(fdr_ovs,format='e',digits=1) }
        if (fdr_pfs == 0) { f_pfs = '1e-3' } ; if (fdr_pfs != 0) { f_pfs = formatC(fdr_pfs,format='e',digits=1) }

        m_ovs= paste(dis_, 'p=', formatC(as.numeric(surv_ovs[[11]]$pval_surv),format='e',digits=1), 
                     ' | fdr<', f_ovs,' | HR=', signif(as.numeric(surv_ovs[[11]]$HZR_cox), 2))
        m_ovs= paste(dis_, 'p=', formatC(as.numeric(surv_ovs[[11]]$pval_surv),format='e',digits=1), 
                     ' | HR=', signif(as.numeric(surv_ovs[[11]]$HZR_cox), 2))
        m_sub_ovs = paste('OVS (', num_genes_in_module,' genes)  ',name_module_trim, sep='')
        m_pfs= paste(dis_, 'p=', formatC(as.numeric(surv_pfs[[11]]$pval_surv),format='e',digits=1), 
                     ' | fdr<', f_pfs,' | HR=', signif(as.numeric(surv_pfs[[11]]$HZR_cox), 2))
        m_pfs= paste(dis_, 'p=', formatC(as.numeric(surv_pfs[[11]]$pval_surv),format='e',digits=1), 
                     ' | HR=', signif(as.numeric(surv_pfs[[11]]$HZR_cox), 2))
        m_sub_pfs = paste('PFS (', num_genes_in_module,' genes). ' ,name_module_trim, sep='')
       
        pl_ovs = ggsurvplot(surv_ovs[[2]], surv_ovs[[3]], pval=F, title=m_ovs, font.title=12, censor.shape=124,censor.size=2,
                            subtitle=m_sub_ovs, font.subtitle=10,surv.median.line='hv', palette =col_surv, risk.table=F)
    
        pl_pfs = ggsurvplot(surv_pfs[[2]], surv_pfs[[3]], pval=F, title=m_pfs, font.title=12, censor.shape=124,censor.size=2,
                            subtitle=m_sub_pfs, font.subtitle=10,surv.median.line='hv', palette =col_surv, risk.table=F)
        pl_o = pl_ovs$plot ; pl_p = pl_pfs$plot
        if(as.character(unique(all_clin$disease_type)) == 'AML') { pl_p = nullGrob() }
    }
    val_ret = list('ovs' = pl_o, 'pfs' = pl_p)
}

get_geneExp_info <- function(clin_data) {
    all_clin = clin_data
    dis_ = as.character(unique(all_clin$disease_type))
    s_p = as.character(all_clin$sample_analyzed[which(all_clin$MPS_groups == 'MPS+')])
    s_n = as.character(all_clin$sample_analyzed[which(all_clin$MPS_groups == 'MPS-')])
    path_id = as.character(unique(all_clin$module_id))
    name_module = unique(as.character(all_clin$module_name))
    name_module_trim = strtrim(gsub('at','-',gsub('OUTPUT','',gsub('-','',gsub('_', '', unique(as.character(all_clin$module_name)))))),17)
    genes_in_module = sort(unique(unlist(strsplit(as.character(all_clin$genes_in_module), '\\|'))))
    num_genes_in_module = unique(as.numeric(all_clin$num_genes_in_module))
    
    d_genev = vroom(paste('data/disease_gene_expression/', dis_, '_tumdat_zscore', sep=''), delim = "\t",col_names=T)

    d_gene_module = as.data.frame(d_genev[which(d_genev$gene %in% genes_in_module),])
    rownames(d_gene_module) = as.character(d_gene_module$gene)
    d_gene_module = d_gene_module[,-grep('gene',colnames(d_gene_module))]
    genes_in_module = intersect(rownames(d_gene_module), genes_in_module)
    dg_p = d_gene_module[,s_p] ; dg_n = d_gene_module[,s_n]
    pv_vec = sapply(genes_in_module, function(g_) { 
                                                    v_ab = as.numeric(dg_p[g_,]) ; 
                                                    v_be = as.numeric(dg_n[g_,]) ; 
                                                    w_ = wilcox.test(v_ab, v_be) ; w_$p.value 
                                                })
    pv_vec[is.na(pv_vec)] = 1 ; fdr_vec = p.adjust(pv_vec, 'fdr')
    df_g = data.frame(gene = names(pv_vec), pv = pv_vec, fdr = fdr_vec)
    
    d_genev = readRDS(paste('data/disease_gene_expression/', dis_, '_log2dat_rds', sep=''))
    #intersect(as.character(d_genev, genes_in_module))
    d_gene_module = d_genev[intersect(as.character(d_genev$gene), genes_in_module), grep('TCGA',colnames(d_genev))]
    #median_gene_module = as.numeric(apply(d_gene_module[,c(s_p, s_n)], 2, mean))

    dg_p = d_gene_module[,s_p] ; dg_n = d_gene_module[,s_n]
    lR_vec = (apply(dg_p, 1, median) - apply(dg_n, 1, median))
    tmp_ = data.frame(gene = names(lR_vec), lR = as.numeric(lR_vec)) ; df_g = merge(df_g, tmp_)
    df_g$path_id = path_id ; df_g$cohort = dis_ ; rm(tmp_)
    eps_th = min(c(1e-16, min(df_g$fdr[df_g$fdr !=0]))) ; df_g$fdr[which(df_g$fdr == 0)] = eps_th
    val_ret = df_g
}


plot_geneExp <- function(gene_DE_data) {
    df_g = gene_DE_data
    pl_volc = nullGrob()
    if (dim(df_g)[1] > 25) {
        pl_volc = ggplot(df_g, aes(lR, -log10(fdr))) + geom_point(alpha=0.25) + theme_classic() +
                    geom_hline(yintercept = -log10(0.01),col=rgb(1,0,0,0.2)) + geom_vline(xintercept=0, col=rgb(1,0,0,0.2)) +
                    xlab('log2-ratio MPS+ vs. MPS-') + ylab('-log10(fdr)')

    }
    val_ret = list(pl_volc)
}

getSurv_locus <- function(gene_name, disease_name, rand_iter = 1000) {
    out_gene_dir = paste('results/single_locus/',sep='') ; system(gsub('DIR', out_gene_dir, cmd_mkdir))
    ovs_surv_gene = pfs_surv_gene = 'NA'
    gene_sym = as.character(gene_map[which(tolower(gene_map$V1) == tolower(gene_name)),]$V2)
    if (length(gene_sym) != 0) {
        dis_clin = readRDS(paste('data/disease_clinical/', disease_name, '_clin_rds',sep=''))
        ## SNV
        dis_snv = readRDS(paste('data/disease_gene_snv/', disease_name, '_snv_rds',sep='')) ; g_snv = rownames(dis_snv) ; samp_snv = colnames(dis_snv)
        ## CNA
        dis_cna = readRDS(paste('data/disease_gene_snv/', disease_name, '_cna_rds',sep='')) ; g_cna = rownames(dis_cna) ; samp_cna = colnames(dis_cna)
        ## expression
        dis_exp = as.data.frame(vroom(paste('data/disease_gene_expression/', disease_name, '_tumdat_zscore', sep=''), delim = "\t",col_names=T))
        rownames(dis_exp) = as.character(dis_exp$gene) ; dis_exp = dis_exp[,-grep('gene', colnames(dis_exp))] ; g_exp = rownames(dis_exp) ; samp_exp = colnames(dis_exp)
        ovs_snv = pfs_snv = ovs_cna = pfs_cna = ovs_exp = pfs_exp = c()
        if (gene_sym %in% g_snv) {
            mut_samp = samp_snv[which(as.numeric(dis_snv[gene_sym,]) == 1)] ; wt_samp = samp_snv[which(as.numeric(dis_snv[gene_sym,]) != 1)]
            clin_snv = merge(rbind(data.frame(sample_TCGA = mut_samp, MPS_groups = rep('mutated', length(mut_samp))), 
                                   data.frame(sample_TCGA = wt_samp, MPS_groups = rep('not_mutated', length(wt_samp)))), dis_clin)
            ovs_surv_snv = getSurv(clin_snv, 'OS', rand_iter, 'mutated', 'not_mutated') ; pfs_surv_snv = getSurv(clin_snv, 'PFI', rand_iter, 'mutated', 'not_mutated')
            ovs_snv = ovs_surv_snv[['out_df']] ; ovs_snv$gene_symbol = gene_sym ; ovs_snv$gene = gene_name
            ovs_snv$cohort = disease_name ; ovs_snv$survival_type = 'OS' ; ovs_snv$fdr = ovs_surv_snv[['fdr']] ; ovs_snv$category = 'SNV'
            pfs_snv = ovs_surv_snv[['out_df']] ; pfs_snv$gene_symbol = gene_sym ; pfs_snv$gene = gene_name
            pfs_snv$cohort = disease_name ; pfs_snv$survival_type = 'PFI' ; pfs_snv$fdr = pfs_surv_snv[['fdr']] ; pfs_snv$category = 'SNV'
        }
        if (gene_sym %in% g_cna) {
            pos_samp = samp_cna[which(as.numeric(dis_cna[gene_sym,]) > 0)] ; neg_samp = samp_cna[which(as.numeric(dis_cna[gene_sym,]) < 0)]
            clin_cna = merge(rbind(data.frame(sample_TCGA = pos_samp, MPS_groups = rep('cna_positive', length(pos_samp))), 
                                   data.frame(sample_TCGA = neg_samp, MPS_groups = rep('cna_negative', length(neg_samp)))), dis_clin)
            ovs_surv_cna = getSurv(clin_cna, 'OS', rand_iter, 'cna_positive', 'cna_negative') ; pfs_surv_cna = getSurv(clin_cna, 'PFI', rand_iter, 'cna_positive', 'cna_negative')
            ovs_cna = ovs_surv_cna[['out_df']] ; ovs_cna$gene_symbol = gene_sym ; ovs_cna$gene = gene_name
            ovs_cna$cohort = disease_name ; ovs_cna$survival_type = 'OS' ; ovs_cna$fdr = ovs_surv_cna[['fdr']] ; ovs_cna$category = 'CNA'
            pfs_cna = ovs_surv_cna[['out_df']] ; pfs_cna$gene_symbol = gene_sym ; pfs_cna$gene = gene_name
            pfs_cna$cohort = disease_name ; pfs_cna$survival_type = 'PFI' ; pfs_cna$fdr = pfs_surv_cna[['fdr']] ; pfs_cna$category = 'CNA'
        }
        if (gene_sym %in% g_exp) {
            pos_samp = samp_exp[which(as.numeric(dis_exp[gene_sym,]) > 0)] ; neg_samp = samp_exp[which(as.numeric(dis_exp[gene_sym,]) < 0)]
            clin_exp = merge(rbind(data.frame(sample_analyzed = pos_samp, MPS_groups = rep('expression_positive', length(pos_samp))), 
                                   data.frame(sample_analyzed = neg_samp, MPS_groups = rep('expression_negative', length(neg_samp)))), dis_clin)
            ovs_surv_exp = getSurv(clin_exp, 'OS', rand_iter, 'expression_positive', 'expression_negative') ; pfs_surv_exp = getSurv(clin_exp, 'PFI', rand_iter, 'expression_positive', 'expression_negative')
            ovs_exp = ovs_surv_exp[['out_df']] ; ovs_exp$gene_symbol = gene_sym ; ovs_exp$gene = gene_name
            ovs_exp$cohort = disease_name ; ovs_exp$survival_type = 'OS' ; ovs_exp$fdr = ovs_surv_exp[['fdr']] ; ovs_exp$category = 'expression'
            pfs_exp = ovs_surv_exp[['out_df']] ; pfs_exp$gene_symbol = gene_sym ; pfs_exp$gene = gene_name
            pfs_exp$cohort = disease_name ; pfs_exp$survival_type = 'PFI' ; pfs_exp$fdr = pfs_surv_exp[['fdr']] ; pfs_exp$category = 'expression'
        }
        ovs_surv_gene = rbind(ovs_snv, ovs_cna, ovs_exp) ; pfs_surv_gene = rbind(pfs_snv, pfs_cna, pfs_exp)
    }
    write.table(ovs_surv_gene, paste(out_gene_dir, disease_name,'_', gene_name,'_OVS.xls',sep=''),sep='\t', quote=F,row.names=F)
    write.table(pfs_surv_gene, paste(out_gene_dir, disease_name,'_', gene_name,'_PFS.xls',sep=''),sep='\t', quote=F,row.names=F)
    val_ret = list('ovs' = ovs_surv_gene, 'pfs' = pfs_surv_gene)
}

plot_Immune <- function(clin_data) {
    all_clin = clin_data ; all_clin = all_clin[! is.na(all_clin$MPS_groups),]
    dis_ = as.character(unique(all_clin$disease_type))
    path_id = as.character(unique(all_clin$module_id))
    name_module = unique(as.character(all_clin$module_name))
    name_module_trim = strtrim(gsub('at','-',gsub('OUTPUT','',gsub('-','',gsub('_', '', unique(as.character(all_clin$module_name)))))),17)
    genes_in_module = sort(unique(unlist(strsplit(as.character(all_clin$genes_in_module), '\\|'))))
    num_genes_in_module = unique(as.numeric(all_clin$num_genes_in_module))
	d_imm = readRDS(paste('data/disease_immune/', dis_, '_immune_rds', sep=''))
    pl_imm = nullGrob()
    if (sum(table(all_clin$MPS_groups) > 10) == length(unique(all_clin$MPS_groups))) {
        df_imm = merge(d_imm, all_clin[,c('sample_analyzed', 'MPS_groups')])
        w_m2 = wilcox.test(df_imm[df_imm$MPS_groups == 'MPS+',]['IMMUNE_Macrophages.M2'][,1], 
                    df_imm[df_imm$MPS_groups == 'MPS-',]['IMMUNE_Macrophages.M2'][,1])$p.value
        w_cd8 = wilcox.test(df_imm[df_imm$MPS_groups == 'MPS+',]['IMMUNE_T.cells.CD8'][,1], 
                    df_imm[df_imm$MPS_groups == 'MPS-',]['IMMUNE_T.cells.CD8'][,1])$p.value
        
        colnames(df_imm) = gsub('IMMUNE_T.cells.CD8', paste('CD8+ T-cells (',formatC(w_cd8,format='e',digits=0),')',sep=''),
                           gsub('IMMUNE_Macrophages.M2', paste('M2 macrophages (',formatC(w_m2,format='e',digits=0),')',sep=''), 
                           colnames(df_imm)))
        melt_df_imm = reshape2::melt(df_imm)
        pl_imm = ggplot(melt_df_imm, aes(x=MPS_groups, y=value, fill=MPS_groups)) + scale_fill_manual(values=col_surv) + 
                    geom_violin(alpha=0.5, draw_quantiles=c(0.25, 0.5, 0.75)) + #scale_fill_hue(l=40, c=35) + 
                    labs(title='') +
                    facet_wrap(~variable) + theme_bw() #theme_classic()
    }
    val_ret = list(pl_imm)
}


get_commonHistopath <- function(disease_name) {
    dis_ = disease_name
    hist_v = c('histological_type', 'coarse_pathologic_stage', 'age_groups')
    if (dis_ == 'breast') {
        hist_v = c('hormone_receptor_status', 'HER2_status', 'TNBC_status',
                   'histological_type', 'coarse_pathologic_stage', 'age_groups')
    }
    if (dis_ == 'prostate') {
        hist_v = c('gleason_score','PSA_value','histological_type', 
                   'coarse_pathologic_stage', 'age_groups')
    }
    if (dis_ == 'AML') {
        hist_v = c('mol_test_status','histological_type', 
                   'coarse_pathologic_stage', 'age_groups')
    }
    if (dis_ == 'head_neck') {
        hist_v = c('HPV_status','histological_type', 
                   'coarse_pathologic_stage', 'age_groups')
    }
    if (dis_ == 'cervical') {
        hist_v = c('HPV_status','histological_type', 
                   'coarse_pathologic_stage', 'age_groups')
    }
    if (dis_ == 'colon') {
        hist_v = c('MSI_status','tumor_side', 'histological_type', 
                   'coarse_pathologic_stage', 'age_groups')
    }
    val_ret = hist_v
}

plot_Histopath <- function(clin_data) {
    all_clin = clin_data ; all_clin = all_clin[! is.na(all_clin$MPS_groups),]
    dis_ = as.character(unique(all_clin$disease_type)) ; path_id = as.character(unique(all_clin$path_id))
    var_vec = c(get_commonHistopath(dis_)) #c('histological_type', 'coarse_pathologic_stage', 'age_groups')

    ix_rm = c() ; var_vec_ex = c()
    for(v_v in var_vec) {
        all_clin[v_v][,1] = as.character(all_clin[v_v][,1])
        ix_na = which(is.na(all_clin[v_v][,1]))
        if (dim(all_clin)[1] - length(ix_na) < 30) { var_vec_ex = c(var_vec_ex, v_v) }
        if (dim(all_clin)[1] - length(ix_na) >= 30) {
            ix_rm = union(ix_rm, which(is.na(all_clin[v_v][,1]))) #; print(ix_rm) ; print(v_v)
        }
        t_ = table(as.character(all_clin[v_v][,1])) ; cat_n = names(which(t_ <= 10))
        if ((length(t_) - length(cat_n)) > 1 & length(cat_n) >= 1) {
            for(n_ in cat_n) { ix_rm = union(ix_rm, which(all_clin[v_v][,1] == n_)) }
        }
        if ((length(t_) - length(cat_n)) <= 1) { var_vec_ex = c(var_vec_ex, v_v) }
    }
    var_vec = setdiff(var_vec, var_vec_ex) ; all_clin = all_clin[setdiff(seq(dim(all_clin)[1]), ix_rm), ]
    pl_hist = nullGrob()
    if (length(var_vec) > 0) {
    dis_hist = all_clin[,c(var_vec, 'MPS_groups')]
    if (sum(table(dis_hist$MPS_groups) > 20) == length(unique(dis_hist$MPS_groups))) {
        pl_hist = list()
        for(h_v in sort(var_vec)) {
            if (h_v != 'PSA_value') {
                tmp_h = dis_hist[,union('MPS_groups', h_v)] ; colnames(tmp_h)[which(colnames(tmp_h) == h_v)] = 'var'
                t_ = table(as.character(tmp_h$var)) ; tmp_h$var = as.character(tmp_h$var)
                if (length(t_) > 1) {
                    for(i in names(t_)) { ix_ = which(tmp_h$var == i) ; tmp_h$var[ix_] = paste(i, ' (', as.numeric(t_[i]), ')',sep='') }
                    pl_ = ggplot(tmp_h, aes(var, fill=MPS_groups)) + geom_bar(stat = "count",position = "fill", width=0.25) + 
                            scale_y_continuous(labels=scales::percent) + theme_bw() + xlab('') + #coord_flip() +
                            scale_fill_manual(values=col_surv) + ggtitle(gsub('_', ' ', h_v)) + ylab('percent') + coord_flip()#+ facet_wrap(~ hist_var, scales = "free_x")
                    pl_hist[[h_v]] = pl_
                }
            }
            if (h_v == 'PSA_value') {
                tmp_h = dis_hist[,union('MPS_groups', h_v)] ; colnames(tmp_h)[which(colnames(tmp_h) == h_v)] = 'var'
                pl_ = ggplot(tmp_h, aes(MPS_groups, log2(as.numeric(var)), fill=MPS_groups)) + theme_bw() + xlab('') + ylab(paste('log2(', gsub('_', ' ', h_v), ') a.u')) + 
                        geom_violin(alpha=0.5, draw_quantiles=c(0.25, 0.5, 0.75)) + scale_fill_manual(values=col_surv) + 
                        ggtitle(gsub('_', ' ', h_v)) + coord_flip()
                pl_hist[[h_v]] = pl_
            }
        }
    }}
    val_ret = pl_hist
}

plot_Histopath_multivariate <- function(clin_data) {
    all_clin = clin_data ; rownames(all_clin) = all_clin$sample_analyzed ; all_clin = all_clin[! is.na(all_clin$MPS_groups),]
    dis_ = as.character(unique(all_clin$disease_type)) ; path_id = as.character(unique(all_clin$module_id))
    var_vec = c('MPS_groups', get_commonHistopath(dis_)) #c('histological_type', 'coarse_pathologic_stage', 'age_groups')

    ix_rm = c() ; var_vec_ex = c()
    for(v_v in var_vec) {
        all_clin[v_v][,1] = as.character(all_clin[v_v][,1])
        ix_na = which(is.na(all_clin[v_v][,1]))
        if (dim(all_clin)[1] - length(ix_na) < 30) { var_vec_ex = c(var_vec_ex, v_v) }
        if (dim(all_clin)[1] - length(ix_na) >= 30) {
            ix_rm = union(ix_rm, which(is.na(all_clin[v_v][,1]))) #; print(ix_rm) ; print(v_v)
        }
        t_ = table(as.character(all_clin[v_v][,1])) ; cat_n = names(which(t_ <= 10))
        if ((length(t_) - length(cat_n)) > 1 & length(cat_n) >= 1) {
            for(n_ in cat_n) { ix_rm = union(ix_rm, which(all_clin[v_v][,1] == n_)) }
        }
        if ((length(t_) - length(cat_n)) <= 1) { var_vec_ex = c(var_vec_ex, v_v) }
    }
    var_vec = setdiff(var_vec, var_vec_ex) ; all_clin = all_clin[setdiff(seq(dim(all_clin)[1]), ix_rm), ]
    p_for_ovs = p_for_pfs = nullGrob()
    if (length(setdiff(var_vec,'MPS_groups')) > 0) {
    #all_clin$mod_censor_time = all_clin$mod_censor_time/30
    all_clin_num = all_clin ; for(v_v in var_vec) { all_clin_num[v_v][,1] = as.numeric(factor(all_clin_num[v_v][,1])) }
    ss_ovs = Surv(as.numeric(as.character(all_clin$OS.time)),all_clin$OS)
    ss_pfs = Surv(as.numeric(as.character(all_clin$PFI.time)),all_clin$PFI)
    val_ret = list(NA, NA, NA) ; p_for_ovs = p_for_pfs = nullGrob()
    if (! (sum(is.na(ss_ovs)) == length(is.na(ss_ovs)))) {
        su_ovs = Surv((all_clin$OS.time), (all_clin$OS))
        s_cox_ovs = coxph(as.formula(paste('su_ovs ~', paste(var_vec, collapse='+'))), data=all_clin_num)
        p_for_ovs = ggforest(s_cox_ovs,all_clin_num, main='OVS',refLabel = "reference") + theme_bw()
        
    }
    if (! (sum(is.na(ss_pfs)) == length(is.na(ss_pfs)))) {
        su_pfs = Surv((all_clin$PFI.time), (all_clin$PFI))
        s_cox_pfs = coxph(as.formula(paste('su_pfs ~', paste(var_vec, collapse='+'))), data=all_clin_num)
        p_for_pfs = ggforest(s_cox_pfs,all_clin_num, main='PFS',refLabel = "reference") + theme_bw()
    }}
    val_ret = list(p_for_ovs, p_for_pfs)
}

plot_Histopath_subset <- function(clin_data) {#1 = 'group', var2 = 'histology') {
    all_clin = clin_data ; rownames(all_clin) = all_clin$sample_analyzed ; all_clin = all_clin[! is.na(all_clin$MPS_groups),]
    dis_ = as.character(unique(all_clin$disease_type)) ; path_id = as.character(unique(all_clin$module_id))
    var_vec = c('MPS_groups', get_commonHistopath(dis_)) #c('histological_type', 'coarse_pathologic_stage', 'age_groups')

    ix_rm = c() ; var_vec_ex = c()
    for(v_v in var_vec) {
        all_clin[v_v][,1] = as.character(all_clin[v_v][,1])
        ix_na = which(is.na(all_clin[v_v][,1]))
        if (dim(all_clin)[1] - length(ix_na) < 30) { var_vec_ex = c(var_vec_ex, v_v) }
        if (dim(all_clin)[1] - length(ix_na) >= 30) {
            ix_rm = union(ix_rm, which(is.na(all_clin[v_v][,1]))) #; print(ix_rm) ; print(v_v)
        }
        t_ = table(as.character(all_clin[v_v][,1])) ; cat_n = names(which(t_ <= 10))
        if ((length(t_) - length(cat_n)) > 1 & length(cat_n) >= 1) {
            for(n_ in cat_n) { ix_rm = union(ix_rm, which(all_clin[v_v][,1] == n_)) }
        }
        if ((length(t_) - length(cat_n)) <= 1) { var_vec_ex = c(var_vec_ex, v_v) }
    }
    var_vec = setdiff(var_vec, var_vec_ex) ; all_clin = all_clin[setdiff(seq(dim(all_clin)[1]), ix_rm), ]
    c_ = 1 ; pl_list_ovs = pl_list_pfs = list()
    if (length(setdiff(var_vec, 'MPS_groups')) > 0) {
    	for(v_v in setdiff(var_vec, 'MPS_groups')) {
        	var_cat = as.character(unique(all_clin[v_v][,1]))
	        for(v_c in var_cat) {
    	        print(v_c)
    	        all_clin_var = all_clin[which(all_clin[v_v][,1] == v_c),]
        	    pl_o = pl_p = nullGrob()
            	if ((sum(table(all_clin_var$MPS_groups) > 5) >= 2) & (sum(table(all_clin_var$MPS_groups) > 5) == length(unique(all_clin_var$MPS_groups)))) {
                	surv_ovs = getSurv(all_clin_var, 'OS',1) ; surv_pfs = getSurv(all_clin_var, 'PFI',1)
	                m_ovs= paste('p=', formatC(as.numeric(surv_ovs[[11]]$pval_surv),format='e',digits=1), 
    	                     ' | HR=', signif(as.numeric(surv_ovs[[11]]$HZR_cox), 2))
        	        m_sub_ovs = paste(v_v, ' (', v_c,')', sep='')
            		
            		m_pfs= paste('p=', formatC(as.numeric(surv_pfs[[11]]$pval_surv),format='e',digits=1), 
                    	     ' | HR=', signif(as.numeric(surv_pfs[[11]]$HZR_cox), 2))
	                m_sub_pfs = m_sub_ovs
        			legend_labs = c('MPS-', 'MPS+')
    	            pl_ovs = ggsurvplot(surv_ovs[[2]], surv_ovs[[3]], pval=F, title=m_ovs, font.title=9, censor.shape=124,censor.size=1,
            	                legend.labs = legend_labs,xlab= 'Time (months)', ylab='Survival fraction (OVS)',
                	            subtitle=m_sub_ovs, font.subtitle=7,surv.median.line='hv', palette =col_surv, risk.table=F)
    				pl_pfs = ggsurvplot(surv_pfs[[2]], surv_pfs[[3]], pval=F, title=m_pfs, font.title=9, censor.shape=124,censor.size=1,
                    	        legend.labs = legend_labs, xlab= 'Time (months)', ylab='Survival fraction (PFS)',
                        	    subtitle=m_sub_pfs, font.subtitle=7,surv.median.line='hv', palette =col_surv, risk.table=F)
    				pl_o = pl_ovs$plot ; pl_p = pl_pfs$plot ; pl_list_ovs[[c_]] = pl_o ; pl_list_pfs[[c_]] = pl_p ; c_ = c_ + 1
            	}
	        }
    	}
    }
    val_ret = list(pl_list_ovs, pl_list_pfs)         
}

get_MPS_existingModule <- function(disease_name = 'pancreatic', module_identifier = 'MSigDBONC_52', number_rand = 10000, MPS_thresh = 0) {
	t1 = Sys.time()
	dis_ = disease_name ; path_id = module_identifier
	mod_info = readRDS(paste('data/module_info/', path_id,'_rds', sep=''))
	g_set = as.character(readRDS('data/common_genes_TCGA_rds')[,1])
	inp_g = unlist(strsplit(as.character(mod_info$genes_with_alias), '\\|'))
	genes_in_mod = intersect(inp_g, g_set) ; div_f = as.numeric(length(genes_in_mod))
	genes_in_mod_vec = paste(unique(as.character(sort(genes_in_mod))), collapse='|')
    set.seed(108)
    num_rand = number_rand
    mod_name = paste(as.character(mod_info$pathway_name), collapse=',')
    out_mod_dir = paste('results/', path_id,'/',sep='')
    system(gsub('DIR', out_mod_dir, cmd_mkdir))
    saveRDS(genes_in_mod, paste(out_mod_dir, 'binary_input_genes_in_module', sep=''))
#    saveRDS(v_rand, paste(out_mod_dir, 'binary_null_distribution_', num_rand, sep=''))
    write.table(data.frame(gene = genes_in_mod), paste(out_mod_dir,'input_genes_in_module',sep=''),sep='\t', quote=F,row.names=F)
#    write.table(data.frame(null_dist = v_rand), paste(out_mod_dir,'null_distribution_module',sep=''),sep='\t', quote=F,row.names=F)
	dis_clin = readRDS(paste('data/disease_clinical/',dis_, '_clin_rds',sep=''))
	mps_inp = readRDS(paste('data/disease_mps/', dis_, '/', path_id,'mps1000_rds',sep=''))
	v_z = as.numeric(mps_inp[path_id][,1]) ; s_set = rownames(mps_inp)
	tmp_ = data.frame(sample_analyzed = s_set, MPS = v_z) ; rownames(tmp_) = s_set
	tmp_ = merge(tmp_, dis_clin) ; tmp_$module_name = mod_name
    tmp_$module_id = path_id ; tmp_$genes_in_module = genes_in_mod_vec
    tmp_$num_genes_in_module = div_f
    tmp_$MPS_groups = NA ; 
    tmp_$MPS_groups[which(tmp_$MPS > MPS_thresh)] = 'MPS+'
    tmp_$MPS_groups[which(tmp_$MPS < (-1*MPS_thresh))] = 'MPS-'
    saveRDS(tmp_, paste(out_mod_dir, 'binary_', dis_,'_MPS', sep=''))
    write.table(tmp_, paste(out_mod_dir,dis_,'_MPS',sep=''), sep='\t', quote=F, row.names=F)
    t2 = Sys.time() ; print(paste('done reading MPS for ', path_id, 'in', dis_, difftime(t2,t1)))
    all_clin = tmp_ ; all_clin = all_clin[! is.na(all_clin$MPS_groups),]
    all_clin$OS.time = all_clin$OS.time/30 ; all_clin$PFI.time = all_clin$PFI.time/30
    plot_all_plots(all_clin, out_mod_dir)
}

get_MPS_newModule <- function(path_to_input, select_cohorts = '', module_name = '', number_bins = 10, pval_thresh = 0.01, MPS_significance = 0, number_rand = 10000, MPS_thresh = 0) {
    if (select_cohorts == 'full') { coh_list = names(tcga_key) }
    if (select_cohorts != 'full') { coh_list = intersect(select_cohorts, names(tcga_key)) }
    
    #cont_genes = disc_genes = data.frame(gene = as.character(gene2module$Approved.Symbol))
    T_1 = Sys.time()
    input_genes = read.csv(path_to_input, sep='\t', header=T)
    inp_g = as.character(input_genes[,1])
    g_set = as.character(readRDS('data/common_genes_TCGA_rds')[,1])
    #inp_g = as.character(read.csv('~/Desktop/shiny_test/test_input_genes',sep='\t',header=T)[,1])
    genes_in_mod = intersect(inp_g, g_set) ; div_f = as.numeric(length(genes_in_mod))
    genes_in_mod_vec = paste(unique(as.character(sort(genes_in_mod))), collapse='|')
    set.seed(108)
    num_rand = number_rand
    total_genes = length(g_set)
    total_genes_bin = rmultinom(1, total_genes, rep((1/number_bins), number_bins))
    num_in_path = div_f
    path_genes_bin = rmultinom(num_rand, num_in_path, rep((1/number_bins), number_bins))
    c_vec_r = path_genes_bin ; nc_vec_r = total_genes_bin[,1] - c_vec_r
    v_rand = sapply(seq(num_rand), function(i) calculateMI_v2(c_vec_r[,i], nc_vec_r[,i]))
    mu_r = mean(v_rand) ; sd_r = sd(v_rand)
    v_rand_q = as.numeric(quantile(v_rand, (1-pval_thresh)))
    len_rand = length(v_rand)

    mod_name = module_name
    if (module_name == '') { mod_name = rev(unlist(strsplit(path_to_input,'\\/')))[1] }
    #paste(module_name,'OUTPUT_',gsub(':','.',gsub(' ','_at_', Sys.time())),sep='')
    out_mod_dir = paste('results/', mod_name,'/',sep='')
    system(gsub('DIR', out_mod_dir, cmd_mkdir))
    saveRDS(genes_in_mod, paste(out_mod_dir, 'binary_input_genes_in_module', sep=''))
    saveRDS(v_rand, paste(out_mod_dir, 'binary_null_distribution_', num_rand, sep=''))
    write.table(data.frame(gene = genes_in_mod), paste(out_mod_dir,'input_genes_in_module',sep=''),sep='\t', quote=F,row.names=F)
    write.table(data.frame(null_dist = v_rand), paste(out_mod_dir,'null_distribution_module',sep=''),sep='\t', quote=F,row.names=F)
    
    for(dis_ in coh_list) {
        con_g = as.data.frame(vroom(paste('data/disease_gene_expression/', dis_, '_tumdat_zscore', sep=''), delim = "\t",col_names=T))
        dis_g = as.data.frame(vroom(paste('data/disease_gene_expression/', dis_, '_tumdata_zscore_bins10', sep=''), delim = "\t",col_names=T))
        rownames(con_g) = as.character(con_g$gene) ; con_g = con_g[,-grep('gene', colnames(con_g))]
        rownames(dis_g) = as.character(dis_g$gene) ; dis_g = dis_g[,-grep('gene', colnames(dis_g))]
        dis_clin = readRDS(paste('data/disease_clinical/',dis_, '_clin_rds',sep=''))
        g_row = intersect(g_set, intersect(rownames(con_g), rownames(dis_g)))
        s_set = intersect(colnames(con_g), colnames(dis_g))
        con_g = con_g[g_row, s_set] ; dis_g = dis_g[g_row, s_set]
        
        Nrow = dim(con_g)[1] ; per_bin = round(Nrow/number_bins)
        t1 = Sys.time()
        registerDoMC(cores=4)
        out_list <- foreach(sam_ = s_set) %dopar% {
        #for(sam_ in s_set) {
            vv = con_g[sam_] ; vv$bin = 0 ; vv[genes_in_mod,]$bin = 1 ; sign_fact = sign(cor(vv)[1,2])
            vv = dis_g[sam_]
            l_ = lapply(seq(number_bins), function(i) rownames(vv)[vv[sam_][,1] == i])
            c_vec = unlist(lapply(lapply(l_, intersect, genes_in_mod), length))
            nc_vec= as.numeric(sapply(l_, length) - c_vec )
            mi_ = calculateMI_v2(c_vec, nc_vec) ; mi_sign = sign_fact*mi_
        }
        mi_vec = unlist(out_list) ; abs_mi = abs(mi_vec) ; ix_0 = which(abs_mi < v_rand_q)
        v_z = (abs_mi - mu_r)/sd_r ; v_z[v_z < 0] = 0 ; v_z[is.na(v_z)] = 0
        if (! MPS_significance) { v_z[ix_0] = 0 } ## MPS_significance == 0 will NOT impose the p-value threshold (flag 'pval_thresh') for MPS
        v_z = v_z*sign(mi_vec)
        tmp_ = data.frame(sample_analyzed = s_set, MPS = v_z) ; rownames(tmp_) = s_set
        tmp_ = merge(tmp_, dis_clin) ; tmp_$module_name = mod_name
        tmp_$module_id = paste('user_module', mod_name) ; tmp_$genes_in_module = genes_in_mod_vec
        tmp_$num_genes_in_module = div_f
        tmp_$MPS_groups = NA ; 
        tmp_$MPS_groups[which(tmp_$MPS > MPS_thresh)] = 'MPS+'
        tmp_$MPS_groups[which(tmp_$MPS < (-1*MPS_thresh))] = 'MPS-'
        saveRDS(tmp_, paste(out_mod_dir, 'binary_', dis_,'_MPS', sep=''))
        write.table(tmp_, paste(out_mod_dir,dis_,'_MPS',sep=''), sep='\t', quote=F, row.names=F)
        t2 = Sys.time() ; print(paste('done MPS calculation in', dis_, difftime(t2,t1)))
        all_clin = tmp_ ; all_clin = all_clin[! is.na(all_clin$MPS_groups),]
        all_clin$OS.time = all_clin$OS.time/30 ; all_clin$PFI.time = all_clin$PFI.time/30
        plot_all_plots(all_clin, out_mod_dir)
    }
    T_2 = Sys.time() ; print(paste('done ', difftime(T_2,T_1)))
}

plot_all_plots <- function(clin_data, output_dir) {
	out_mod_dir = output_dir
	all_clin = clin_data
	all_clin = all_clin[! is.na(all_clin$MPS_groups),]
    dis_ = as.character(unique(all_clin$disease_type))
    path_id = as.character(unique(all_clin$module_id))
    name_module = unique(as.character(all_clin$module_name))
    name_module_trim = strtrim(gsub('at','-',gsub('OUTPUT','',gsub('-','',gsub('_', '', unique(as.character(all_clin$module_name)))))),17)
    genes_in_module = sort(unique(unlist(strsplit(as.character(all_clin$genes_in_module), '\\|'))))
    num_genes_in_module = unique(as.numeric(all_clin$num_genes_in_module))
	
	out_dir = paste(out_mod_dir, 'plots_', dis_,'/', sep='')# = paste('results/', mod_name,'/',sep='')
    system(gsub('DIR', out_dir, cmd_mkdir))

    pdf(paste(out_dir,'plot_survival.pdf', sep=''), width=5, height=5) ; print(plot_Surv(all_clin)) ; dev.off()

	pdf(paste(out_dir,'plot_volcano.pdf', sep=''), width=4, height=4) ; print(plot_geneExp(get_geneExp_info(all_clin))) ; dev.off()

	pdf(paste(out_dir,'plot_histopath.pdf', sep=''), width=6, height=4) ; print(plot_Histopath(all_clin)) ; dev.off()

	pdf(paste(out_dir,'plot_histopath_multivariate.pdf', sep=''), width=7, height=5) ; print(plot_Histopath_multivariate(all_clin)) ; dev.off()

	pdf(paste(out_dir,'plot_histopath_subset.pdf', sep=''), width=4, height=4) ; print(plot_Histopath_subset(all_clin)) ; dev.off()

	pdf(paste(out_dir,'plot_immune.pdf', sep=''), width=5, height=3) ; print(plot_Immune(all_clin)) ; dev.off()

}

