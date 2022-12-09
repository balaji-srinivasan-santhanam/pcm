library('randomForestSRC')
library('caret')
library('survival')
library('gridExtra')
library('grid')

cohort = 'pancreatic' ; d_t = 'paad'


getSurv_rfsrc <- function(all_clin_df, surv_type = 'OS', rand_iter = 1000, samp_name = 'high_risk', ctrl_name = 'low_risk') {
    all_clin = all_clin_df
    s_name = samp_name ; c_name = ctrl_name
    all_clin$OS.time = all_clin$OS.time/30 ; all_clin$PFI.time = all_clin$PFI.time/30
    if (surv_type == 'OS') { all_clin$mod_censor_time = as.numeric(all_clin$OS.time) ; all_clin$death_event_binary = as.numeric(all_clin$OS) }
    if (surv_type =='PFI') { all_clin$mod_censor_time =as.numeric(all_clin$PFI.time) ; all_clin$death_event_binary = as.numeric(all_clin$PFI) }
    len_s = sum(all_clin$risk_group == samp_name) ; len_c = sum(all_clin$risk_group == ctrl_name)
    s_all = survfit(Surv(as.numeric(all_clin$mod_censor_time), all_clin$death_event_binary) ~ 1)
    ss = Surv(as.numeric(as.character(all_clin$mod_censor_time)),all_clin$death_event_binary)
    val_ret = NA
    pv = s_fit = pv_cox = hz_ratio = s_coxph = z = c_ix = NA
    if (! (sum(is.na(ss)) == length(is.na(ss)))) {
        s_fit = survfit(Surv(mod_censor_time, death_event_binary) ~ risk_group, data = all_clin)
        s_diff = tryCatch({ 
            survdiff(Surv(mod_censor_time, death_event_binary) ~ risk_group, data = all_clin)
            }, error = function(err) { list(chisq=NA, n=dim(all_clin)[1]) } )
        pv <- ifelse(is.na(s_diff),1,(round(1 - pchisq(s_diff$chisq, length(s_diff$n) - 1),50)))[[1]]
        
        s_coxph = coxph(Surv(mod_censor_time, death_event_binary) ~ risk_group, data = all_clin)
        if (! is.na(s_coxph[['coefficients']])) { 
            s_c = summary(s_coxph)
            hz_ratio = as.numeric(data.frame(s_c$coefficients)['exp.coef.'])
            pv_cox = as.numeric(data.frame(s_c$coefficients)['Pr...z..'])
            z = as.numeric(data.frame(s_c$coefficients)['z'])
            c_ix = as.numeric(summary(s_coxph)$concordance['C'])
        }

        pv_r = c()
        for(r_i in seq(rand_iter)) {
            all_clin_r = all_clin ; prob_v = as.numeric(prop.table(table(all_clin$risk_group)))
            all_clin_r$risk_group = sample(c(s_name, c_name), dim(all_clin)[1], replace=T, prob=prob_v)
            s1_r = tryCatch({ 
            survdiff(Surv(mod_censor_time, death_event_binary) ~ risk_group, data = all_clin_r)
            }, error = function(err) { list(chisq=NA, n=dim(all_clin_r)[1]) } )
            pv_r = c(pv_r, round(1 - pchisq(s1_r$chisq, length(s1_r$n) - 1),50))
        }
        fp = (sum(pv_r < pv)/rand_iter)
        median_surv = summary(s_fit)$table[,'median']
        med_ab = round(as.numeric(median_surv[paste('risk_group=',s_name,sep='')]),1)
        med_re = round(as.numeric(median_surv[paste('risk_group=',c_name,sep='')]),1)
        m_sub = paste(s_name,' ', len_s, '(', med_ab,' mo)',
                                    '| ', c_name,' ', len_c, '(', med_re, ' mo)', sep='')
    
        out_df = data.frame(pval_surv = pv, pval_cox = pv_cox, C_index = c_ix, z_cox = z,HZR_cox = hz_ratio, info_survival = m_sub)
        val_ret = list(KM_pv = pv, KM_fit = s_fit, clin_data = all_clin, cox_pv = pv_cox, hzr = hz_ratio, 
                                cox_fit = s_coxph, avg_fit = s_all, rand_pv = pv_r, z=z,
                                concordance = c_ix, out_df = out_df, fdr = fp)
    }
}


rand_int = sample(unique(as.integer(floor(runif(15, 50, 50000)))), 10)
surv_type_list = list('OS' = 'ovs', 'PFI' = 'pfs')
#surv_type = 'OS' ; surv_ = surv_type_list[[surv_type]]
surv_type = 'PFI' ; surv_ = surv_type_list[[surv_type]]

folds = 10 ; repeat_iterations = 10 ; rm_list_dis = list()
rf_inp_dir = 'data/rsf/'
rf_out_dir = 'results_rsf/'
cmd_mkdir = "mkdir -p DIR"
system(gsub("DIR", paste(rf_out_dir, cohort, '/',sep=''), cmd_mkdir))


immune_df = read.csv(paste(rf_inp_dir, 'TCGA.Kallisto.fullIDs.cibersort.relative.tsv', sep=''), sep='\t', header=T)
immune_df$sample_analyzed = gsub('\\.', '-', as.character(immune_df$SampleID))
immune_df = immune_df[,-grep('SampleID|CancerType|P.value|Correlation|RMSE', colnames(immune_df))]
immune_df = immune_df[which(! duplicated(immune_df$sample_analyzed)),]
rownames(immune_df) = as.character(immune_df$sample_analyzed) ; immune_df = immune_df[,-grep('sample',colnames(immune_df))]
colnames(immune_df) = paste('IMMUNE_', colnames(immune_df), sep='')
immune_df$sample_analyzed = rownames(immune_df)
col_immune = colnames(immune_df)[grep('IMMUNE_', colnames(immune_df))]

dis_clin = readRDS(paste('data/disease_clinical/',cohort, '_clin_rds',sep=''))
dis_clin$status = dis_clin[surv_type][,1] ; dis_clin$time = dis_clin[paste(surv_type,'.time',sep='')][,1]/30
dis_clin_time = dis_clin[,c('sample_analyzed', 'time', 'status')]

if (surv_type == 'OS') { pid_list = readRDS(paste(rf_inp_dir, 'PCMs_paad_ovs',sep='')) }
if (surv_type == 'PFI') { pid_list = readRDS(paste(rf_inp_dir, 'PCMs_paad_pfs',sep='')) }

dis_mps = readRDS(paste(rf_inp_dir, 'MPS_paad', sep=''))
dis_snv = readRDS(paste(rf_inp_dir, 'prominent_SNV_paad', sep=''))
c_n_snv = colnames(dis_snv)[grep('SNV_', colnames(dis_snv))]
dis_cna = readRDS(paste(rf_inp_dir, 'prominent_CNA_paad', sep=''))
c_n_cna = colnames(dis_cna)[grep('CNA_', colnames(dis_cna))]

d_snv = merge(dis_clin[,c('sample_analyzed', 'sample_TCGA')], dis_snv, all.x=T)
d_cna = merge(dis_clin[,c('sample_analyzed', 'sample_TCGA')], dis_cna, all.x=T)
d_mps = dis_mps ; d_mps$sample_analyzed = as.character(rownames(d_mps))

d_snv = d_snv[,-grep('sample_TCGA', colnames(d_snv))]
d_cna = d_cna[,-grep('sample_TCGA', colnames(d_cna))]
for(v_v in c_n_snv) { d_snv[v_v][,1] = factor(d_snv[v_v][,1]) }
for(v_v in c_n_cna) { d_cna[v_v][,1] = factor(d_cna[v_v][,1]) }
d_gen = merge(d_snv, d_cna, all.x=T)

hist_v = c('histological_type', 'ajcc_pathologic_tumor_stage','coarse_pathologic_stage', 'age_groups')
dis_hist = dis_clin[,c('sample_analyzed', hist_v)]

if (cohort == 'breast') {
	hist_v = c('estrogen_receptor_status','progesterone_receptor_status','TNBC_status',
	                   'her2_receptor_status','histological_type', 'ajcc_pathologic_tumor_stage',
	                   'coarse_pathologic_stage', 'age_groups')
	dis_hist = dis_clin[,c('sample_analyzed', hist_v)]
}
if (cohort== 'prostate') {
    hist_v = c('gleason_score','PSA_value','histological_type', 'ajcc_pathologic_tumor_stage',
	                   'coarse_pathologic_stage', 'age_groups')
    dis_hist = dis_clin[,c('sample_analyzed', hist_v)]
}
if (cohort == 'AML') {
    hist_v = c('mol_test_status','histological_type', 'ajcc_pathologic_tumor_stage',
	                   'coarse_pathologic_stage', 'age_groups')
    dis_hist = dis_clin[,c('sample_analyzed', hist_v)]
}
if (cohort == 'head_neck') {
    hist_v = c('HPV_status','histological_type', 'ajcc_pathologic_tumor_stage',
	                   'coarse_pathologic_stage', 'age_groups')
	dis_hist = dis_clin[,c('sample_analyzed', hist_v)]
}
if (cohort == 'cervical') {
	hist_v = c('HPV_status','histological_type', 'ajcc_pathologic_tumor_stage',
	                   'coarse_pathologic_stage', 'age_groups')
	dis_hist = dis_clin[,c('sample_analyzed', hist_v)]
}
if (cohort == 'colon') {
	hist_v = c('MSI_status','tumor_side', 'histological_type', 'ajcc_pathologic_tumor_stage',
	                   'coarse_pathologic_stage', 'age_groups')
	dis_hist = dis_clin[,c('sample_analyzed', hist_v)]
}

for(v_v in hist_v) { dis_hist[v_v][,1] = factor(dis_hist[v_v][,1]) }
rm_v = c()
t_ = d_gen[,-grep('sample',colnames(d_gen))]
if (sum(rowSums(is.na(t_)) != dim(t_)[2]) < 0.25*dim(t_)[1]) { c_n = c() ; rm_v = c(rm_v, 'gen') } ; rm(t_)
t_ = dis_hist[,-grep('sample',colnames(dis_hist))]
if (sum(rowSums(is.na(t_)) != dim(t_)[2]) < 0.25*dim(t_)[1]) { hist_v = c() ; rm_v = c(rm_v, 'std') } ; rm(t_)
rm_list_dis[[cohort]] = rm_v

d_mpsonly = merge(dis_clin_time, d_mps, all.x=T)# ; rownames(d_mpsonly) = as.character(d_mpsonly$sample_analyzed)
d_std = merge(dis_clin_time, dis_hist, all.x=T)#; rownames(d_std) = as.character(d_std$sample_analyzed)
d_stdgen = merge(merge(dis_clin_time, dis_hist, all.x = T), d_gen, all.x=T)
d_all = merge(merge(d_mpsonly, d_stdgen, all.x=T), immune_df, all.x=T)
d_all = d_all[! duplicated(d_all$sample_analyzed),]
rownames(d_all) = as.character(d_all$sample_analyzed)
var_mpsonly = c('time', 'status', pid_list)
var_std = c('time', 'status', hist_v)
var_snv = c('time', 'status', c_n_snv)
var_cna = c('time', 'status', c_n_cna)
var_gen = c('time', 'status', c_n_cna, c_n_snv)
var_stdsnv = c('time', 'status', hist_v, c_n_snv)
var_stdcna = c('time', 'status', hist_v, c_n_cna)
var_stdgen = c('time', 'status', hist_v, c_n_cna, c_n_snv)
var_immune = c('time', 'status', col_immune)
var_std_immune = c('time', 'status', hist_v, col_immune)
var_std_immune_mps = c('time', 'status', hist_v, col_immune, pid_list)
var_stdsnv_mps = c('time', 'status', hist_v, c_n_snv, pid_list)
var_stdcna_mps = c('time', 'status', hist_v, c_n_cna, pid_list)
var_stdgen_mps = c('time', 'status', hist_v, c_n_cna, c_n_snv, pid_list)
var_std_cna_immune_mps = c('time', 'status', hist_v, c_n_cna, col_immune, pid_list)
var_cnamps = c('time', 'status', c_n_cna, pid_list)
var_snvmps = c('time', 'status', c_n_snv, pid_list)
var_immunemps = c('time', 'status', col_immune, pid_list)
var_stdmps = c('time', 'status', hist_v, pid_list)


var_list = list('_MPSonly' = var_mpsonly, '_standard' = var_std,'_std_MPS' = var_stdmps)
iter_ = 1
t1 = Sys.time()
while( iter_ <= repeat_iterations) {
	set.seed(rand_int[iter_])
    cvIndex <- createFolds(factor(d_all$status), folds, returnTrain = T)
    saveRDS(cvIndex, paste(rf_out_dir, cohort, '/',d_t,'_indexFile_list_iter', iter_, '_of_', repeat_iterations, '_CV', folds,sep=''))
    fi_cv = paste(rf_out_dir, cohort, '/',d_t,'_indexFile_list_iter', iter_, '_of_', repeat_iterations, '_CV', folds,sep='')
   	cvIndex = readRDS(fi_cv)

   	for(var_name in names(var_list)) {
	    col_n = var_list[[var_name]]
	    d = d_all[,col_n] ; rf_out_list = list()
	    for(i in seq(folds)) {
	        train_data <- d[cvIndex[[i]],] ; eval_data <- d[-cvIndex[[i]],]
	        rf_tr = tryCatch({ rf_tr = rfsrc(Surv(time, status) ~ ., train_data, mtry = ceiling(log(dim(d)[2]))) }, error = function(err) { list(NA) } )
	        rf_pred = list(NA) ; if (! is.na(rf_tr)) { rf_pred = predict.rfsrc(rf_tr, eval_data, mtry = ceiling(log(dim(d)[2]))) }
	        rf_out_list[['train']][[i]] = rf_tr ; rf_out_list[['prediction']][[i]] = rf_pred ; rm(rf_tr) ; rm(rf_pred)
	    }
	    print(paste('done RF ', 'iteration number ', iter_, ' of ', repeat_iterations, ' for ', var_name, ' with ', folds,'-fold cross-validation in ', cohort, sep=''))
	    fi_ = paste(rf_out_dir, cohort, '/', d_t,'_', surv_type, var_name,'_iter', iter_, '_of_', repeat_iterations, '_CV',folds,'_2020-05-07',sep='')
	    saveRDS(rf_out_list, fi_) ; rm(rf_out_list)
	}
	iter_ = iter_ + 1
}
t2 = Sys.time()

frac_vec = c(1, 0.666, 0.5, 0.4) ; q_v = c(25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)
covariate_categ = names(var_list)
surv_out_df = c()
len_f_vec = c()
for(var_ in covariate_categ) {
   	iter_ = sample(seq(repeat_iterations))[1] ; i = sample(seq(folds))[1] ; 
    tmp_ = readRDS(paste(rf_out_dir, cohort,'/', d_t, "_", surv_type, var_, "_iter", iter_, "_of_10_CV10_2020-05-07", sep=''))
    tr_ = tmp_[['train']] ; pr_ = tmp_[['prediction']]
	len_f_vec = c(len_f_vec, (dim(tr_[[i]]$yvar)[1] + dim(pr_[[i]]$yvar)[1]))
}
len_f = min(len_f_vec) ; rm(iter_) ; rm(i) ; rm(tr_) ; rm(pr_)
for(var_ in covariate_categ) {
	for(iter_ in seq(repeat_iterations)) {
		tmp_ = readRDS(paste(rf_out_dir, cohort,'/', d_t, "_", surv_type, var_, "_iter", iter_, "_of_10_CV10_2020-05-07", sep=''))
        tr_ = tmp_[['train']] ; pr_ = tmp_[['prediction']]
	    mort_v = c() ; samp_v = c()
    	for(i in seq(folds)) { mort_v = c(mort_v, pr_[[i]]$predicted) ; samp_v = c(samp_v, rownames(pr_[[i]]$xvar)) }
        mort_v_ = mort_v[! is.na(mort_v)] ; samp_v_ = samp_v[! is.na(mort_v)]
        if (length(mort_v_) > 0) {
        	df_ = data.frame(sample_analyzed = samp_v_, mortality = mort_v_)
	        s_dec = as.character(df_[sort(df_$mortality, decreasing=T, index.return=T)$ix,]$sample_analyzed)
    	    s_inc = as.character(df_[sort(df_$mortality, decreasing=F, index.return=T)$ix,]$sample_analyzed)
            for(frac_ in frac_vec) {
        	    len_ = frac_*len_f ; q_ = q_v[which.min(abs(floor(len_ / 50)*25 - q_v))]
				samp_high_risk = s_dec[1:q_] ; samp_low_risk = s_inc[1:q_]
				clin_df = rbind(data.frame(sample_analyzed = samp_high_risk, risk_group = rep('high_risk', length(samp_high_risk))),
								data.frame(sample_analyzed = samp_low_risk, risk_group = rep('low_risk', length(samp_low_risk))))
				clin_df = merge(clin_df, dis_clin)
                s_ = getSurv_rfsrc(clin_df, surv_type, rand_iter = 1000, samp_name = 'high_risk', ctrl_name = 'low_risk')
                su_df = s_[[11]] ; su_df$cohort = cohort ; su_df$iteration = iter_ ; su_df$q = q_ ; su_df$frac = frac_
            	su_df$var_category = var_
	            surv_out_df = rbind(surv_out_df, su_df)
    	        print(paste('done: ', cohort,   iter_, ' of 10 in ', var_, q_))
        	}
	    }
    }
}
print(paste('done: ', surv_type, cohort))
ab_v = sapply(strsplit(as.character(surv_out_df$info_survival), '\\|'), '[', 1)
be_v = sapply(strsplit(as.character(surv_out_df$info_survival), '\\|'), '[', 2)
surv_out_df$num_ab = as.numeric(gsub('high_risk ', '', sapply(strsplit(ab_v, '\\('), '[', 1)))
surv_out_df$num_be = as.numeric(gsub('low_risk ', '', sapply(strsplit(be_v, '\\('), '[', 1)))
write.table(surv_out_df, paste(rf_out_dir, 'survival_rf_', cohort, '_', surv_type, '_',Sys.Date(), sep=''), sep='\t', quote=F, row.names=F)
pdf(paste(rf_out_dir, 'survival_rf_', cohort, '_', surv_type, '_',Sys.Date(),'.pdf', sep=''),width=10, height=7)
for(f_ in unique(surv_out_df$frac)) {
	df_plot = surv_out_df[which(surv_out_df$frac == f_), ] ; q_ = as.numeric(unique(df_plot$num_ab))
	df_plot$hazard_ratio = 2^abs(log2(df_plot$HZR_cox))
	title_text = paste('n = ', q_, ' ( frac = ', f_, ' )',sep='')
	p1 = ggplot(df_plot, aes(cohort, hazard_ratio, col=var_category)) +
		geom_jitter(width=0.05,alpha=0.5, size=2) + 
		coord_flip() + ggtitle(title_text) + theme_bw() + ylab('hazard ratio') + xlab('')
	p2 = ggplot(df_plot, aes(cohort, hazard_ratio, col=var_category)) +
		geom_jitter(width=0.05,alpha=0.5, size=2) +
		coord_flip() + ggtitle(title_text) + theme_bw() + ylab('hazard ratio') + xlab('')
	grid.arrange(p1, p2, nullGrob(), nullGrob())
}
t3 = Sys.time()

print('done')


