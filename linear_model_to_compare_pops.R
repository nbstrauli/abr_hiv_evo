check_for_sig_ab_artifacts = function(ab_stat, sample_stats_filepath, covariance_plots_dirpath){
	#ab_stat is the name of the antibody repertoire statisitic to test for correlation
	t = read.table(sample_stats_filepath, header=TRUE)
	t$patient_id = as.factor(t$patient_id)
	t$ab_barcode_id = as.factor(t$ab_barcode_id)
	#remove 4th time-point from patient 3, outlier for Ab data
	t_complete = subset(t, sample_ID!='3_4')

	artifact_data = data.frame(patient_id.=t_complete$patient_id, time_point.=t_complete$time_point, abr_library_creation_date.=t_complete$abr_library_creation_date, abr_rna_extraction_date.=t_complete$abr_rna_extraction_date, ab_sequencing_date.=t_complete$ab_sequencing_date, hiv_rna_extraction_date.=t_complete$hiv_rna_extraction_date, hiv_cdna_synthesis_date.=t_complete$hiv_cdna_synthesis_date, hiv_1st_round_pcr_date.=t_complete$hiv_1st_round_pcr_date, hiv_2nd_round_pcr_date.=t_complete$hiv_2nd_round_pcr_date, abr_initial_library_concentration.=t_complete$abr_initial_library_concentration, abr_rna_sample_concentration.=t_complete$abr_rna_sample_concentration, hiv_read_depth.=t_complete$hiv_read_depth, hiv_sequencing_date.=t_complete$hiv_sequencing_date, hiv_initial_library_concentration.=t_complete$hiv_initial_seq_lib_concentration)#add a '.' to end of coefficient names for easy parsing later on
	ab_value = t_complete[[ab_stat]]
	ab_artifact_data = data.frame(artifact_data, ab_value)

	ab_artifacts = c('abr_initial_library_concentration.', 'abr_rna_sample_concentration.')
	for (i in ab_artifacts){
		single_artifact_data = data.frame(patient_id=artifact_data$patient_id., ab_artifact=artifact_data[[i]], ab_value=ab_value, time_point=artifact_data$time_point., abr_rna_sample_concentration=artifact_data$abr_rna_sample_concentration.)
		single_artifact_data = single_artifact_data[complete.cases(single_artifact_data),]
		#if the experimental artifact variable is continuous then don't model is as a random effect (duh)
		if ((i == 'abr_initial_library_concentration.')|(i == 'abr_rna_sample_concentration.')){
			single_ab_artifact_model = lmer(ab_value ~ 1 + (1|patient_id) + ab_artifact, data=single_artifact_data, REML=FALSE)
		} else{
			single_ab_artifact_model = lmer(ab_value ~ 1 + (1|patient_id) + (1|ab_artifact), data=single_artifact_data, REML=FALSE)
		}
		null_model = lmer(ab_value ~ (1|patient_id), data=single_artifact_data, REML=FALSE)
		a = anova(null_model, single_ab_artifact_model)
		LRT_pval = a[2,8]
		if (LRT_pval < 0.05){
			print ('#############################')
			print (paste('AbR artifact ', i, ' model', sep=''))
			print (LRT_pval)
			print ('#############################')
			print ('')
			#plot Ab covariance plot
			library(GGally)
			pm = ggpairs(single_artifact_data, columns=c(1,2,3,4,5), cardinality_threshold=20)
			# p = getPlot(pm, 2, 2)
			# p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
			dir.create(file.path(covariance_plots_dirpath, 'abr', ab_stat), showWarnings=FALSE)
			pdf(paste(covariance_plots_dirpath, 'abr/', ab_stat, '/', i, '_plot.pdf', sep=''))
			print (pm)
			dev.off()
		}
	}
}

check_for_sig_hiv_artifacts = function(hiv_stat, sample_stats_filepath, covariance_plots_dirpath){
	#hiv_stat is the name of the HIV population statisitic to test for correlation
	t = read.table(sample_stats_filepath, header=TRUE)
	t$patient_id = as.factor(t$patient_id)
	t$ab_barcode_id = as.factor(t$ab_barcode_id)
	t_complete = t
	#remove outlier viral load sample, if testing hiv 'viral_load'
	if (hiv_stat == 'viral_load'){
		t_complete = subset(t_complete, sample_ID!='7_1')
	}

	artifact_data = data.frame(patient_id.=t_complete$patient_id, time_point.=t_complete$time_point, abr_library_creation_date.=t_complete$abr_library_creation_date, abr_rna_extraction_date.=t_complete$abr_rna_extraction_date, ab_sequencing_date.=t_complete$ab_sequencing_date, hiv_rna_extraction_date.=t_complete$hiv_rna_extraction_date, hiv_cdna_synthesis_date.=t_complete$hiv_cdna_synthesis_date, hiv_1st_round_pcr_date.=t_complete$hiv_1st_round_pcr_date, hiv_2nd_round_pcr_date.=t_complete$hiv_2nd_round_pcr_date, abr_initial_library_concentration.=t_complete$abr_initial_library_concentration, abr_rna_sample_concentration.=t_complete$abr_rna_sample_concentration, hiv_read_depth.=t_complete$hiv_read_depth, hiv_sequencing_date.=t_complete$hiv_sequencing_date, hiv_initial_library_concentration.=t_complete$hiv_initial_seq_lib_concentration)#add a '.' to end of coefficient names for easy parsing later on
	hiv_value = t_complete[[hiv_stat]]
	hiv_artifact_data = data.frame(artifact_data, hiv_value)

	#test individual HIV artifacts for signifigance
	hiv_artifacts = c('hiv_read_depth.', 'hiv_initial_library_concentration.')
	for (i in hiv_artifacts){
		single_artifact_data = data.frame(patient_id=artifact_data$patient_id., hiv_artifact=artifact_data[[i]], hiv_value=hiv_value, time_point=artifact_data$time_point.)
		single_artifact_data = single_artifact_data[complete.cases(single_artifact_data),]
		#if the experimental artifact variable is continuous then don't model is as a random effect (duh)
		if ((i == 'hiv_initial_library_concentration.')|(i == 'hiv_read_depth.')){
			single_hiv_artifact_model = lmer(hiv_value ~ 1 + (1|patient_id) + hiv_artifact, data=single_artifact_data, REML=FALSE)
		} else{
			single_hiv_artifact_model = lmer(hiv_value ~ 1 + (1|patient_id) + (1|hiv_artifact), data=single_artifact_data, REML=FALSE)
		}
		null_model = lmer(hiv_value ~ 1 + (1|patient_id), data=single_artifact_data, REML=FALSE)
		a = anova(null_model, single_hiv_artifact_model)
		LRT_pval = a[2,8]
		if (LRT_pval < 0.05){
			print ('#############################')
			print (paste('HIV artifact ', i, ' model', sep=''))
			print (LRT_pval)
			print ('#############################')
			print ('')
			#plot Ab covariance plot
			library(GGally)
			pm = ggpairs(single_artifact_data, columns=c(1,2,3,4), cardinality_threshold=20)
			# p = getPlot(pm, 2, 2)
			# p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
			dir.create(file.path(covariance_plots_dirpath, 'hiv', hiv_stat), showWarnings=FALSE)
			pdf(paste(covariance_plots_dirpath, 'hiv/', hiv_stat, '/', i, '_plot.pdf', sep=''))
			print (pm)
			dev.off()
		}
	}
}

fit_model = function(ab_stat, hiv_stat, sample_stats_filepath, heatmaps_dirpath, covariance_plots_dirpath){
	#ab_stat is the name of the antibody repertoire statisitic to test for correlation
	#hiv_stat is the name of the HIV population statisitic to test for correlation

	#### pltting parameters ####
	line_width = 1
	point_size = 2
	plot_height = 2.25
	plot_width = 3
	############################

	t = read.table(sample_stats_filepath, header=TRUE)
	t$patient_id = as.factor(t$patient_id)
	t$ab_barcode_id = as.factor(t$ab_barcode_id)
	#remove 4th time-point from patient 3, outlier for Ab data
	t_complete = subset(t, sample_ID!='3_4')
	#remove outlier viral load sample, if testing hiv 'viral_load'
	if (hiv_stat == 'viral_load'){
		t_complete = subset(t_complete, sample_ID!='7_1')
	}

	artifact_data = data.frame(patient_id.=t_complete$patient_id, time_point.=t_complete$time_point, abr_library_creation_date.=t_complete$abr_library_creation_date, abr_rna_extraction_date.=t_complete$abr_rna_extraction_date, ab_sequencing_date.=t_complete$ab_sequencing_date, hiv_rna_extraction_date.=t_complete$hiv_rna_extraction_date, hiv_cdna_synthesis_date.=t_complete$hiv_cdna_synthesis_date, hiv_1st_round_pcr_date.=t_complete$hiv_1st_round_pcr_date, hiv_2nd_round_pcr_date.=t_complete$hiv_2nd_round_pcr_date, abr_initial_library_concentration.=t_complete$abr_initial_library_concentration, abr_rna_sample_concentration.=t_complete$abr_rna_sample_concentration, hiv_read_depth.=t_complete$hiv_read_depth, hiv_sequencing_date.=t_complete$hiv_sequencing_date, hiv_initial_library_concentration.=t_complete$hiv_initial_seq_lib_concentration)#add a '.' to end of coefficient names for easy parsing later on
	hiv_value = t_complete[[hiv_stat]]
	ab_value = t_complete[[ab_stat]]

	final_data = data.frame(ab_value, hiv_value, patient_id=artifact_data$patient_id., abr_initial_library_concentration=artifact_data$abr_initial_library_concentration., abr_rna_sample_concentration=artifact_data$abr_rna_sample_concentration.)
	final_data = final_data[complete.cases(final_data),]

	#fit model
	final_model = lmer(ab_value ~ hiv_value + (1|patient_id), data=final_data, REML=FALSE)
	final_null_model = lmer(ab_value ~ 1 + (1|patient_id), data=final_data, REML=FALSE)

	final_model_slope = coef(summary(final_model))[2,1]
	a = anova(final_null_model, final_model)
	final_LRT_pval = a[2,8]
	if (final_LRT_pval < 0.05){
		print ('#############################')
		print (paste('Final lin. model:', ab_stat, 'vs.', hiv_stat, sep=' '))
		print (final_LRT_pval)
		print ('#############################')
	}

	#if the final model p value is significant then make a scatter plot b/t the ab values and hiv values, except first regress out the effect of patient_id (and potentially any artifcat variables) on the ab_values.
	if (final_LRT_pval <= 0.05){
		ab_artifact_model = lmer(ab_value ~ 1 + (1|patient_id), data=final_data, REML=FALSE)
		hiv_artifact_model = lmer(hiv_value ~ 1 + (1|patient_id), data=final_data, REML=FALSE)
		#assign unique color foreach patient
		patients_sorted = sort(as.numeric(levels(final_data$patient_id)))
		num_patients = length(patients_sorted)
		colors = rainbow(num_patients * 1.3)[1:num_patients]
		pat_col_list = list()
		index = 0
		for (i in patients_sorted){
			index = index + 1
			pat_col_list[[i]] = colors[index]
		}
		cols = c()
		for (i in final_data$patient_id){
			cols = c(cols, pat_col_list[[as.numeric(i)]])
		}

		plot_data = data.frame(HIV_values=resid(hiv_artifact_model), AbR_values=resid(ab_artifact_model), patients=final_data$patient_id)
		p = ggplot(data=plot_data, aes(x=HIV_values, y=AbR_values))
		p = p + geom_abline(intercept=0, slope=final_model_slope, linetype=2, size=line_width)
		p = p + geom_point(aes(colour=patients), size=point_size)
		p = p + theme_bw()
		p = p + theme(axis.title=element_blank(), legend.key.size=unit(0.4, "cm"))
		pdf(paste(heatmaps_dirpath, ab_stat, '_vs_', hiv_stat, '_plot.pdf', sep=''), height=plot_height, width=plot_width)
		plot(p)
		dev.off()

		# pdf(paste(heatmaps_dirpath, ab_stat, '_vs_', hiv_stat, '_plot.pdf', sep=''))
		# plot(c(), xlim=range(resid(hiv_artifact_model)), ylim=range(resid(ab_artifact_model)), main=paste('Residual of ', hiv_stat, ' vs ', ab_stat, sep=''), xlab='Residual of HIV Stat', ylab='Residual of AbR Stat')
		# abline(coef=c(0,final_model_slope), lty=2, lwd=line_width)
		# points(resid(hiv_artifact_model), resid(ab_artifact_model), col=cols, pch=point_type, cex=point_size)
		# legend('bottomright', legend=patients_sorted, col=colors, pch=point_type, pt.cex=point_size, title='Patient ID')
		# dev.off()

		#as a sanity check that the above scatter plot with the line is indeed a good way to visualize the relationship b/t the ab_value and hiv_value, we did a simple linear regression on the residuals (that are plotted above) to make sure that the line that results from this is (very) similar to the one that is plotted above. If one were to plot the line from 'lin_mod' below, they should see that it is almost identical to the line plotted above.
		df = data.frame(ab=resid(ab_artifact_model), hiv=resid(hiv_artifact_model))
		lin_mod = lm(ab ~ hiv, data=df)
		# print (coef(lin_mod))
	}

	return(final_LRT_pval)
}

test_for_cors = function(sample_stats_filepath, heatmaps_dirpath, covariance_plots_dirpath){
	library(lme4)
	ab_sample_stats_to_test = c('ab_diversity_pi', 'abr_baseline_mean_sigma_CDR', 'abr_baseline_mean_sigma_FWR', 'ab_divergence')
	hiv_sample_stats_to_test = c('viral_load', 'hiv_diversity_pi', 'hiv_synonymous_divergence', 'hiv_non_synonymous_divergence', 'hiv_selection_dN_dS')

	for (ab_stat in ab_sample_stats_to_test){
		check_for_sig_ab_artifacts(ab_stat=ab_stat, sample_stats_filepath=sample_stats_filepath, covariance_plots_dirpath=covariance_plots_dirpath)
	}

	for (hiv_stat in hiv_sample_stats_to_test){
		check_for_sig_hiv_artifacts(hiv_stat=hiv_stat, sample_stats_filepath=sample_stats_filepath, covariance_plots_dirpath=covariance_plots_dirpath)
	}

	#One must visually check to see the output of each of these artifact models to see if any of the artifact varaiables were indeed significant. If they are, then manually include them in the final model below. 

	names = list(ab_sample_stats_to_test, hiv_sample_stats_to_test)
	pvals = matrix(nrow=length(ab_sample_stats_to_test), ncol=length(hiv_sample_stats_to_test), dimnames=names)
	pvals_raw = matrix(nrow=length(ab_sample_stats_to_test), ncol=length(hiv_sample_stats_to_test), dimnames=names)
	for (ab_stat_index in 1:length(ab_sample_stats_to_test)){
		for (hiv_stat_index in 1:length(hiv_sample_stats_to_test)){
			print ('#########################')
			print ('#########################')
			print (paste('Ab stat: ', ab_sample_stats_to_test[ab_stat_index], '; HIV stat: ', hiv_sample_stats_to_test[hiv_stat_index]), sep='')
			print ('#########################')
			print ('#########################')
			print ('')
			print ('')
			print ('')
			print ('')
			print ('')
			pval = fit_model(ab_sample_stats_to_test[ab_stat_index], hiv_sample_stats_to_test[hiv_stat_index], sample_stats_filepath, heatmaps_dirpath, covariance_plots_dirpath)
			pvals_raw[ab_stat_index, hiv_stat_index] = pval
			pvals[ab_stat_index, hiv_stat_index] = -log10(min(c(1, (pval * length(pvals)))))#bonfr. correction and then negative log transform p vals
		}
	}
	library(gplots)
	print (pvals_raw)
	print ('')
	print (pvals)
	pdf(paste(heatmaps_dirpath, 'ab_stats_vs_hiv_stats.pdf'))
	ab_names = c('pi', 'selection\n(CDR)', 'selection\n(FWR)', 'divergence')
	hiv_names = c('viral load', 'pi', 'syn.\ndivergence', 'non-syn.\ndivergence', 'dN/dS')
	heatmap.2(pvals, trace='none', margin=c(10,10), labRow=ab_names, labCol=hiv_names, denscol='black', key.title=NA, key.xlab='-log10(p value)', xlab='HIV Summary Statistics', ylab='AbR Summary Statistics')
	dev.off()
	print (pvals)
}

#fit_model(ab_stat='ab_baseline_selection_CDR', hiv_stat='hiv_selection_dN_dS', sample_stats_filepath='/Users/nstrauli/data/abr_hiv_coevo/sample_stats/changeo/sample_stats.txt')
# test_for_cors(sample_stats_filepath='/Users/nstrauli/data/abr_hiv_coevo/sample_stats/changeo/sample_stats.txt', heatmaps_dirpath='/Users/nstrauli/data/abr_hiv_coevo/sample_stats/changeo/sum_stat_corrs/', covariance_plots_dirpath='/Users/nstrauli/data/abr_hiv_coevo/sample_stats/changeo/covariance_plots/')
