import os
import sys

#This translates the patient IDs in the options database to the patient IDs that we gave.
options_id_dic = {'305':'1', '377':'2', '667':'3', '682':'4', '1313':'5', '1439':'6', '1465':'7', '1686':'8', '1721':'9', '1723':'10'}

def gather_stats_changeo_processed_data(sample_names_key_filepath, cd4_count_viral_load_filepath, hiv_read_depth_dirpath, hiv_pi_diversity_dirpath, ab_library_conc_filepath, ab_rna_conc_filepath, ab_pi_diversity_dirpath, samp_seqDate_and_barcodeID_filepath, irep_barcode_seqs_filepath, selection_abr_baseline_dirpath, divergence_abr_dirpath, divergence_hiv_dirpath, selection_hiv_dirpath, selection_abr_baseline_mean_dirpath, hiv_sample_prep_dirpath, hiv_seq_run_info_filepath, art_status_filepath, output_filepath):
    """
    This script takes in various inputs, parses them, then outputs them into a large (ggplot friendly) table.
    'sample_names_key_filepath' = Path to the file that has a key that maps all the different sample names to each other. This also has the time-point (in days since infection) information for each sample. The path for this should be '/Users/nstrauli/data/abr_hiv_coevo/sample_names_key/key.txt'
    'cd4_count_viral_load_filepath' = Path to the file that contains the CD4 count and viral load measurements for each sample. Originally this path was '/Users/nstrauli/data/options_cohort_data/CD4_viral_load/cd4_viral_load_data'
    'hiv_read_depth_dirpath' = This is the path to the directory that contians the read depth for each file. The path to this directory should be: '/Users/nstrauli/data/abr_hiv_coevo/read_depth_foreach_sample/hiv'.
    'hiv_pi_diversity_dirpath' = The path to the directory that contains the diversity, in pi, for each of the HIV samples.
    'ab_library_conc_filepath' = This is the path to the file in the sample database that records the info for all the AbR libraries for each sample. The path for this file should be '/Users/nstrauli/data/options_cohort_sample_processing/sample_database/PBMCs/seq_libraries.txt'.
    'ab_rna_conc_filepath' = This is the path to the file in the sample database that records the info for each of the RNA extration samples. It has the RNA concentration measurements. THe path to the file should be '/Users/nstrauli/data/options_cohort_sample_processing/sample_database/PBMCs/rna_extractions.txt'
    'ab_pi_diversity_dirpath' = This is the path to the directory that contains the pi diversity values for each of the AbR samples. The path for this should be '/Users/nstrauli/data/abr_hiv_coevo/diversity/abr/pi'
    'samp_seqDate_and_barcodeID_filepath' = This is the path to the file that contains the BCR sequencing date as well as barcode ID for each of the PBMC samples in the sample set. The path for this should be '/Users/nstrauli/data/abr_hiv_coevo/sample_seq_date_and_barcode_id/sample_seq_date_and_barcode_id.txt'
    'irep_barcode_seqs_filepath' = This is the path to the file that contains the actual nucleotide sequences for each of the iRepertoire barcode IDs (i.e. maps barcode ID to a barcode sequence). The path for this file should be '/Users/nstrauli/data/abr_hiv_coevo/abr_barcode_analysis/irepertoire_barcod_seqs.txt'
    'selection_abr_baseline_dirpath' = This is the path to the directory that contains the overall selection values for each abr sample.
    'divergence_abr_dirpath' = This is the path to the directory that contains the divergence information as calculated from the changeo processed data.
    'divergence_hiv_dirpath' = This is the path to the directory that contains the divergence information for the HIV data
    'selection_hiv_dirpath' = This is the path to the directory that has mean dN/dS values for the seqs in the data.
    'selection_abr_baseline_mean_dirpath' = This is the path to the directory that contains the mean sigma selection values for each of the samples in the data. It is similar to 'selection_abr_baseline_dirpath', but instead of using BASELINE's convolution technique to get one summary stat for the sample, we took the mean (weighting by seq count) of each of the sequences baseline sigma values (for CDR and FWR regions).
    'hiv_seq_run_info_filepath' - This is the filepath that contains the sequencing date as well as the original library concentration for each sample.
    'output_filepath' = path to the output file.
    """
    #this will be a dictionary, where each index is a sample, and is
    #recorded as 'patient_id'_'time-point'. So, day 275 of patient 7
    #would be encoded as '7_275'
    sample_stats_dic = {}
    filein = open(sample_names_key_filepath, "r")
    filein.readline()
    #this will be a dic that maps a time-point rank (or order) to
    #that time-points day (as in days after infection). So if we
    #encounter a sample that is day 254 from patient 4 then we can 
    #look at this dictionary to see that that day corresponds to
    #the 4th time-point for patient 4
    tpoint_rank_dic = {}
    for i in filein:
        line = i.split('\t')
        patient_id = line[0].split('_')[0]
        tpoint = line[0].split('_')[1]
        tpoint_rank = line[1].split('_')[2]
        sample_id = '%s_%s' % (patient_id, tpoint_rank)
        #currently, the ordering of each of the elements in these lists
        #is as follows:
        #0: sample ID
        #1: CD4 count
        #2: viral load
        #3: hiv sequencing depth
        #4: hiv diversity statistic, pi
        #5: ab initial library concentration
        #6: ab initial library concentration, duplicate
        #7: ab library creation date
        #8: ab library creation date, duplicate
        #9: ab RNA sample concentration
        #10: ab RNA sample concentration, duplicate
        #11: ab RNA extration date
        #12: ab RNA extration date, duplicate
        #13: ab diversity statistic, pi
        #14: ab diversity statistic, pi, duplicate
        #15: ab barcode ID
        #16: ab barcode ID, duplicate
        #17: ab barcode sequence
        #18: ab barcode sequence, duplicate
        #19: ab sequencing date
        #20: ab sequencing date, duplicate
        #21: ab 'baseline' selection value, CDR
        #22: ab 'baseline' selection value, CDR lower confidence interval
        #23: ab 'baseline' selection value, CDR upper confidence interval
        #24: ab 'baseline' selection value, CDR p value
        #25: ab 'baseline' selection value, FWR
        #26: ab 'baseline' selection value, FWR lower confidence interval
        #27: ab 'baseline' selection value, FWR upper confidence interval
        #28: ab 'baseline' selection value, FWR p value
        #29: ab 'baseline' selection value, CDR, duplicate
        #30: ab 'baseline' selection value, CDR lower confidence interval, duplicate
        #31: ab 'baseline' selection value, CDR upper confidence interval, duplicate
        #32: ab 'baseline' selection value, CDR p value, duplicate
        #33: ab 'baseline' selection value, FWR, duplicate
        #34: ab 'baseline' selection value, FWR lower confidence interval, duplicate
        #35: ab 'baseline' selection value, FWR upper confidence interval, duplicate
        #36: ab 'baseline' selection value, FWR p value, duplicate
        #37: ab mean divergence
        #38: ab mean divergence, duplicate
        #39: hiv mean synonymous divergence
        #40: hiv mean non-synonymous divergence
        #41: hiv mean dN/dS
        #42: abr mean baseline sigma, CDR
        #43: abr mean baseline sigma, CDR, duplicate
        #44: abr mean baseline sigma, FWR
        #45: abr mean baseline sigma, FWR, duplicate
        #46: hiv RNA extraction date
        #47: hiv cDNA synthesis date
        #48: hiv 1st round PCR date
        #49: hiv 2nd round PCR date
        #50: hiv sequencing date
        #51: hiv initial seq library concentration
        #52: ART status
        sample_stats_dic['%s_%s' % (patient_id, tpoint)] = [sample_id, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']
        #fill in tpoint rank dic
        tpoint_rank_dic[sample_id] = tpoint
    filein.close()
    #fill in the descriptor for each stat of a sample. These will
    #be the header line in the output file (in the order filled in)
    header = ['patient_id', 'time_point', 'sample_ID']

    #get CD4 count and viral load data
    filein = open(cd4_count_viral_load_filepath, "r")
    filein.readline()
    for i in filein:
        line = i[:-1].split('\t')
        patient_id = options_id_dic[line[2].split('_')[1]]
        tpoint = line[6]
        cd4_count = line[7]
        viral_load = line[8]
        try:
            sample_stats_dic['%s_%s' % (patient_id, tpoint)][1] = cd4_count
            sample_stats_dic['%s_%s' % (patient_id, tpoint)][2] = viral_load
        except KeyError:
            #There should be 2 samples found in this data that were 
            #not found in the time-point dataset because they are the
            #two samples that did not get shipped, and were thus removed
            #from the time-point dataset. That samples are patient 1, day 4183
            #and patient 5, day 219
            print 'These samples were found in the CD4 count/viral load data, but not in the time-point data:'
            print 'patient ID:', patient_id
            print 'time-point:', tpoint
    filein.close()
    #update the header
    header += ['CD4_count', 'viral_load']

    #get HIV sequencing depth info
    if hiv_read_depth_dirpath[-1] != '/':
        hiv_read_depth_dirpath += '/'
    for i in os.listdir(hiv_read_depth_dirpath):
        if i[0] == '.' or i[:6] == 'README':
            continue
        input_filepath = hiv_read_depth_dirpath + i
        filein = open(input_filepath, "r")
        h = filein.readline()
        for j in filein:
            line = j[:-1].split('\t')
            tpoint_id = line[0]
            tpoint = tpoint_rank_dic[i + '_' + tpoint_id]
            read_depth = line[1]
            sample_stats_dic['%s_%s' % (i, tpoint)][3] = read_depth
        filein.close()
    #update header
    header.append('hiv_read_depth')

    #get HIV pi diversity stats
    ##### diversity specific parameters #####
    #this is the number of downsampling trials
    #used when calculating diversity
    diversity_calc_trials = 10
    ##### diversity specific parameters #####
    if hiv_pi_diversity_dirpath[-1] != '/':
        hiv_pi_diversity_dirpath += '/'
    for i in os.listdir(hiv_pi_diversity_dirpath):
        if i[0] == '.' or i[:6] == 'README':
            continue
        patient_id = i
        for j in os.listdir(hiv_pi_diversity_dirpath + i):
            if j[0] == '.' or j[:6] == 'README':
                continue
            tpoint = j
            input_filepath = hiv_pi_diversity_dirpath + i + '/' + j
            filein = open(input_filepath, "r")
            filein.readline()
            diversities = []
            trials = []
            for k in filein:
                line = k[:-1].split('\t')
                trial = line[0]
                trials.append(trial)
                diversity = float(line[1])
                diversities.append(diversity)
            filein.close()
            #check if the mean has already been calculated
            if 'mean' in trials:
                diversity = str(diversities[trials.index('mean')])
                #check to see if all trials are present. There should be 10 trials
                if len(diversities) != diversity_calc_trials+1:
                    print 'strange number of diversity values in:'
                    print 'patient:', patient_id
                    print 'time-point:', tpoint
            else:
                diversity = str(sum(diversities) / len(diversities))
                #check to see if all trials are present. There should be 10 trials
                if len(diversities) != diversity_calc_trials:
                    print 'strange number of diversity values in:'
                    print 'patient:', patient_id
                    print 'time-point:', tpoint
            sample_stats_dic['%s_%s' % (patient_id, tpoint)][4] = diversity
    #update header
    header.append('hiv_diversity_pi')

    #get AbR library concentration measurements
    filein = open(ab_library_conc_filepath, "r")
    filein.readline()
    for i in filein:
        line = i.split('\t')
        sample_id = line[0]
        sample_info = sample_id.split('_')
        patient_id = sample_info[1]
        tpoint_id = sample_info[2]
        tpoint = tpoint_rank_dic['%s_%s' % (patient_id, tpoint_id)]
        concentration = line[4]
        date = line[3]
        if ';' in concentration:
            #these are the samples that have duplicates
            #Need to make sure to fill in approp. column
            #for the duplicates
            if sample_id == '1_2_17' or sample_id == '1_3_8' or sample_id == '1_10_10':
                concentration_dup = concentration.split(';')[1]
                date_dup = date.split(';')[1]
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][6] = concentration_dup
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][8] = date_dup
                concentration = concentration.split(';')[0]
                date = date.split(';')[0]
            #these are the samples that failed to
            #amplify the first time, so we only take the
            #2nd value
            else:
                concentration = concentration.split(';')[1]
                date = date.split(';')[1]
        elif concentration == '':
            continue
        sample_stats_dic['%s_%s' % (patient_id, tpoint)][5] = concentration
        if len(date.split('/')[-1]) == 2:#if the full year isin't written, add it
            date = '%s20%s' % (date[:-2], date[-2:])
        sample_stats_dic['%s_%s' % (patient_id, tpoint)][7] = date
    filein.close()
    #update header
    header += ['abr_initial_library_concentration', 'abr_initial_library_concentration_duplicate', 'abr_library_creation_date', 'abr_library_creation_date_duplicate']

    #get AbR RNA concentration after extraction
    filein = open(ab_rna_conc_filepath, "r")
    filein.readline()
    for i in filein:
        line = i.split('\t')
        sample_id = line[0]
        sample_info = sample_id.split('_')
        patient_id = sample_info[1]
        tpoint_id = sample_info[2]
        tpoint = tpoint_rank_dic['%s_%s' % (patient_id, tpoint_id)]
        concentration = line[3]
        date = line[2]
        if ';' in concentration:
            #these are the samples that have duplicates.
            #Need to make sure to fill in approp. column
            #for the duplicates
            if sample_id == '1_2_17' or sample_id == '1_3_8' or sample_id == '1_10_10':
                concentration_dup = concentration.split(';')[1]
                date_dup = date.split(';')[1]
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][10] = concentration_dup
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][12] = date_dup
                concentration = concentration.split(';')[0]
                date = date.split(';')[0]
            #these are the samples that failed to
            #amplify the first time, so we take the
            #2nd value
            else:
                concentration = concentration.split(';')[1]
                date = date.split(';')[1]
        elif concentration == '':
            continue
        if len(date.split('/')[-1]) == 2:#if the full year isin't written, add it
            date = '%s20%s' % (date[:-2], date[-2:])
        sample_stats_dic['%s_%s' % (patient_id, tpoint)][9] = concentration
        sample_stats_dic['%s_%s' % (patient_id, tpoint)][11] = date
    filein.close()
    #update header
    header += ['abr_rna_sample_concentration', 'abr_rna_sample_concentration_duplicate', 'abr_rna_extraction_date', 'abr_rna_extraction_date_duplicate']

    #get AbR pi stats foreach sample
    if ab_pi_diversity_dirpath[-1] != '/':
        ab_pi_diversity_dirpath += '/'
    for i in os.listdir(ab_pi_diversity_dirpath):
        if i[0] == '.' or i[:6] == 'README' or i == 'unknown':
            continue
        patient_id = i
        for j in os.listdir(ab_pi_diversity_dirpath + i):
            if j[0] == '.' or j[:6] == 'README':
                continue
            #get pi
            filein = open(ab_pi_diversity_dirpath+i+'/'+j, "r")
            filein.readline()
            for k in filein:
                line = k[:-1].split('\t')
                if line[0] == 'mean':
                    pi = line[1]
                    break
            filein.close()
            tpoint = j
            #if t point is a duplicate
            if '.' in tpoint:
                tpoint = tpoint.split('.')[0]
                #then assign pi to 14th column
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][14] = pi
            else:
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][13] = pi
    #update header
    header += ['ab_diversity_pi', 'ab_diversity_pi_duplicate']

    #get the Ab sequencing barcode ID, barcode sequence, and Ab sequencing date
    filein = open(irep_barcode_seqs_filepath, 'r')
    filein.readline()
    #this will be a dic that maps barcode ID to barcode seq
    irep_barcode_dic = {}
    for i in filein:
        line = i[:-1].split('\t')
        irep_barcode_dic[line[0]] = line[1]
    filein.close()
    filein = open(samp_seqDate_and_barcodeID_filepath, 'r')
    filein.readline()
    for i in filein:
        replicate_indicator = False
        line = i[:-1].split('\t')
        samp_id = line[0].split('_')
        pat_id = samp_id[1]
        tpoint_id = samp_id[2]
        #if this is a duplicate
        if '(' in tpoint_id:
            #then fix time-point id and set the indicator variable
            tpoint_id = tpoint_id.split('(')[0]
            replicate_indicator = True
        tpoint = tpoint_rank_dic['%s_%s' %(pat_id, tpoint_id)]
        seq_date = line[1]
        barcode_id = line[2]
        barcode_seq = irep_barcode_dic[barcode_id]
        if replicate_indicator:
            sample_stats_dic['%s_%s' % (pat_id, tpoint)][16] = barcode_id
            sample_stats_dic['%s_%s' % (pat_id, tpoint)][18] = barcode_seq
            sample_stats_dic['%s_%s' % (pat_id, tpoint)][20] = seq_date
        else:
            sample_stats_dic['%s_%s' % (pat_id, tpoint)][15] = barcode_id
            sample_stats_dic['%s_%s' % (pat_id, tpoint)][17] = barcode_seq
            sample_stats_dic['%s_%s' % (pat_id, tpoint)][19] = seq_date
    filein.close()
    #update header
    header += ['ab_barcode_id', 'ab_barcode_id_duplicate', 'ab_barcode_seq', 'ab_barcode_seq_duplicate', 'ab_sequencing_date', 'ab_sequencing_date_duplicate']

    #get baseline selection values
    if selection_abr_baseline_dirpath[-1] != '/':
        selection_abr_baseline_dirpath += '/'
    for i in os.listdir(selection_abr_baseline_dirpath):
        if i[0] == '.' or i[:6] == 'README' or i == 'unknown' or i[-4:] == '.pdf':
            continue
        patient_id = i
        for j in os.listdir(selection_abr_baseline_dirpath + i):
            if j[0] == '.' or j[:6] == 'README' or j[-4:] == '.pdf':
                continue
            tpoint = j
            filein = open(selection_abr_baseline_dirpath + i + '/' + j, 'r')
            filein.readline()
            selection_sum_stats = {}
            for k in filein:
                line = k[:-1].split('\t')
                region = line[1]
                sigma = line[2]
                ci_low = line[3]
                ci_high = line[4]
                p_val = line[5]
                selection_sum_stats[region + '_sigma'] = sigma
                selection_sum_stats[region + '_ci_low'] = ci_low
                selection_sum_stats[region + '_ci_high'] = ci_high
                selection_sum_stats[region + '_p_val'] = p_val
            filein.close()
            if not '.' in tpoint:
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][21] = selection_sum_stats['CDR_sigma']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][22] = selection_sum_stats['CDR_ci_low']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][23] = selection_sum_stats['CDR_ci_high']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][24] = selection_sum_stats['CDR_p_val']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][25] = selection_sum_stats['FWR_sigma']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][26] = selection_sum_stats['FWR_ci_low']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][27] = selection_sum_stats['FWR_ci_high']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][28] = selection_sum_stats['FWR_p_val']
            else:
                tpoint = tpoint.split('.')[0]
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][29] = selection_sum_stats['CDR_sigma']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][30] = selection_sum_stats['CDR_ci_low']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][31] = selection_sum_stats['CDR_ci_high']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][32] = selection_sum_stats['CDR_p_val']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][33] = selection_sum_stats['FWR_sigma']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][34] = selection_sum_stats['FWR_ci_low']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][35] = selection_sum_stats['FWR_ci_high']
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][36] = selection_sum_stats['FWR_p_val']
    #update header
    header += ['ab_baseline_selection_CDR', 'ab_baseline_selection_CDR_confidence_int_low', 'ab_baseline_selection_CDR_confidence_int_high', 'ab_baseline_selection_CDR_p_value', 'ab_baseline_selection_FWR', 'ab_baseline_selection_FWR_confidence_int_low', 'ab_baseline_selection_FWR_confidence_int_high', 'ab_baseline_selection_FWR_p_value', 'ab_baseline_selection_CDR_duplicate', 'ab_baseline_selection_CDR_confidence_int_low_duplicate', 'ab_baseline_selection_CDR_confidence_int_high_duplicate', 'ab_baseline_selection_CDR_p_value_duplicate', 'ab_baseline_selection_FWR_duplicate', 'ab_baseline_selection_FWR_confidence_int_low_duplicate', 'ab_baseline_selection_FWR_confidence_int_high_duplicate', 'ab_baseline_selection_FWR_p_value_duplicate']

    #get abr divergence values
    if divergence_abr_dirpath[-1] != '/':
        divergence_abr_dirpath += '/'
    for i in os.listdir(divergence_abr_dirpath):
        if i[0] == '.' or i[:6] == 'README' or i[-4:] == '.pdf' or i == 'unknown':
            continue
        patient_id = i
        filein = open(divergence_abr_dirpath + i, "r")
        filein.readline()
        for j in filein:
            line = j[:-1].split('\t')
            tpoint = line[0]
            divergence = line[1]
            if '.' in tpoint:
                tpoint = tpoint.split('.')[0]
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][38] = divergence
            else:
                sample_stats_dic['%s_%s' % (patient_id, tpoint)][37] = divergence
        filein.close()
    #update header
    header += ['ab_divergence', 'ab_divergence_duplicate']

    #get hiv divergence values
    if divergence_hiv_dirpath[-1] != '/':
        divergence_hiv_dirpath += '/'
    for i in os.listdir(divergence_hiv_dirpath):
        if i[0] == '.' or i[:6] == 'README' or i[-4:] == '.pdf':
            continue
        patient = i
        filein = open(divergence_hiv_dirpath + i, "r")
        filein.readline()
        for j in filein:
            line = j[:-1].split('\t')
            tpoint = line[0]
            syn_divg = line[1]
            nonsyn_divg = line[2]
            sample_stats_dic['%s_%s' % (patient, tpoint)][39] = syn_divg
            sample_stats_dic['%s_%s' % (patient, tpoint)][40] = nonsyn_divg
        filein.close()
    #update header
    header += ['hiv_synonymous_divergence', 'hiv_non_synonymous_divergence']

    #get hiv selection values
    if selection_hiv_dirpath[-1] != '/':
        selection_hiv_dirpath += '/'
    for i in os.listdir(selection_hiv_dirpath):
        if i[0] == '.' or i[:6] == 'README' or i[-4:] == '.pdf':
            continue
        patient = i
        filein = open(selection_hiv_dirpath + i, "r")
        filein.readline()
        for j in filein:
            line = j[:-1].split('\t')
            tpoint = line[0]
            dnds = line[1]
            sample_stats_dic['%s_%s' % (patient, tpoint)][41] = dnds
        filein.close()
    #update header
    header += ['hiv_selection_dN_dS']

    #get abr selection values, but get the weighted mean of the baseline sigma values (not the sigmas that were calculated by convolution)
    if selection_abr_baseline_mean_dirpath[-1] != '/':
        selection_abr_baseline_mean_dirpath += '/'
    for i in os.listdir(selection_abr_baseline_mean_dirpath):
        patient = i
        filein = open(selection_abr_baseline_mean_dirpath+i, "r")
        filein.readline()
        for j in filein:
            line = j[:-1].split('\t')
            tpoint = line[0]
            fwr_sigma = line[1]
            cdr_sigma = line[2]
            if tpoint[-2:] == '.5':
                tpoint = tpoint[:-2]
                sample_stats_dic['%s_%s' % (patient, tpoint)][43] = cdr_sigma
                sample_stats_dic['%s_%s' % (patient, tpoint)][45] = fwr_sigma
            else:
                sample_stats_dic['%s_%s' % (patient, tpoint)][42] = cdr_sigma
                sample_stats_dic['%s_%s' % (patient, tpoint)][44] = fwr_sigma
        filein.close()
    #update header
    header += ['abr_baseline_mean_sigma_CDR', 'abr_baseline_mean_sigma_CDR_duplicate', 'abr_baseline_mean_sigma_FWR', 'abr_baseline_mean_sigma_FWR_duplicate']

    #get hiv sample prep batch info
    #below is a function that will get the desired info from all the sample type files
    def get_dates_and_info(line):
        line = line[:-1].split('\t')
        sample_info = line[0].split('_')
        patient_id = sample_info[1]
        tpoint_rank = sample_info[2]
        date = line[3]
        if date == '':
            return patient_id, tpoint_rank, 'NA'
        if ';' in date:
            date = date.split(';')[-1]
        if len(date.split('/')[-1]) == 2:#if the full year isin't written, add it
            date = '%s20%s' % (date[:-2], date[-2:])
        return patient_id, tpoint_rank, date
    if hiv_sample_prep_dirpath[-1] != '/':
        hiv_sample_prep_dirpath += '/'
    #get RNA extraction date
    fileout = open(hiv_sample_prep_dirpath+'rna_extractions.txt', "r")
    fileout.readline()
    for i in fileout:
        patient_id, tpoint_rank, date = get_dates_and_info(i)
        sample_stats_dic['%s_%s' % (patient_id, tpoint_rank_dic['%s_%s' % (patient_id, tpoint_rank)])][46] = date
    fileout.close()
    #get cDNA synthesis date
    fileout = open(hiv_sample_prep_dirpath+'cDNA.txt', "r")
    fileout.readline()
    for i in fileout:
        patient_id, tpoint_rank, date = get_dates_and_info(i)
        sample_stats_dic['%s_%s' % (patient_id, tpoint_rank_dic['%s_%s' % (patient_id, tpoint_rank)])][47] = date
    fileout.close()
    #get 1st round PCR date
    fileout= open(hiv_sample_prep_dirpath+'PCR_product_1st_round.txt', "r")
    fileout.readline()
    for i in fileout:
        patient_id, tpoint_rank, date = get_dates_and_info(i)
        sample_stats_dic['%s_%s' % (patient_id, tpoint_rank_dic['%s_%s' % (patient_id, tpoint_rank)])][48] = date
    fileout.close()
    #get 2nd round PCR date
    fileout= open(hiv_sample_prep_dirpath+'PCR_product_2nd_round.txt', "r")
    fileout.readline()
    for i in fileout:
        patient_id, tpoint_rank, date = get_dates_and_info(i)
        sample_stats_dic['%s_%s' % (patient_id, tpoint_rank_dic['%s_%s' % (patient_id, tpoint_rank)])][49] = date
    fileout.close()
    #update header
    header += ['hiv_rna_extraction_date', 'hiv_cdna_synthesis_date', 'hiv_1st_round_pcr_date', 'hiv_2nd_round_pcr_date']

    #get HIV seq date and seq library concentration
    filein = open(hiv_seq_run_info_filepath, 'r')
    filein.readline()
    for i in filein:
        line = i[:-1].split('\t')
        sample_info = line[0].split('_')
        patient_id = sample_info[1]
        tpoint_rank = sample_info[2]
        tpoint = tpoint_rank_dic['%s_%s' % (patient_id, tpoint_rank)]
        seq_date = line[1]
        lib_conc = line[2]
        sample_stats_dic['%s_%s' % (patient_id, tpoint)][50] = seq_date
        sample_stats_dic['%s_%s' % (patient_id, tpoint)][51] = lib_conc
    filein.close()
    #update header
    header += ['hiv_sequencing_date', 'hiv_initial_seq_lib_concentration']

    #get ART status
    filein = open(art_status_filepath, "r")
    filein.readline()
    for i in filein:
        line = i[:-1].split('\t')
        patient_id = line[0]
        tpoint = line[2]
        art_status = line[5]
        sample_stats_dic['%s_%s' % (patient_id, tpoint)][52] = art_status
    filein.close()
    #update header
    header += ['ART_niave']

    #Now write sample_stats_dic to output file
    fileout = open(output_filepath, "w")
    fileout.write('\t'.join(header) + '\n')
    for i in sorted(sample_stats_dic):
        #get patient ID and time-point
        patient_id = i.split('_')[0]
        tpoint = i.split('_')[1]
        #add to this list as more stats are added to table.
        #each element must be a string though.
        output_line = [patient_id, tpoint]
        output_line += sample_stats_dic[i]
        fileout.write('\t'.join(output_line) + '\n')
    fileout.close()
    return


if __name__ == '__main__':
    gather_stats_changeo_processed_data(sample_names_key_filepath='/Users/nstrauli/data/abr_hiv_coevo/sample_names_key/key.txt', cd4_count_viral_load_filepath='/Users/nstrauli/data/options_cohort_data/CD4_viral_load/cd4_viral_load_data', hiv_read_depth_dirpath='/Users/nstrauli/data/abr_hiv_coevo/read_depth_foreach_sample/hiv', hiv_pi_diversity_dirpath='/Users/nstrauli/data/abr_hiv_coevo/diversity/hiv/pi', ab_library_conc_filepath='/Users/nstrauli/data/options_cohort_sample_processing/sample_database/PBMCs/seq_libraries.txt', ab_rna_conc_filepath='/Users/nstrauli/data/options_cohort_sample_processing/sample_database/PBMCs/rna_extractions.txt', ab_pi_diversity_dirpath='/Users/nstrauli/data/abr_hiv_coevo/diversity/abr/changeo/pi', samp_seqDate_and_barcodeID_filepath='/Users/nstrauli/data/abr_hiv_coevo/sample_seq_date_and_barcode_id/sample_seq_date_and_barcode_id.txt', irep_barcode_seqs_filepath='/Users/nstrauli/data/abr_hiv_coevo/abr_barcode_analysis/irepertoire_barcod_seqs.txt', selection_abr_baseline_dirpath='/Users/nstrauli/data/abr_hiv_coevo/selection/abr/changeo/baseline/summary_stats_for_samples', divergence_abr_dirpath='/Users/nstrauli/data/abr_hiv_coevo/divergence/abr/changeo/mean_divergence_trajectories/', divergence_hiv_dirpath='/Users/nstrauli/data/abr_hiv_coevo/divergence/hiv/counting_relative_to_most_numerous_1st_tpoint_seq/mean_divergence_trajectories', selection_hiv_dirpath='/Users/nstrauli/data/abr_hiv_coevo/selection/hiv/dN_dS_by_counting_changes_relative_to_most_numerouse_1st_tpoint_seq/mean_selection_trajectories', selection_abr_baseline_mean_dirpath='/Users/nstrauli/data/abr_hiv_coevo/selection/abr/changeo/baseline/mean_selection_trajectories', hiv_sample_prep_dirpath='/Users/nstrauli/data/options_cohort_sample_processing/sample_database/plasma', hiv_seq_run_info_filepath='/Users/nstrauli/data/abr_hiv_coevo/hiv_seq_run_info/info.txt', art_status_filepath='/Users/nstrauli/data/abr_hiv_coevo/sample_art_status/ART_status.txt', output_filepath='/Users/nstrauli/data/abr_hiv_coevo/sample_stats/changeo/sample_stats.txt')
