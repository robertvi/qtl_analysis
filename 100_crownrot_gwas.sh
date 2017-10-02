#!/bin/bash

#
# fit glmnet model to crownrot data
#

# working in ~/octoploid_mapping/crot_gwas
# using clean-GWA-data_relaxed/stringent as input from github:harrisonlab/popgen/snp/crown_rot_qtl_markers_gwas_1.sh
#

set -eu

export PATH=~/rjv_mnt/cluster/git_repos/qtl_analysis:${PATH}
export PATH=~/programs/plink-1.90beta:${PATH}

NAME=clean-GWA-data_relaxed

#see https://www.cog-genomics.org/plink/1.9/data#recode
#recode to vcf -problem: we loose the phenotype values!
#plink --bfile ${NAME} --recode vcf | vcf-fid|vcf-iid --out ${NAME}_
#plink --bfile ${NAME} --recode vcf-iid --out ${NAME}
#recode to ped
#plink --bfile ${NAME} --recode ped --out ${NAME}
#recode to A (Additive)
#plink --bfile ${NAME} --recode A --out ${NAME}_A
#recode to AD (additive plus dominant)
#plink --bfile ${NAME} --recode AD --out ${NAME}_AD


#recode to A
#alleles recoded to numerical (additive, one line per marker)
#plink --bfile ${NAME} --recode A --out ${NAME}

#impute missing data using knn
impute_missing_data.R ${NAME}.raw ${NAME}_imputed.raw

#fit elasticnet, select alpha and lambda using cross validation
#run_glmnet_cv_alpha_lambda_Araw.R ${NAME}_imputed.raw 10 10 0.1 ${NAME}_glmnet


if false ; then ########################################################

#convert from affycodes to numerical codes
convert_allele_codes.sh all_genotypes_affycodes.csv > all_genotypes_numerical.csv

#single marker kruskal wallis tests
calc_kruskalwallis.R ./uniq_genotypes_numerical.csv ./phenotypes.csv score uniq_kw.csv &

#using either unique or unique split genotypes, use crossvalidation to find best elasticnet lambda and alpha
run_glmnet_cv_alpha_lambda.R ./uniq_genotypes_numerical.csv ./phenotypes.csv score 10 10 0.1 uniq_glmnetcvalpha &


run_stepwise_AIC.R uniq_genotypes_numerical.csv       phenotypes.csv score uniq_kw.csv      0.004 uniq_aic

predict_phenotypes.R uniq_genotypes_numerical.csv       uniq_aic_score_AICstep.csv      uniq_aic_preds.csv

predict_phenotypes.R uniq_genotypes_numerical.csv       uniq_glmnetcvalpha_score_glmnet_cvalpha.csv      uniq_cvalpha_preds.csv


plot_predictions.R phenotypes.csv predictions.png \
    uniq_cvalpha_preds.csv               score cvalpha \
    uniqsplit_cvalpha_preds.csv          score cvalpha_split \
    uniq_ebglmnetlasso_preds.csv         score eblasso \
    uniqsplit_ebglmnetlasso_preds.csv    score eblasso_split \
    uniq_ebglmnetlassoNEG_preds.csv      score eblassoNEG \
    uniqsplit_ebglmnetlassoNEG_preds.csv score eblassoNEG_split \
    uniq_ebglmnetelastic_net_preds.csv      score ebelastic \
    uniqsplit_ebglmnetelastic_net_preds.csv score ebelastic_split \
    uniq_aic_preds.csv                   score aic \
    uniqsplit_aic_preds.csv              score aic_split

plot_predictions.R phenotypes.csv predictions_ebelastic.png \
    uniq_ebglmnetelastic_net_preds.csv      score ebelastic \
    uniqsplit_ebglmnetelastic_net_preds.csv score ebelastic_split

plot_predictions.R phenotypes.csv predictions_eblassoNEG.png \
    uniq_ebglmnetlassoNEG_preds.csv      score eblassoNEG \
    uniqsplit_ebglmnetlassoNEG_preds.csv score eblassoNEG_split

plot_predictions.R phenotypes.csv predictions_eblasso.png \
    uniq_ebglmnetlasso_preds.csv         score eblasso \
    uniqsplit_ebglmnetlasso_preds.csv    score eblasso_split

plot_predictions.R phenotypes.csv predictions_aic.png \
    uniq_aic_preds.csv                   score aic \
    uniqsplit_aic_preds.csv              score aic_split

plot_predictions.R phenotypes.csv predictions_cvalpha.png \
    uniq_cvalpha_preds.csv               score cvalpha \
    uniqsplit_cvalpha_preds.csv          score cvalpha_split

plot_predictions.R phenotypes.csv predictions_cvalpha_aic.png \
    uniq_aic_preds.csv                   score aic \
    uniqsplit_aic_preds.csv              score aic_split \
    uniq_cvalpha_preds.csv               score cvalpha \
    uniqsplit_cvalpha_preds.csv          score cvalpha_split

plot_qtl_vs_mapposn.R \
    vescax4_emxfeonly.csv \
    glmnet_aic_kw.png \
    uniq_kw.csv                                      marker pvalue KW           pvalue \
    uniqsplit_kw.csv                                 marker pvalue KWsplit      pvalue_split \
    uniq_glmnetcvalpha_score_glmnet_cvalpha.csv      marker value  glmnet       effect \
    uniqsplit_glmnetcvalpha_score_glmnet_cvalpha.csv marker value  glmnetsplit  effect_split \
    uniq_aic_score_AICstep.csv                       marker value  AIC          effect \
    uniqsplit_aic_score_AICstep.csv                  marker value  AICsplit     effect_split

fi #####################################################################
