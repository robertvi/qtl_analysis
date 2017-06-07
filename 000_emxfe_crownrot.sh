#!/bin/bash

#
# fit various models to EMxFE crownrot data
#

# working in ~/octoploid_mapping/crot_analysis (formerly ~/octoploid_mapping/qmsim)

set -eu

export PATH=~/rjv_mnt/cluster/git_repos/qtl_analysis:${PATH}

if false ; then ########################################################

#get unique and complete, population ordered genotypes with correct column headers (drop the phase column)
cat ../consensus_map4/popn_EMxFE/map/uniq_affycodes.csv | cut -d, -f3,4,5 --complement > uniq_genotypes_affycodes.csv
cat ../consensus_map4/popn_EMxFE/map/all_affycodes.csv | cut -d, -f3,4,5 --complement > all_genotypes_affycodes.csv

#convert from affycodes to numerical codes
convert_allele_codes.sh uniq_genotypes_affycodes.csv > uniq_genotypes_numerical.csv
convert_allele_codes.sh all_genotypes_affycodes.csv > all_genotypes_numerical.csv

#convert to numerical codes, but split hkxhk markers into separate maternal and paternal pseudo markers
convert_allele_codes_splitshared.sh uniq_genotypes_affycodes.csv > uniq_genotypes_numerical_split.csv
convert_allele_codes_splitshared.sh all_genotypes_affycodes.csv > all_genotypes_numerical_split.csv

#single marker kruskal wallis tests
calc_kruskalwallis.R ./uniq_genotypes_numerical.csv ./phenotypes.csv score uniq_kw.csv &
calc_kruskalwallis.R ./uniq_genotypes_numerical_split.csv ./phenotypes.csv score uniqsplit_kw.csv &

#using either unique or unique split genotypes, use crossvalidation to find best elasticnet lambda and alpha
run_glmnet_cv_alpha_lambda.R ./uniq_genotypes_numerical.csv ./phenotypes.csv score 10 10 0.1 uniq_glmnetcvalpha &
run_glmnet_cv_alpha_lambda.R ./uniq_genotypes_numerical_split.csv ./phenotypes.csv score 10 10 0.1 uniqsplit_glmnetcvalpha &


for model in elasticnet lasso lassoNEG
do
    nice run_ebglmnet.R ./uniq_genotypes_numerical.csv       ./phenotypes.csv score ${model} 10 uniq_ebglmnet &
    nice run_ebglmnet.R ./uniq_genotypes_numerical_split.csv ./phenotypes.csv score ${model} 10 uniqsplit_ebglmnet &
done

plot_qtl_vs_mapposn.R \
    vescax4_emxfeonly.csv \
    glmnet_kw.png \
    uniq_kw.csv                                      marker pvalue KW       pvalue \
    uniqsplit_kw.csv                                 marker pvalue KWsplit  pvalue_split \
    uniq_glmnetcvalpha_score_glmnet_cvalpha.csv      marker value glmnet    effect \
    uniqsplit_glmnetcvalpha_score_glmnet_cvalpha.csv marker value glmnetsplit    effect_split

predict_phenotypes.R uniq_genotypes_numerical.csv       uniq_glmnetcvalpha_score_glmnet_cvalpha.csv      uniq_cvalpha_preds.csv
predict_phenotypes.R uniq_genotypes_numerical_split.csv uniqsplit_glmnetcvalpha_score_glmnet_cvalpha.csv uniqsplit_cvalpha_preds.csv

for model in lasso lassoNEG #elasticnet
do
    predict_phenotypes.R uniq_genotypes_numerical.csv       uniq_ebglmnet_score_ebglmnet${model}_model.csv      uniq_ebglmnet${model}_preds.csv
    predict_phenotypes.R uniq_genotypes_numerical_split.csv uniqsplit_ebglmnet_score_ebglmnet${model}_model.csv uniqsplit_ebglmnet${model}_preds.csv
done

fi #####################################################################

plot_predictions.R phenotypes.csv predictions.png \
    uniq_cvalpha_preds.csv               score cvalpha \
    uniqsplit_cvalpha_preds.csv          score cvalpha_split \
    uniq_ebglmnetlasso_preds.csv         score eblasso \
    uniqsplit_ebglmnetlasso_preds.csv    score eblasso_split \
    uniq_ebglmnetlassoNEG_preds.csv      score eblassoNEG \
    uniqsplit_ebglmnetlassoNEG_preds.csv score eblassoNEG_split
