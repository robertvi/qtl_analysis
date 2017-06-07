#!/bin/bash

#
# fit various models to EMxFE crownrot data
#

# working in ~/octoploid_mapping/crot_analysis (formerly ~/octoploid_mapping/qmsim)

export PATH=~/rjv_mnt/cluster/git_repos/qtl_analysis:${PATH}

#get unique and complete, population ordered genotypes with correct column headers (drop the phase column)
cat ../consensus_map4/popn_EMxFE/map/uniq_affycodes.csv | cut -d, -f3,4,5 --complement > uniq_genotypes_affycodes.csv
cat ../consensus_map4/popn_EMxFE/map/all_affycodes.csv | cut -d, -f3,4,5 --complement > all_genotypes_affycodes.csv

#convert from affycodes to numerical codes
convert_allele_codes.sh uniq_genotypes_affycodes.csv > uniq_genotypes_numerical.csv
convert_allele_codes.sh all_genotypes_affycodes.csv > all_genotypes_numerical.csv

#convert to numerical codes, but split hkxhk markers into separate maternal and paternal pseudo markers
convert_allele_codes_splitshared.sh uniq_genotypes_affycodes.csv > uniq_genotypes_numerical_split.csv
convert_allele_codes_splitshared.sh all_genotypes_affycodes.csv > all_genotypes_numerical_split.csv

if false ; then ########################################################

#for alpha in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1
#do
#    run_glmnet.R ./genotypes.csv ./phenotypes.csv score ${alpha} 10 10 testglmnet &
#done

#run_glmnet_cv_alpha_lambda.R ./genotypes.csv ./phenotypes.csv score 10 100 0.05 testcvalpha &

#for model in elasticnet lasso lassoNEG
#do
#    run_ebglmnet.R ./genotypes.csv ./phenotypes.csv score ${model} 10 testebglmnet &
#done

#calc_kruskalwallis.R genotypes.csv phenotypes.csv score testkw.csv &


plot_qtl_vs_mapposn.R \
    vescax4_emxfeonly.csv \
    testplot.png \
    testkw.csv                                  marker pvalue KW         pvalue \
    testglmnet_score_glmnet1.csv                marker value  glmnet1    effect \
    testglmnet_score_glmnet0.5.csv              marker value  glmnet0.5  effect \
    testebglmnet_score_ebglmnetelastic_net.csv  marker pvalue ebelastic  pvalue \
    testebglmnet_score_ebglmnetelastic_net.csv  marker beta   ebelastic  effect \
    testebglmnet_score_ebglmnetlasso.csv        marker pvalue eblasso    pvalue \
    testebglmnet_score_ebglmnetlasso.csv        marker beta   eblasso    effect

plot_qtl_vs_mapposn.R \
    vescax4_emxfeonly.csv \
    testplot.png \
    testkw.csv                                  marker pvalue KW         pvalue \
    $(for alpha in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 ; do echo testglmnet_score_glmnet${alpha}.csv marker value glmnet${alpha} effect ; done)

plot_qtl_vs_mapposn.R \
    vescax4_emxfeonly.csv \
    glmnet_ebglmnet_kw.png \
    testkw.csv                                  marker pvalue KW       pvalue \
    cv_alpha_score_glmnet_cvalpha.csv           marker value glmnet    effect \
    testebglmnet_score_ebglmnetelastic_net.csv  marker pvalue ebelastic  pvalue \
    testebglmnet_score_ebglmnetelastic_net.csv  marker beta   ebelastic  effect \
    testebglmnet_score_ebglmnetlasso.csv        marker pvalue eblasso    pvalue \
    testebglmnet_score_ebglmnetlasso.csv        marker beta   eblasso    effect


plot_qtl_vs_mapposn.R \
    vescax4_emxfeonly.csv \
    glmnet_kw.png \
    testkw.csv                                  marker pvalue KW       pvalue \
    cv_alpha_score_glmnet_cvalpha.csv           marker value glmnet    effect
fi #####################################################################
