#!/bin/bash

#
# fit glmnet model to crownrot data
#

# working in ~/octoploid_mapping/crot_gwas
# using clean-GWA-data_relaxed/stringent as input from github:harrisonlab/popgen/snp/crown_rot_qtl_markers_gwas_1.sh
#

set -eu

export PATH=~/rjv_mnt/cluster/git_repos/qtl_analysis:${PATH}
export PATH=~/rjv_mnt/cluster/git_repos/axiom_strawberry/mysql_sample_database:${PATH}
export PATH=~/git_repos/axiom_strawberry/mysql_sample_database:${PATH}
export PATH=~/git_repos/axiom_strawberry/qtl_analysis:${PATH}
export PATH=~/rjv_mnt/cluster/programs/plink-1.90beta:${PATH}
export SCRIPT_PATH=~/git_repos/qtl_analysis

NAME=clean-GWA-data_relaxed

#recode to A (additive numerical encoding)
#alleles recoded to numerical (additive, one line per marker)
plink --bfile ${NAME} --recode A --out ${NAME}_A

#impute missing data using knn
impute_missing_data.R ${NAME}_A.raw ${NAME}_Aimputed.raw

#add dominance info to imputed values (plink cannot reload the .raw file to add them itself)
add_dominance_info.R ${NAME}_Aimputed.raw ${NAME}_ADimputed.raw

#fit elasticnet, select alpha and lambda using cross validation
run_glmnet_cv_alpha_lambda_Araw.R   ${NAME}_Aimputed.raw  10  10  0.1  ${NAME}_Aimputed
run_glmnet_cv_alpha_lambda_Araw.R   ${NAME}_ADimputed.raw  10  10  0.1  ${NAME}_ADimputed

#predict values for the training set only
predict_phenotypes_Araw.R   ${NAME}_ADimputed.raw  ${NAME}_ADimputed_glmnet_cvalpha.csv  ${NAME}_ADimputed_glmnet_pred.csv

#insert list of all samples we want to include into tmp_sample_list
insert_sample_ids.py 

#dump all genotype data for samples requiring predictions
mysql -B -h mongo -u vicker -p$(cat /home/vicker/passwords/mysql_mongo_vicker) -D strawberry_samples  <<XXX | gzip > selections.tsv.gz
select
    m.name,a.marker_id,g.sample_id,g.genotype,g.pipeline_id,m.ref,m.alt
from
    genotype g
        join alias a on g.alias_id=a.id
        join marker m on a.marker_id=m.id
where
    g.sample_id in (select l.sample_id from tmp_sample_list l)
    and
    m.id not in (select t.marker_id from marker_tag t where t.tag='multiform');
XXX

#selections only
#g.sample_id in (2194,2013,2061,2163,1943,1991,2087,2218,2185,1949,2047,2049,1907,2230,2053,2206,2208,1931,2027,2123,2171)

#convert dump into LGEN file
dump2lgen3.py   selections.tsv.gz   selections_lgen

#convert into plink's recodeA raw format
plink --lfile selections_lgen   --recode A   --out selections_raw

#impute missing data using knn
impute_missing_data.R selections_raw.raw selections_imputed.raw

#add dominance info to imputed values (plink cannot reload the .raw file to add them itself)
add_dominance_info.R selections_imputed.raw selections_ADimputed.raw

#predict values for the selections
predict_phenotypes_Araw.R   selections_ADimputed.raw  ${NAME}_ADimputed_glmnet_cvalpha.csv  selections_ADimputed_glmnet_pred.csv
