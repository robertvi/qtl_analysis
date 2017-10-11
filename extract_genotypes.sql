#
# dump genotypes of all samples as a list to csv
#

SELECT
    g.sample_id, a.snp_id, g.genotype, m.ref, m.alt, g.pipeline_id, a.platform, a.probe_id, a.marker_id, 
FROM
    genotype g JOIN
    alias a ON g.alias_id = a.id JOIN
    marker m ON a.marker_id = m.id;
