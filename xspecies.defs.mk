##
## cross-species mapping
##
asm_names = \
    ponAbe-GCA_028885655.2 \
    ponPyg-GCA_028885625.2 \
    gorGor-GCA_029281585.2 \
    panPan-GCA_029289425.2 \
    panTro-GCA_028858775.2 \
    symSyn-GCA_028878055.2

asmname_to_acc = $(word 2,$(subst -, ,${1}))

##
# asmname=
# asmacc=
##i
src_asmname = hg38
src_genome = /hive/data/genomes/hg38/hg38.2bit
src_chromSizes = /hive/data/genomes/hg38/chrom.sizes

xspecies_url = https://public.gi.ucsc.edu/~pnhebbar/t2t_autosomes/assemblyHub
xspecies_rootdir=/hive/users/markd/gencode/projs/cls/cls-t2t-primate/hg38/data/assemblies

xspecies_catbb_url = ${xspecies_url}/${asmacc}/tm.bb
xspecies_genome = ${xspecies_rootdir}/${asmname}/${asmname}.2bit
xspecies_chromsizes_url = ${xspecies_url}/${asmacc}/chrom.sizes
src_xspecies_chains = ${xspecies_rootdir}/${asmname}/${src_asmname}-${asmname}.chain.gz

xspecies_dir = ${root}/bigtest/xspecies/${asmname}
xspecies_chromsizes = ${xspecies_dir}/${asmname}.sizes
xspecies_cat_pre = ${xspecies_dir}/${asmname}.cat
xspecies_cat_bb = ${xspecies_cat_pre}.bb
xspecies_cat_psl = ${xspecies_cat_pre}.psl
xspecies_cat_bed = ${xspecies_cat_pre}.bed
xspecies_cat_fa = ${xspecies_cat_pre}.fa

src_xspecies_rna_aln = ${xspecies_dir}/${src_asmname}-${asmname}.rna.psl
xspecies_mapannot_pre = ${xspecies_dir}/${src_asmname}-${asmname}.mapAnnots
xspecies_mapannot_psl = ${xspecies_mapannot_pre}.psl
xspecies_mapannot_ref_tsv = ${xspecies_mapannot_pre}.ref.tsv
xspecies_mapannot_problems_tsv = ${xspecies_mapannot_pre}.problems.tsv

xspecies_decorators_pre = ${xspecies_dir}/${src_asmname}-${asmname}.decorators
xspecies_decorators_bed = ${xspecies_decorators_pre}.bed
