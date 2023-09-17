#!/bin/bash

# dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition
dirPATH=${HOME}/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/RegionEnrichment/Repeats
pyPATH=/mnt/Storage/home/yanghui/scripts/python
shPATH=/mnt/Storage/home/wangyiman/bin/utilities
# shPATH=/mnt/Storage/home/yanghui/scripts/shell
binPATH=/mnt/Storage/home/yanghui/scripts/bin
rPATH=/mnt/Storage/home/yanghui/scripts/R
rPATH_me=/mnt/Storage/home/wangyiman/bin/utilities
anPATH=/mnt/Storage/home/yanghui/annotations

CHMbedPATH=/mnt/Storage/home/wangyiman/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific


enrich_unenrich_LTR_subF_ele(){
	cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/RegionEnrichment/Repeats/LTR_Subfamily
	for regionT in Universal Universal_complementarySet;do
		cat ${regionT}.EnrichScore_LTR_*.txt | grep -v "EnrichScore" | \
			awk 'BEGIN{FS=OFS="\t";}{if($2!="-inf" && $2>=1){print $1;}}' | sort -u > ${regionT}.subfml_enriched.txt # 50 129
		
		cat ${regionT}.EnrichScore_LTR_*.txt | grep -v "EnrichScore" | \
			awk 'BEGIN{FS=OFS="\t";}{if($2=="-inf" || $2<1){print $1;}}' | sort -u | grep -wvf ${regionT}.subfml_enriched.txt - > ${regionT}.subfml_nonenriched.txt # 592 513
	done
}
# enrich_unenrich_LTR_subF_ele


obtain_corresponding_anno(){
	cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/RegionEnrichment/Repeats/LTR_Subfamily
	for regionT in Universal Universal_complementarySet;do
		for group in enriched nonenriched
		do
			for subfml in $(cut -f 1 ${regionT}.subfml_${group}.txt)
			do
				cat ${anPATH}/mm10/Repeats/LTR/*/mm10.${subfml}.bed | \
					intersectBed -a - -b ${anPATH}/mm10/mm10_euch.chrom.limits | \
					cut -f 1-3 >> ${regionT}.subfml_${group}.bed # 41,236 75483 (enriched) || 837,636 803389 (nonenriched)
			done # for subfml end
		done # for group end
	done # for regionT end
}
# obtain_corresponding_anno

ave_conservation_score(){
	cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/RegionEnrichment/Repeats/LTR_Subfamily
	echo -e "xValue\tyValue\tgroup" > phastCons60wayPlacental.ggplot_basic.txt
	for regionT in Universal Universal_complementarySet;do
		# echo -e "xValue\tyValue\tgroup" > ${regionT}_phastCons60wayPlacental.ggplot_basic.txt
		for group in enriched nonenriched
		do
			if [ -f ${regionT}.subfml_${group}.bed ]
			then
				bash ${shPATH}/averageSignalInRegion.sh \
					${regionT}.subfml_${group}.bed \
					${anPATH}/mm10/Conservation/mm10.60way.phastCons60wayPlacental.bw \
					1 \
					${regionT}.subfml_${group}.phastConsPlacental.txt

				# awk -v XVALUE=${group} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "LTR";}' ${regionT}.subfml_${group}.phastConsPlacental.txt >> ${regionT}_phastCons60wayPlacental.ggplot_basic.txt
				awk -v XVALUE=${group} 'BEGIN{FS=OFS="\t"}{print XVALUE, $4, "LTR";}' ${regionT}.subfml_${group}.phastConsPlacental.txt >> phastCons60wayPlacental.ggplot_basic.txt
			fi
		done # for group end
	done # for regionT end
}
ave_conservation_score

boxplot(){
	rPATH=/mnt/Storage/home/yanghui/scripts/R
	cd ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/RegionEnrichment/Repeats/LTR_Subfamily
	# for regionT in Universal Universal_complementarySet;do
	# 	Rscript ${rPATH}/ggplot_boxplot.r \
	# 		${regionT}_phastCons60wayPlacental.ggplot_basic.txt \
	# 		${regionT}_phastCons60wayPlacental.ggplot_boxplot.pdf \
	# 		LTR \
	# 		"" \
	# 		"Conservation score" \
	# 		0 \
	# 		0.4 \
	# 		0.2 \
	# 		"enriched;nonenriched" \
	# 		"LTR" \
	# 		3.5 \
	# 		4.5
	# done

	Rscript ${rPATH}/ggplot_boxplot.r \
		phastCons60wayPlacental.ggplot_basic.txt \
		phastCons60wayPlacental.ggplot_boxplot.pdf \
		LTR \
		"" \
		"Conservation score" \
		0 \
		0.4 \
		0.2 \
		"enriched;nonenriched" \
		"LTR" \
		3.5 \
		4.5
}
boxplot

##### universal CHMs and non-universal CHMs combinely
# Warning messages:
# 1: Removed 49262 rows containing non-finite values (stat_boxplot).
# 2: Removed 49262 rows containing non-finite values (stat_boxplot).
# 3: Removed 1 rows containing missing values (geom_hline).
# null device
#           1
# [1] "wilcox.test of significance between different xValues"
# # A tibble: 1 × 8
#   .y.    group1   group2          p p.adj p.format p.signif method
#   <chr>  <chr>    <chr>       <dbl> <dbl> <chr>    <chr>    <chr>
# 1 yValue enriched nonenriched     0     0 <2e-16   ****     Wilcoxon
# [1] "t.test of significance between different xValues"
# # A tibble: 1 × 8
#   .y.    group1   group2          p p.adj p.format p.signif method
#   <chr>  <chr>    <chr>       <dbl> <dbl> <chr>    <chr>    <chr>
# 1 yValue enriched nonenriched     0     0 <2e-16   ****     T-test


##### universal CHMs and non-universal CHMs separately

# Warning messages:
# 1: Removed 24631 rows containing non-finite values (stat_boxplot).
# 2: Removed 24631 rows containing non-finite values (stat_boxplot).
# 3: Removed 1 rows containing missing values (geom_hline).
# null device
#           1
# [1] "wilcox.test of significance between different xValues"
# # A tibble: 1 × 8
#   .y.    group1   group2          p p.adj p.format p.signif method
#   <chr>  <chr>    <chr>       <dbl> <dbl> <chr>    <chr>    <chr>
# 1 yValue enriched nonenriched     0     0 <2e-16   ****     Wilcoxon
# [1] "t.test of significance between different xValues"
# # A tibble: 1 × 8
#   .y.    group1   group2              p    p.adj p.format p.signif method
#   <chr>  <chr>    <chr>           <dbl>    <dbl> <chr>    <chr>    <chr>
# 1 yValue enriched nonenriched 5.42e-158 5.4e-158 <2e-16   ****     T-test
# Warning messages:
# 1: Removed 24631 rows containing non-finite values (stat_boxplot).
# 2: Removed 24631 rows containing non-finite values (stat_boxplot).
# 3: Removed 1 rows containing missing values (geom_hline).
# null device
#           1
# [1] "wilcox.test of significance between different xValues"
# # A tibble: 1 × 8
#   .y.    group1   group2          p p.adj p.format p.signif method
#   <chr>  <chr>    <chr>       <dbl> <dbl> <chr>    <chr>    <chr>
# 1 yValue enriched nonenriched     0     0 <2e-16   ****     Wilcoxon
# [1] "t.test of significance between different xValues"
# # A tibble: 1 × 8
#   .y.    group1   group2              p     p.adj p.format p.signif method
#   <chr>  <chr>    <chr>           <dbl>     <dbl> <chr>    <chr>    <chr>
# 1 yValue enriched nonenriched 1.50e-288 1.50e-288 <2e-16   ****     T-test

soft_link_box_plot(){
	cd ~/CHMsInOtherContexts/figures
	# for regionT in Universal Universal_complementarySet;do
	# 	ln -sf ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/RegionEnrichment/Repeats/LTR_Subfamily/${regionT}_phastCons60wayPlacental.ggplot_boxplot.pdf make8_box_LTR_conservation_score_${regionT}_phastCons60wayPlacental.pdf
	# done
	ln -sf ~/CHMsInOtherContexts/CellStateTransition/CHMOrganization/Universal_specific/Features/RegionEnrichment/Repeats/LTR_Subfamily/phastCons60wayPlacental.ggplot_boxplot.pdf make8_box_LTR_conservation_score_phastCons60wayPlacental.pdf
}
soft_link_box_plot