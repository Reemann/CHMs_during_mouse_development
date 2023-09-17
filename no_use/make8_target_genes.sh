#!/bin/bash

    targetGenes(){
        # D: Apr-20-2022 12:27 Wed
        # N: Genes within 2kb of stable CHMs
        mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes
        genesWithin2kb(){
            # T: executed in 1.43s, finished 12:37:11 2022-04-20
            mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb
            manage(){
                for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific DCSpecific
                do
                     windowBed -w 2000 -a ${dirPATH}/CHMOrganization/Universal_specific/${group}.bed -b ${anPATH}/mm10/mm10.refGene.withGeneSymbol.bed | \
                         cut -f 8 | sort -u > ${group}.targetGenes.txt # 477, 158, 85, 117, 333, 6694
                done # for group end
                
                sort -R DCSpecific.targetGenes.txt | head -3000 > DCSpecific.targetGenes_rand1.txt 
                sort -R DCSpecific.targetGenes.txt | head -3000 > DCSpecific.targetGenes_rand2.txt 
                sort -R DCSpecific.targetGenes.txt | head -3000 > DCSpecific.targetGenes_rand3.txt 
            }
            # manage
            
            DAVID_2021(){
                prepare(){
                    # BP_ALL
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650429689442.txt" -O Universal.BP_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430110573.txt" -O EarlyEmbryoSpecific.BP_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430255416.txt" -O PGCSpecific.BP_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430412492.txt" -O SpermSpecific.BP_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430733111.txt" -O RetinalSpecific.BP_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650431134297.txt" -O DCSpecific_rand1.BP_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433251835.txt" -O DCSpecific_rand2.BP_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433433267.txt" -O DCSpecific_rand3.BP_ALL.txt
                    
                    # CC_ALL
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650429894357.txt" -O Universal.CC_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430135432.txt" -O EarlyEmbryoSpecific.CC_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430284696.txt" -O PGCSpecific.CC_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430437186.txt" -O SpermSpecific.CC_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430760375.txt" -O RetinalSpecific.CC_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650431234916.txt" -O DCSpecific_rand1.CC_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433292623.txt" -O DCSpecific_rand2.CC_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433458403.txt" -O DCSpecific_rand3.CC_ALL.txt
                    
                    # MF_ALL
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650429964292.txt" -O Universal.MF_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430160923.txt" -O EarlyEmbryoSpecific.MF_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/tr211E66E8CC5961650430319114.txt" -O PGCSpecific.MF_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430465227.txt" -O SpermSpecific.MF_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430784349.txt" -O RetinalSpecific.MF_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650431268625.txt" -O DCSpecific_rand1.MF_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433318768.txt" -O DCSpecific_rand2.MF_ALL.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433478106.txt" -O DCSpecific_rand3.MF_ALL.txt

                    # KEGG
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650429988439.txt" -O Universal.KEGG.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430183836.txt" -O EarlyEmbryoSpecific.KEGG.txt
                    wget -c "https://david.ncifcrf.gov/data/download/tr211E66E8CC5961650430498930.txt" -O SpermSpecific.KEGG.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430806571.txt" -O RetinalSpecific.KEGG.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433010184.txt" -O DCSpecific_rand1.KEGG.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433341879.txt" -O DCSpecific_rand2.KEGG.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433497196.txt" -O DCSpecific_rand3.KEGG.txt

                    # BIOGRID_INTERACTION
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430019286.txt" -O Universal.BIOGRID_INTERACTION.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430213512.txt" -O EarlyEmbryoSpecific.BIOGRID_INTERACTION.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_11E66E8CC5961650430833732.txt" -O RetinalSpecific.BIOGRID_INTERACTION.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433038448.txt" -O DCSpecific_rand1.BIOGRID_INTERACTION.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433374195.txt" -O DCSpecific_rand2.BIOGRID_INTERACTION.txt
                    wget -c "https://david.ncifcrf.gov/data/download/t2t_69CB942DB5841650433519683.txt" -O DCSpecific_rand3.BIOGRID_INTERACTION.txt
                }
                # prepare
                
                plot(){
                    # D: Apr-20-2022 13:48 Wed
                    # for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific DCSpecific DCSpecific_rand1 DCSpecific_rand2 DCSpecific_rand3
                    for group in DCSpecific_rand1 DCSpecific_rand2 DCSpecific_rand3
                    do
                        for type in BP_ALL CC_ALL MF_ALL KEGG BIOGRID_INTERACTION
                        do
                            if [ -f ${group}.${type}.txt ]
                            then
                                python ${binPATH}/GO.py ${group}.${type}.txt
                                Rscript ${rPATH}/GO.r ${group}.${type}.plot.txt ${group}.${type}.plot.pdf ${group} 6 3
                            fi
                        done # for type end
                    done # for group end
                }
                # plot
                
                special_KEGG(){
                    wget -c "https://david.ncifcrf.gov/data/download/tr2C2BD5B1F15B71650452360159.txt" -O RetinalSpecific.KEGG.txt
                }
                # special_KEGG
                
                enriched_GO_Term_genome_enrichment(){
                    # D: May-16-2022 22:55 Mon
                    GO_0007155_cell_adhesion(){
                        # cut -f 2 GO_0007155_cell_adhesion.txt | sort -u > GO_0007155_cell_adhesion.genes.txt # 1,488
                        # python ${pyPATH}/ExtractSubset.py ExtractSubset GO_0007155_cell_adhesion.genes.txt ${anPATH}/mm10/mm10.refGene.withGeneSymbol.bed 4 GO_0007155_cell_adhesion.genes.bed # 3,028
                        # cut -f 1-3 GO_0007155_cell_adhesion.genes.bed | sort -u | sort -k1,1 -k2,2n | mergeBed -i - > GO_0007155_cell_adhesion.genes.merged.bed # 1,453
                          
                        bash ${shPATH}/regionEnrichment_userDefined.sh \
                           "${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb" \
                           "mm10_euch" \
                           "${anPATH}/mm10/mm10.chrom.sizes" \
                           "${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb/GO_0007155_cell_adhesion.genes.merged.bed" \
                           "${anPATH}/mm10/RegionEnrichment/mm10_euch.promoter.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.geneBody.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intergenic.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.exon.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intron.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI_Unmasked.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.SINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LTR.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.PAD.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.potentialEnhancer.merged.bed;" \
                           "GO_0007155_cell_adhesion.genes"
                    }
                    GO_0007155_cell_adhesion
                    
                    GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules(){
                        # D: May-16-2022 23:39 Mon
                        cut -f 2 GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules.txt | sort -u > GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules.genes.txt # 216
                        python ${pyPATH}/ExtractSubset.py ExtractSubset GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules.genes.txt ${anPATH}/mm10/mm10.refGene.withGeneSymbol.bed 4 GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules.genes.bed # 467
                        cut -f 1-3 GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules.genes.bed | sort -u | sort -k1,1 -k2,2n | mergeBed -i - > GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules.genes.merged.bed # 209
                          
                        bash ${shPATH}/regionEnrichment_userDefined.sh \
                           "${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb" \
                           "mm10_euch" \
                           "${anPATH}/mm10/mm10.chrom.sizes" \
                           "${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb/GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules.genes.merged.bed" \
                           "${anPATH}/mm10/RegionEnrichment/mm10_euch.promoter.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.geneBody.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intergenic.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.exon.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intron.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI_Unmasked.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.SINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LTR.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.PAD.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.potentialEnhancer.merged.bed;" \
                           "GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules.genes"
                    }
                    # GO_0098742_cell-cell_adhesion_via_plasma-membrane_adhesion_molecules
                    
                    GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules(){
                        # D: May-16-2022 23:46 Mon
                        cut -f 2 GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules.txt | sort -u > GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules.genes.txt # 117
                        python ${pyPATH}/ExtractSubset.py ExtractSubset GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules.genes.txt ${anPATH}/mm10/mm10.refGene.withGeneSymbol.bed 4 GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules.genes.bed # 212
                        cut -f 1-3 GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules.genes.bed | sort -u | sort -k1,1 -k2,2n | mergeBed -i - > GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules.genes.merged.bed # 111
                          
                        bash ${shPATH}/regionEnrichment_userDefined.sh \
                           "${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb" \
                           "mm10_euch" \
                           "${anPATH}/mm10/mm10.chrom.sizes" \
                           "${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb/GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules.genes.merged.bed" \
                           "${anPATH}/mm10/RegionEnrichment/mm10_euch.promoter.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.geneBody.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intergenic.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.exon.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intron.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI_Unmasked.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.SINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LTR.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.PAD.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.potentialEnhancer.merged.bed;" \
                           "GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules.genes"
                    }
                    # GO_0007156_homophilic_cell_adhesion_via_plasma_membrane_adhesion_molecules
                    
                    all_enriched_GO_terms(){
                        # D: May-17-2022 09:56 Tue
                        cat *.BP_ALL.plot.txt *.MF_ALL.plot.txt *.CC_ALL.plot.txt | \
                            awk -F '[~\t]' 'BEGIN{OFS="\t"}{if($3>=1.30103){print $1;}}' | \
                            sort -u > all_enriched_GO_terms.txt # 711 => too much => only showed
                    }
                    # all_enriched_GO_terms
                    
                    showed_enriched_GO_terms(){
                        # # D: May-17-2022 10:27 Tue
                        # for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific
                        # do
                        #     for category in BP MF CC
                        #     do
                        #         sort -u ${group}.${category}_ALL.plot.txt | \
                        #             awk -F '[~\t]' 'BEGIN{OFS="\t"}{if($3>=1.30103){print $1;}}' | \
                        #             sort -t$'\t' -k3,3nr | head -10 | cut -f 1 >> showed_enriched_GO_terms.txt
                        #     done # for category end
                        # done # for group end
                        # 
                        # # DCSpecific
                        # for category in BP MF CC
                        # do
                        #     sort -u DCSpecific_rand1.${category}_ALL.plot.txt | \
                        #         awk -F '[~\t]' 'BEGIN{OFS="\t"}{if($3>=1.30103){print $1;}}' | \
                        #         sort -t$'\t' -k3,3nr | head -10 | cut -f 1 >> showed_enriched_GO_terms.txt
                        # done # for category end
                        # 
                        # sort -u showed_enriched_GO_terms.txt > showed_enriched_GO_terms.txt.tmp && mv showed_enriched_GO_terms.txt.tmp showed_enriched_GO_terms.txt # 103
                        
                        # # GO term enrichment
                        # echo -e "xValue\tyValue\tgroup\tcol\tlabel" > showed_enriched_GO_terms_used.ggplot_sina.txt
                        # for GOTerm in $(cat showed_enriched_GO_terms_used.txt)
                        # do
                        #     cut -f 2 ${GOTerm}.txt | sort -u > ${GOTerm}.genes.txt
                        #     python ${pyPATH}/ExtractSubset.py ExtractSubset ${GOTerm}.genes.txt ${anPATH}/mm10/mm10.refGene.withGeneSymbol.bed 4 ${GOTerm}.genes.bed
                        #     cut -f 1-3 ${GOTerm}.genes.bed | sort -u | sort -k1,1 -k2,2n | mergeBed -i - > ${GOTerm}.genes.merged.bed
                        #     bash ${shPATH}/regionEnrichment_userDefined.sh \
                        #        "${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb" \
                        #        "mm10_euch" \
                        #        "${anPATH}/mm10/mm10.chrom.sizes" \
                        #        "${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb/${GOTerm}.genes.merged.bed" \
                        #        "${anPATH}/mm10/RegionEnrichment/mm10_euch.promoter.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.geneBody.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intergenic.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.exon.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.intron.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.CGI_Unmasked.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.SINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LINE.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.LTR.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.PAD.merged.bed;${anPATH}/mm10/RegionEnrichment/mm10_euch.potentialEnhancer.merged.bed;" \
                        #        "${GOTerm}.genes"
                        # 
                        #     # Prepare for Sina plot
                        #     
                        #     tail -n +2 ${GOTerm}.genes.EnrichScore.txt | \
                        #         awk 'BEGIN{FS=OFS="\t"}{if($2!="-inf"){print $1, $2, "All", "All", "All";}}' >> showed_enriched_GO_terms_used.ggplot_sina.txt
                        # done # for GOTerm end
                        
                        Sina_plot(){
                            Rscript ${rPATH}/ggplot_sina_colorfull.r \
                                "showed_enriched_GO_terms_used.ggplot_sina.txt" \
                                "showed_enriched_GO_terms_used.ggplot_sina.pdf" \
                                "Enrichment score" \
                                "-8;8;2" \
                                "promoter;geneBody;intergenic;exon;intron;CGI_Unmasked;SINE;LINE;LTR" \
                                "All" \
                                "All" \
                                "Black" \
                                "0.2" \
                                "0.75" \
                                "1;3" \
                                "5" \
                                "5"
                        }
                        # Sina_plot
                        
                        violin_plot(){
                            # D: May-17-2022 12:04 Tue
                            echo -e "xValue\tyValue\tgroup" > showed_enriched_GO_terms_used.ggplot_violin.txt
                            for GOTerm in $(cat showed_enriched_GO_terms_used.txt)
                            do
                                tail -n +2 ${GOTerm}.genes.EnrichScore.txt | \
                                    awk 'BEGIN{FS=OFS="\t"}{if($2!="-inf"){print $1, $2, "All";}}' >> showed_enriched_GO_terms_used.ggplot_violin.txt
                            done # for GOTerm
                            
                            Rscript ${rPATH}/ggplot_violin.r \
                                "showed_enriched_GO_terms_used.ggplot_violin.txt" \
                                "showed_enriched_GO_terms_used.ggplot_violin.pdf" \
                                "" \
                                "" \
                                "Enrichment score" \
                                "-3" \
                                "3" \
                                "1" \
                                "promoter;geneBody;intergenic;exon;intron;CGI_Unmasked;SINE;LINE;LTR" \
                                "All" \
                                "5" \
                                "5"
                        }
                        # violin_plot
                        
                        jitter_plot(){
                            # D: May-17-2022 12:11 Tue
                            Rscript ${rPATH}/ggplot_jitter.r \
                                "showed_enriched_GO_terms_used.ggplot_violin.txt" \
                                "showed_enriched_GO_terms_used.ggplot_jitter.pdf" \
                                "" \
                                "" \
                                "Enrichment score" \
                                "-2" \
                                "3" \
                                "1" \
                                "promoter;geneBody;intergenic;exon;intron;CGI_Unmasked;SINE;LINE;LTR" \
                                "All" \
                                "6" \
                                "4" \
                                "1"
                        }
                        jitter_plot
                        
                        boxplot_plot(){
                            # D: May-17-2022 13:26 Tue                        
                            Rscript ${rPATH}/ggplot_boxplot.r \
                                "showed_enriched_GO_terms_used.ggplot_violin.txt" \
                                "showed_enriched_GO_terms_used.ggplot_boxplot.pdf" \
                                "" \
                                "" \
                                "Enrichment score" \
                                "-2" \
                                "3" \
                                "1" \
                                "promoter;geneBody;intergenic;exon;intron;CGI_Unmasked;SINE;LINE;LTR" \
                                "All" \
                                "6" \
                                "4" \
                                "1"
                        }
                        # boxplot_plot
                    }
                    # showed_enriched_GO_terms
                }
                # enriched_GO_Term_genome_enrichment
            }
            # DAVID_2021
            
            expr(){
                # D: May-31-2022 14:41 Tue
                mkdir -p ${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb/Expression;cd ${dirPATH}/CHMOrganization/Universal_specific/Features/TargetGenes/Within_2kb/Expression
                prepare_expr(){
                    echo -e "xValue\tyValue\tgroup" > targetGenes.expr.ggplot_basic.txt
                    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific DCSpecific
                    do
                        windowBed -w 2000 -a ${dirPATH}/CHMOrganization/Universal_specific/${group}.bed -b ${anPATH}/mm10/mm10.refGene.withGeneSymbol.bed | \
                            cut -f 7 | sort -u > ${group}.targetGenes_RefSeq.txt # 
                            
                        # Early embryos
                        for stage in 2-cell 8-cell Morula ICM
                        do
                            for rep in 1 2 3 4 5 6
                            do
                                if [ -d ${dataPATH_expr_earlyEmbryo}/${stage}_${rep} ]
                                then
                                    python ${pyPATH}/ExtractSubset.py ExtractSubset ${group}.targetGenes_RefSeq.txt ${dataPATH_expr_earlyEmbryo}/${stage}_${rep}/t_data.ctab 5 ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt
                                    awk -v XVALUE=${stage} -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print XVALUE, log($12+1)/log(2), GROUP;}' ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt >> targetGenes.expr.ggplot_basic.txt
                                fi
                            done # for rep end
                        done # for stage end
                        
                        # PGC development
                        for stage in E10.5 E13.5_female E13.5_male
                        do
                            for rep in 1 2 3 4
                            do
                                if [ -d ${dirPATH}/PGCsDevelopment/${stage}.RNASeq.rep${rep} ]
                                then
                                    python ${pyPATH}/ExtractSubset.py ExtractSubset ${group}.targetGenes_RefSeq.txt ${dirPATH}/PGCsDevelopment/${stage}.RNASeq.rep${rep}/${stage}.RNASeq.rep${rep}.t_data.ctab 5 ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt
                                    awk -v XVALUE=${stage} -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print XVALUE, log($12+1)/log(2), GROUP;}' ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt >> targetGenes.expr.ggplot_basic.txt
                                fi
                            done # for rep end
                        done # for stage end
                        
                        # Spermatogenesis
                        for stage in GS PS RS
                        do
                            for rep in 1 2
                            do
                                if [ -d ${dirPATH}/Spermatogenesis_GSE137744/${stage}.RNASeq.rep${rep} ]
                                then
                                    python ${pyPATH}/ExtractSubset.py ExtractSubset ${group}.targetGenes_RefSeq.txt ${dirPATH}/Spermatogenesis_GSE137744/${stage}.RNASeq.rep${rep}/${stage}.RNASeq.rep${rep}.t_data.ctab 5 ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt
                                    awk -v XVALUE=${stage} -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print XVALUE, log($12+1)/log(2), GROUP;}' ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt >> targetGenes.expr.ggplot_basic.txt
                                fi
                            done # for rep end
                        done # for stage end
                        
                        # Retinal development
                        for stage in E14.5 E17.5 P0 P3 P7 P10 P14 P21
                        do
                            for rep in 1 2
                            do
                                if [ -d ${dirPATH}/RetinalDevelopment/GSE87064_Mouse/${stage}.RNASeq.rep${rep} ]
                                then
                                    python ${pyPATH}/ExtractSubset.py ExtractSubset ${group}.targetGenes_RefSeq.txt ${dirPATH}/RetinalDevelopment/GSE87064_Mouse/${stage}.RNASeq.rep${rep}/${stage}.RNASeq.rep${rep}.t_data.ctab 5 ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt
                                    awk -v XVALUE=${stage} -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print XVALUE, log($12+1)/log(2), GROUP;}' ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt >> targetGenes.expr.ggplot_basic.txt
                                fi
                            done # for rep end
                        done # for stage end
                        
                        # DC differentiation
                        for stage in MPP CDP cDC pDC
                        do
                            for rep in 1
                            do
                                if [ -d ${dirPATH}/DendriticCellDevelopment/${stage}.RNASeq.rep${rep} ]
                                then
                                    python ${pyPATH}/ExtractSubset.py ExtractSubset ${group}.targetGenes_RefSeq.txt ${dirPATH}/DendriticCellDevelopment/${stage}.RNASeq.rep${rep}/${stage}.RNASeq.rep${rep}.t_data.ctab 5 ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt
                                    awk -v XVALUE=${stage} -v GROUP=${group} 'BEGIN{FS=OFS="\t"}{print XVALUE, log($12+1)/log(2), GROUP;}' ${group}.targetGenes_RefSeq.expr_${stage}.${rep}.txt >> targetGenes.expr.ggplot_basic.txt
                                fi
                            done # for rep end
                        done # for stage end
                    done # for group end
                }
                # prepare_expr
                
                ggplot_boxplot(){
                    # D: May-31-2022 15:41 Tue
                    for group in Universal EarlyEmbryoSpecific PGCSpecific SpermSpecific RetinalSpecific DCSpecific
                    do
                        Rscript ${rPATH}/ggplot_boxplot.r \
                            targetGenes.expr.ggplot_basic.txt \
                            targetGenes_${group}.expr.ggplot_boxplot.pdf \
                            ${group} \
                            "" \
                            "log2(FPKM+1)" \
                            "0" \
                            "6" \
                            "2" \
                            "2-cell;8-cell;Morula;ICM;E10.5;E13.5_female;E13.5_male;GS;PS;RS;E14.5;E17.5;P0;P3;P7;P10;P14;P21;MPP;CDP;cDC;pDC" \
                            ${group} \
                            "6" \
                            "4" \
                            "2"
                    done # for group end     
                }
                ggplot_boxplot
            }
            expr
        }
        genesWithin2kb
    }
    targetGenes
    