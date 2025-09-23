
colorlist <- list(
  
  "highlight"=
    c("*"="black",
      "-"="white"),
  
  "type" = c(
    "prim" = "green4",
    "met" = "black",
    "ctrl" = "lightblue3",
    "ctrl tissue" = "lightblue3",
    "primary" = "green4",
    "metastasis" = "black"
  ),
  
  "patient" = c(
    "KSCLC112" = "black",
    "KSCLC127" = "#9D758D",
    "KSCLC147" = "#734631",
    "KSCLC153" = "#7f7f7f",
    "KSCLC155" = "#F1CFA2",
    "KSCLC156" = "#B3956F",
    "KSCLC157" = "#8A91A4",
    "KSCLC158" = "darkred",
    "KSCLC159" = "khaki4",
    "KSCLC163" = "seagreen4",
    "KSCLC172" = "navy",
    "KSCLC176" = "mediumpurple4",
    "KSCLC179" = "purple",
    "KSCLC180" = "violet",
    "KSCLC181" = "khaki1",
    "KSCLC189" = "red",
    "KSCLC205" = "gold",
    "KSCLC206" = "seagreen1",
    "KSCLC207" = "cyan2",
    "KSCLC209" = "mediumpurple",
    "KSCLC210" = "cyan4",
    "KSCLC211" = "pink",
    "KSCLC212" = "orange"
  ),
  
  "origin" = c(
    
    "brain" = "red",
    "brain_stem" = "orange2",
    "cerebellum" = "yellow",
    "liver" = "navy",
    "lymph_node" = "#84d7e1",
    "pericardium" = "#ff95a8",
    "peritoneum" = "blue",
    "lung" = "green4",
    "pleura" = "green2",
    "adrenal_gland" = "purple3",
    "kidney" = "brown",
    "subcutan" = "magenta",
    "spleen" = "lightblue4"
    
  ),
  
  "cohort"=
    c("cohort1"="gold",
      "cohort2"="white"),
  
  "omics_data" = c(
    "nr_proteins" = "#104660",
    " \nNr. of quantified\nproteins per sample, cohort#2\n "= "#8FD7D7",
    "nr_transcripts" = "darkred",
    " \nNr. of quantified\nproteins per sample, cohort#1\n " = "#104660",
    " \nNr. of transcripts\nwith TPM>0 per sample, cohort#1\n " = "darkred",
    "nr_protein_groups" = "#104660"),
  
  "histology" = c(
    "tumor" = "#017a4a",
    "patient" = "black",
    "stroma"= "#FFCE4E",
    "necrosis" = "#a00000",
    "normal" = "lightblue3"
  ),
  
  "quality_issue"=
    c("≤10% tumor"="red",
      "≤10% normal"="blue",
      "10% tumor contamination" = "purple",
      "unknown composition" = "lightgray",
      "no"="white"),
  
  "nr_quant_proteins"=
    circlize::colorRamp2(c(500,1000,3000,5000,7000,9000), 
                         c("white","#EDEDED","#D4D4D4","#8FD7D7",#"#7f7f7f",
                           "#4E95D9","#215F9A")),
  
  "corr" = circlize::colorRamp2(c(0.5,0.6,0.7,0.8,0.9,1), 
                                 c("white","gold","darkorange","firebrick3","darkred","black")),
  
  "ESTIMATE_stromal" = circlize::colorRamp2(c(-2500,0,2500), 
                                            c("darkblue","white", "darkred")),
  "ESTIMATE_immune" = circlize::colorRamp2(c(-2500,0,2500), 
                                           c("darkblue","white", "darkred")),
  
  "z-score" = circlize::colorRamp2(c(-2,0,2), 
                                   c("darkblue","#49b06e","gold")),
  
  "NE_score" =  circlize::colorRamp2(c(-0.2, -0.02, 0.02, 0.2), 
                                     c("#fb8c00","white","white","#7559a2")),
  "nonNE_score" = circlize::colorRamp2(c(-0.2, -0.1, -0.02, 0.02, 0.1,0.2), 
                                       c("#19342a","#2d5253","white","white","#a58462","#b17050")),
  "Park_inflamed_signature" = circlize::colorRamp2(c(-0.3, 0), 
                                                       c("white","pink3")),
  "correlation" = c(
    "not highlighted" = "snow3",
    "no correlation" = "black",
    "positive"= "#a00000",
    "negative" = "#1a80bb"
  ),
  
  "normal" = circlize::colorRamp2(c(0,0.33,0.66,1), 
                                  c("white","#66bbb2","#008d7f","black")),
  
  "cd3_peri" = circlize::colorRamp2(c(0,1,2,3), 
                                    c("white","#fde183","#f09854","#ce5e08")),
  "cd3_intra" = circlize::colorRamp2(c(0,1,2,3), 
                                     c("white","#fde183","#f09854","#ce5e08")),
  
  "frequency_cat" = c(
    "high" = "#215F9A",
    "low" = "white"),
  
  "h_score" = colorRamp2(c(0, 50,75, 100,200,300), 
                         c("black","navyblue","lightblue","white", "orange","red")),
  
  "pd_l1_h_score" = circlize::colorRamp2(c(0,50,300), 
                                         c("white","pink2","pink4")),
  "dll3_h_score" = circlize::colorRamp2(c(0,50,300), 
                                        c("white","lightgreen","darkgreen")),
  "b7_h3_h_score" = circlize::colorRamp2(c(0,50,300), 
                                         c("white","violet","darkviolet")),
  
  "ihc_subtype" = c(
    "A" = "#1E4665",
    "ANP" = "black",
    "AN" = "#794924",
    "AP" = "purple",
    "P" = "#ff363c",
    "N" = "#00985C",
    "ctrl" = "white",
    "QN" = "#3d98d3",
    "QN-Infl." = "lightblue",
    "P-Infl." = "pink1",
    "AP-Infl." = "pink3",
    "A-Infl." = "#7559a2",
    "QN-immune-desert"= "#3d98d3",
    "QN-Immune-desert"= "#3d98d3",
    "QN\nimmune-desert"= "#3d98d3",
    "end" = "white",
    "N/A" = "white"
    
  ),
  
  "cd3_peri_mean" = circlize::colorRamp2(c(0,1,2,3), 
                                         c("white","#fde183","#f09854","#ce5e08")),
  "cd3_intra_mean" = circlize::colorRamp2(c(0,1,2,3), 
                                          c("white","#fde183","#f09854","#ce5e08")),
  
  "pd_l1_h_score_mean" = circlize::colorRamp2(c(0,50,300), 
                                              c("white","pink2","pink4")),
  "dll3_h_score_mean" = circlize::colorRamp2(c(0,50,300), 
                                             c("white","lightgreen","darkgreen")),
  "b7_h3_h_score_mean" = circlize::colorRamp2(c(0,50,300), 
                                              c("white","violet","darkviolet")),
  "log2_lfq" = circlize::colorRamp2(c(15,20,25), 
                                    c("darkblue","#49b06e","gold")),
  
  "rnaseq_consensus_subtype" = c(
    "A" = "#1E4665",
    "I" = "gold",
    "N" = "#00985C",
    "equivocal" = "#7f7f7f",
    "ctrl" = "white",
    
    "A,N" = "#794924",
    "AN" = "#794924",
    "equivocal,N" = "#A5B498",
    "equivocal,I" = "#C1B87E",
    "A,equivocal" = "#AABED6"
    
  ),
  
  "sample_group" = c(

    "primary" ="#1a5354",

    "lymph_node_met" = "#84d7e1",
    "pericardial_met" = "#ff95a8",
    "pericardium_met" = "#ff95a8",
    "liver_met"="#5a9599",
    "pleural_met" ="#ade2d0",
    "pleura_met" ="#ade2d0",
    "peritoneal_met" ="black",
    "peritoneum_met" ="black",
    
    "adrenal_gland_met" ="#ff6f00",
    "brain_met" = "#c71000",
    "brain_stem_met" = "#008ea0",
    "cerebellum_met" = "gold",
    "lung_met" = "#ff6348",
    "end" = "white"
    
  ),
  
  "corr_combined" = circlize::colorRamp2(c(-1,-0.95,-0.9,-0.6,0.6,0.9,0.95,1), 
                                         c("black", "darkred","darkorange", "gold","gold","#49b06e","#104660","black")),
  
  "diagnosis_time_m" = circlize::colorRamp2(c(0,6,12,18), 
                                            c("white","pink1","pink3", "pink4")),
  
  "immune_infiltration_score" = circlize::colorRamp2(c(-5000,-1500,0,1500), 
                                                     c("darkblue","lightblue","white", "#794924")),
  "stromal_infiltration_score" = circlize::colorRamp2(c(-5000,-1500,0,1500), 
                                                      c("darkblue","lightblue","white", "#FFCE4E")),
  "inflamed_score" = circlize::colorRamp2(c(-0.2, -0.02, 0.02, 0.2), 
                                          c("#fb8c00","white","white","#7559a2")),
  "corr2" = circlize::colorRamp2(c(0.6,0.9,0.95,1), 
                                 c("gold","darkorange","darkred","black"))
  
  
  
  
)



