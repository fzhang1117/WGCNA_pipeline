#! /bin/bash

Rscript WGCNA_KEGGenrichment.R "./result/modules/merged_modes.dissTOM.JS.pearson.deepsplit\=2.MEDissThres\=0.10.txt" 35581 "./result/enrichment/KEGG_enrichment/KEGG.JS.pearson.deepsplit\=2.MEDissThres\=0.10.txt"
Rscript WGCNA_KEGGenrichment.R "./result/modules/merged_modes.dissTOM.JC.pearson.deepsplit\=2.MEDissThres\=0.10.txt" 35581 "./result/enrichment/KEGG_enrichment/KEGG.JC.pearson.deepsplit\=2.MEDissThres\=0.10.txt"
Rscript WGCNA_KEGGenrichment.R "./result/modules/merged_modes.dissTOM.HS.pearson.deepsplit\=2.MEDissThres\=0.10.txt" 35581 "./result/enrichment/KEGG_enrichment/KEGG.HS.pearson.deepsplit\=2.MEDissThres\=0.10.txt"
Rscript WGCNA_KEGGenrichment.R "./result/modules/merged_modes.dissTOM.HC.pearson.deepsplit\=2.MEDissThres\=0.10.txt" 35581 "./result/enrichment/KEGG_enrichment/KEGG.HC.pearson.deepsplit\=2.MEDissThre\s=0.10.txt"
