# phospho

Step 1: install devtools:
```
install.packages("devtools")
library(devtools)
```
Step 2: install phospho package:
```
install_github('andyhoss/phospho')
```
Step 3: load the R phospho library and set the path to your raw data file:
```
library(phospho)
path_to_raw_data='path/to/excel/here/'
```
Step 4: download the "final_seq.csv" file from this repo and load this file:
```
final_seqs=data.table::fread('path/to/final_seq.csv')
```
Step 5: prepare the "sites" data:
```
d1 <- data.table::fread(path_to_raw_data) %>% 
  dplyr::filter(EG.ProteinPTMLocations != '') %>%
  dplyr::filter(EG.PTMAssayProbability>=0.75) %>% 
  dplyr::mutate(sample = gsub('_Phos_30SPD', '', gsub('20221205_EXPL2_EvoPRI_ZY_AurA_', '', R.FileName)))

initial_fragments=d1 %>% 
  dplyr::select(PG.ProteinAccessions, PEP.StrippedSequence, 
                FG.LabeledSequence, EG.ProteinPTMLocations, 
                EG.PTMLocalizationProbabilities,EG.PTMAssayProbability, 
                FG.Quantity, sample) %>% 
  dplyr::group_by(PG.ProteinAccessions, FG.LabeledSequence, 
                  EG.ProteinPTMLocations, PEP.StrippedSequence, sample) %>%
  dplyr::summarize(value = sum(FG.Quantity)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(id_cols = c(PG.ProteinAccessions, EG.ProteinPTMLocations,
                                 FG.LabeledSequence, PEP.StrippedSequence),
                     names_from=sample,
                     values_from = value) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sites = getSite(FG.LabeledSequence)) %>% 
  dplyr::select(PG.ProteinAccessions, EG.ProteinPTMLocations, FG.LabeledSequence, PEP.StrippedSequence, sites) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(accessions = stringr::str_split(PG.ProteinAccessions, ';'),
                ptm = stringr::str_split( EG.ProteinPTMLocations, ';')) %>%
  tidyr::unnest(cols = c(accessions, ptm)) %>% 
  dplyr::mutate(ptm =  stringr::str_split( gsub("\\)\\(", "\\);\\(",ptm), ';')) %>%
  tidyr::unnest(cols = ptm) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(site_count = length(unlist(strsplit(gsub('(|)', '', 
                                                         gsub("\\)\\(", ",",ptm)), ',')) %>% 
                                      grep('C|M', ., value=T, invert = T)), 
                residue = unlist(strsplit(gsub('\\(|\\)', '', ptm), ',')) %>% 
                  grep('C|M', ., value=T, invert = T) %>%
                  paste0(., collapse=',') %>%
                  stringr::str_split(., ','),
                sites = stringr::str_split(sites, ';'))  %>% 
  tidyr::unnest(cols = c(residue, sites)) %>% 
  dplyr::filter(residue != "") %>% 
  dplyr::mutate(m=ifelse(site_count>3, 3, site_count), 
                pos = as.numeric(substring(residue, first=2))) %>% 
  dplyr::inner_join(final_seq, by='accessions') %>%
  dplyr::rowwise() %>%
  dplyr::mutate(window = getWindow(pos=pos, peptide=peptide))%>% 
  dplyr::mutate(match = ifelse(nchar(substring(PEP.StrippedSequence, as.numeric(sites)-2, as.numeric(sites)+2))==5 &
                                 substring(PEP.StrippedSequence, as.numeric(sites)-2, as.numeric(sites)+2) == substring(window,6,10), 1, 
                               ifelse(substring(PEP.StrippedSequence, as.numeric(sites)-1, as.numeric(sites)+1) == substring(window,7,9), 1, 
                                      ifelse(sites ==1 & substring(PEP.StrippedSequence, as.numeric(sites), as.numeric(sites)+2) == substring(window,8,10), 1, 
                                             ifelse(pos == nchar(peptide) & substring(PEP.StrippedSequence, as.numeric(sites)-2, as.numeric(sites)) == substring(window,6,8), 1,
                                                    0))))) %>%
  dplyr::ungroup()

#test behavior of (1)(2);(1)(2)(3) syntax
#frags2 %>% dplyr::filter(grepl("\\)\\(", EG.ProteinPTMLocations)) %>% View


###
#fix missing or incorrect 15 AA windows
fixed_fragments <- fixWindows(fragments=initial_fragments)

#dealing with pesky multi-entries
#identify accessions+ptms
##use this to fix FIXED_FRAGMENTS for individual PTM + windows
multi_ids <- fixed_fragments %>% 
  dplyr::group_by(PG.ProteinAccessions, EG.ProteinPTMLocations, site_count) %>% 
  dplyr::summarize(n_ptm=length(unique(ptm))) %>% 
  dplyr::filter(n_ptm>1) %>% 
  dplyr::select(-c(n_ptm, site_count)) %>% 
  dplyr::mutate(ptm = sub("\\)\\(.*","\\)", sub(";.*", '', EG.ProteinPTMLocations))) %>% 
  dplyr::inner_join(fixed_fragments, by=c('PG.ProteinAccessions', 'EG.ProteinPTMLocations', 'ptm'))%>% 
  dplyr::group_by(PG.ProteinAccessions, EG.ProteinPTMLocations, FG.LabeledSequence, residue) %>% 
  dplyr::summarise(window = paste0(unique(window), collapse=',')) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(id=paste0(c(PG.ProteinAccessions, residue), collapse="_")) 

fragments <- fixed_fragments %>% 
  dplyr::filter(!FG.LabeledSequence %in% multi_ids$FG.LabeledSequence) %>%
  dplyr::group_by(PG.ProteinAccessions, EG.ProteinPTMLocations, FG.LabeledSequence, residue) %>% 
  dplyr::summarise(window = paste0(unique(window), collapse=',')) %>%
  dplyr::mutate() %>%
  dplyr::distinct() %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(id=paste0(c(PG.ProteinAccessions, residue), collapse="_")) %>% 
  dplyr::distinct() %>%
  dplyr::bind_rows(multi_ids) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(PG.ProteinAccessions,  EG.ProteinPTMLocations, FG.LabeledSequence) 

ids <- fragments %>%
  dplyr::select(id, FG.LabeledSequence) %>% 
  dplyr::distinct()

sites <- d1 %>% 
  dplyr::select(FG.LabeledSequence, 
                EG.ProteinPTMLocations,
                FG.Quantity, 
                sample) %>% 
  dplyr::filter(EG.ProteinPTMLocations != '') %>%
  dplyr::select(-EG.ProteinPTMLocations) %>%
  dplyr::left_join(ids, by = 'FG.LabeledSequence') %>%
  dplyr::group_by(id, sample) %>%
  dplyr::summarize(value = sum(FG.Quantity)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(id_cols = id,
                     names_from=sample,
                     values_from = value)

```
Step 5: analyze the sites data:
```
ppe <- normalizeData(sites=sites)
QCplots(ppe)
distPlots(ppe)

fitAS_WT <- testPPE(ppe=ppe, contrast='grpsAS_DMSO - grpsWT_DMSO')
fitAS=testPPE(ppe=ppe, contrasts ='grpsAS_1NA - grpsAS_DMSO')
fitMLN=testPPE(ppe=ppe, contrasts ='grpsWT_MLN - grpsWT_DMSO')
fitMK=testPPE(ppe=ppe, contrasts ='grpsWT_MK - grpsWT_DMSO')

#motif
motifAS_WT <- getMotif(fitAS_WT, data=final_sites, 'grpsAS_DMSO - grpsWT_DMSO', LFC=0.6)
motifAS <- getMotif(fitAS, data=final_sites, 'grpsAS_1NA - grpsAS_DMSO', LFC=0)
motifMLN <- getMotif(fitMLN, data=final_sites, 'grpsWT_MLN - grpsWT_DMSO')
motifMK <- getMotif(fitMK, data=final_sites, 'grpsWT_MK - grpsWT_DMSO')
```
