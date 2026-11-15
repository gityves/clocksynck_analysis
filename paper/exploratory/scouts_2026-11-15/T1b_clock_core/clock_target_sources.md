# Direct core-clock target gene sources (scout, 2026-11-15)
No human muscle expression data opened; nothing downloaded. GEO pages were behind a captcha, so GEO details come from OmicsDI mirrors and PMC text; items marked (unverified) are from memory/not confirmed and must be checked on GEO.

| # | Source | Species | Tissue | Definition | n genes | File / URL | Public? size | Muscle-matched? |
|---|---|---|---|---|---|---|---|---|
| A1 | Dyar 2014 Mol Metab (mKO, TA and SOL, Ctrl vs mKO, 3 reps x 6 CT 0-20, every 4 h) | mouse | TA, soleus | KO vs control time course | to derive | GEO GSE43071 (PRJNA184340); processed table size not checked | yes | YES (best: several time points, KO) |
| A2 | Hodge 2015 Skeletal Muscle (iMS-Bmal1-/-, gastrocnemius, CT18-38 every 4 h, pooled RNA, microarray) | mouse | gastrocnemius | KO vs WT reference GSE54652 (1,628 circadian mRNAs, JTK adj P<0.05) | 1,628 rhythmic WT; only 21 up / 12 down listed (Tables 2-3) | PMC4440511; KO accession not given in text (search GEO); supp files are PDFs only | WT yes; KO accession unknown | yes |
| A3 | Kumar 2024 (muscle Bmal1 rescue in Bmal1-KO) | mouse | muscle | rescue | to derive | GSE232957 (design unverified; paper link unverified) | yes | yes, whole-body KO background |
| A4 | Esser lab Bmal1 KO +/- exercise | mouse | gastrocnemius | KO vs ctrl, exercise | - | GSE255280, GSE240653 (not obviously a time series) | yes | partial |
| A5 | Perrin 2018 eLife (siCLOCK, human myotubes, 48 h RNA-seq) | human | primary myotubes | 994 rhythmic, 1,330 siCLOCK-affected genes | 1,330 / 994 | eLife 34114, PMC5902165; GEO is probably GSE108539/GSE106800 = FORBIDDEN here; gene lists are in eLife supplements (check supp files only, not expression) | yes | in vitro, human |
| A6 | Zhang 2019 eLife (MYOD1 clock amplifier) | mouse | adult muscle | MYOD1 ChIP; circadian direct targets | 536 circadian genes MYOD1 targets | elifesciences.org/articles/43017 (supp not read: page failed) | yes | YES but MYOD1 not clock |
| B1 | Koike 2012 Science (BMAL1, CLOCK, NPAS2, PER1/2, CRY1/2 ChIP-seq, liver) | mouse | liver | peaks; BMAL1 ~5,900, CLOCK ~4,600, NPAS2 ~2,300, PER1 ~4,600, PER2 ~7,300, CRY1 ~16,000, CRY2 ~10,000 sites; 1,444 sites shared by 6 factors | per-gene from peak table S2 | PMC3694775 Table S2 (xlsx 25.7 MB master peak lists), S7 exon cycling (1.5 MB); GEO GSE39860 (ChIP GSE39977) | yes | no (liver) |
| B2 | Zhang/Esser BMAL1+CLOCK ChIP-seq muscle | mouse | skeletal muscle | peaks, young/old | not known | GEO GSE143334 (only OmicsDI summary read; samples/peaks unverified) | yes | YES (best ChIP) |
| B3 | Beytebiere 2019 Genes Dev BMAL1 ChIP | mouse | liver, kidney, heart (no muscle) | BMAL1 peaks | - | GEO accession unconfirmed (GSE110604/GSE110603 candidates) | yes | no |
| B4 | GSE186438 muscle GR and BMAL1 ChIP at ZT4 (+BMAL1 KO) | mouse | muscle | single time point | - | omicsdi GSE186438 | yes | yes |
| B5 | Rey 2011 (GSE26602 liver BMAL1 time course), Cho 2012 / Fang 2014 REV-ERB ChIP, Perelis 2015 | mouse | liver, pancreas | not checked (unverified accessions beyond GSE26602) | - | - | yes | no |
| C1 | MSigDB GOBP_CIRCADIAN_REGULATION_OF_GENE_EXPRESSION (70) and GOBP_RHYTHMIC_PROCESS (300), v2024.1 Hs | human | none | curated annotation, not targets, includes clock genes | 70 / 300 | local: reference/c5.go.bp.v2024.1.Hs.symbols.gmt | yes | no |
| C2 | CircaDB, CGDB | mouse | - | not checked | - | - | - | - |

## Muscle clock KO with several ZT
Yes: GSE43071 (Dyar, mKO vs Cre- control, TA and soleus, 6 CT x 3 reps = 72 samples) and Hodge iMS-Bmal1 (gastrocnemius, 6 CT, pooled). Neither ships a ready "lost rhythm" gene table; derive it.

## Local disk (read-only, /mnt/max8tb/data/CoPhaser_data/cancer_circadian/)
- reference/: c5.go.bp, c2.cgp (no circadian set found by grep), h.all hallmark GMTs (human); Homo_sapiens.gene_info.gz.
- mouse_gene_map/Mus_musculus.gene_info.gz: symbol/Entrez/Ensembl only; NO mouse-human orthology table. Human-mouse by symbol case or HomoloGene/Ensembl BioMart (needs download).
- farm_gtex/orthologs/: human_{btaurus,ggallus,oaries,sscrofa}_clock.tsv, tiny clock-gene tables for other species, not mouse.
- manella2021 (GSE159135: liver/lung/WAT, AlbCre vs BLKO, xlsx counts) and jouffe2022_bmal1ko (GSE150102, ZT12 only): KO transcriptomes, not muscle; can define KO-sensitive genes in other tissues as a cross-tissue support.
- muscle_timed/GSE108539 is the forbidden human muscle data. No mouse muscle KO data on disk.

## To derive the list
1. A1: pull GSE43071 processed expression; per muscle JTK/cosinor in Ctrl, require rhythmic in Ctrl and loss (amplitude, differential rhythmicity) in mKO; this gives clock-autonomous rhythmic genes. Check normalized matrix availability.
2. Intersect with B2 (GSE143334 BMAL1/CLOCK promoter peaks, TSS +/-) and/or B1 liver BMAL1 peaks as direct-target filter; add D-box/REV-ERB sets if wanted.
3. Map mouse to human (needs an orthology table, not on disk), then freeze the set before any human scoring.
