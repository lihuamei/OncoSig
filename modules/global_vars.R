PAM50_COLOR <- list(
    'Normal' = '#03018B',
    'LumA' = '#C0F7FF',
    'LumB' = '#99C0CC',
    'Her2' = '#FF68B5',
    'Basal' = '#CC00FF'
)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

PAIRWISE_COLORS <- list(
	'MT' = '#F9C996',
	'PT' = '#90B4C6'
)

CANCER_COLORS <- list(
  `ACC` = '#C6A221',`BLCA`= '#F8CED9',`BRCA`= '#E62D8A',
  `CESC`= '#F7B263',`CHOL`= '#00467B',`COAD`= '#81CDE9',
  `DLBC`= '#254D9C',`ESCA`= '#0A78AC',`GBM`= '#BD4E95',
  `HNSC`= '#76C4A1',`KICH`= '#E72728',`KIRC`= '#F4ACB1',
  `KIRP`= '#EC6D72',`LAML`= '#704668',`LGG`= '#DF9AC2',
  `LIHC`= '#C6CCDC',`LUAD`= '#D7C2DE',`LUSC` = '#A381B9',
  `MESO`= '#592685',`OV` = '#E87C1A',`PAAD` = '#6679A1',
  `PCPG`= '#E5C11B', `PRAD`= '#8B1C21',`READ` ='#D0EBF9' ,
  `SARC` = '#12A496',`SKCM` = '#84C36D',`STAD`= '#22A4DF',
  `TGCT` = '#D01F29',`THCA`= '#F4E628', `THYM`= '#D6AA8E',
  `UCEC`= '#FCE1C6', `UCS`= '#F29119',`UVM`= '#158F3B'
)



