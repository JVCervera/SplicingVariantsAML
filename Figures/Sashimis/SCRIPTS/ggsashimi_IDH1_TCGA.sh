python3 ggsashimi.py \
 -b /TSV/IDH1_TCGA.tsv \
 -c chr2:208245319-208248660 \
 -o /PLOTS/IDH1_TCGA\
 -S minus \
 -g gencode.v36.annotation.gtf \
 -F png \
 -O 3 \
 -C 3 \
 -A mean_j \
 --palette palette_U2AF1.txt \
 --shrink \
 --base-size 15 \
 --ann-height 2  \
 -M 1
