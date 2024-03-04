python3 ggsashimi.py \
 -b /TSV/TP53_chr17-7675238-T-C.tsv \
 -c chr17:7675053-7676272  \
 -o /PLOTS/TP53_chr17-7675238-T-C_c.376-2A-G\
 -S minus \
 -g gencode.v36.annotation.gtf \
 -F png \
 -O 3 \
 -C 3 \
 -A mean_j \
 --palette palette_2sashimis.txt \
 --shrink \
 --fix-y-scale \
 --base-size 15 \
 --ann-height 2  \
 -M 4
