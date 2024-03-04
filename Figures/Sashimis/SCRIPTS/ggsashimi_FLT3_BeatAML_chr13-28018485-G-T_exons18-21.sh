python3 ggsashimi.py \
 -b /TSV/FLT3_BeatAML.tsv \
 -c chr13:28015590-28024943 \
 -o /PLOTS/FLT3_chr13-28018485-G-T_exons18-21_BeatAML \
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
 --ann-height 3 \
 -M 4
