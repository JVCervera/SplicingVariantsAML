python3 ggsashimi.py \
 -b /TSV/PTPN11.tsv \
 -c chr12:112489024-112502256 \
 -o /PLOTS/PTPN11_chr12-112489084-G-T_exons13-14 \
 -S plus \
 -g gencode.v36.annotation.gtf \
 -F png \
 -O 3 \
 -C 3 \
 -A mean_j \
 --palette palette_3sashimis.txt \
 --shrink \
 --fix-y-scale \
 --base-size 15 \
 --ann-height 3 \
 -M 4
