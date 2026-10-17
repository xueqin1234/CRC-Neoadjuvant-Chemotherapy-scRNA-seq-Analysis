#==================== CellPhoneDB shell commands (run in terminal) ====================
source activate cpdb_laoban
cpdb_database="/hanlab/xueqin/cellphoneDB_to_CRC/cellphoneDB_laoban/cellphonedb-data-4.0.0/cellphone.db"

# CRC0
cd /hanlab/xueqinliu/test/cellphoneDB_to_CRC/cellphoneDB_laoban/Intermediate/CRC0
nohup cellphonedb method statistical_analysis Intermediate_CRC_CRC0_stromal_meta.txt Intermediate_fit_CRC0.h5ad --database ${cpdb_database} --counts-data=gene_name --threads=24 &

# CRC2
cd /hanlab/xueqinliu/test/cellphoneDB_to_CRC/cellphoneDB_laoban/Intermediate/CRC2
nohup cellphonedb method statistical_analysis Intermediate_CRC_CRC2_stromal_meta.txt Intermediate_fit_CRC2.h5ad --database ${cpdb_database} --counts-data=gene_name --threads=24 &

# CRC4
cd /hanlab/xueqinliu/test/cellphoneDB_to_CRC/cellphoneDB_laoban/Intermediate/CRC4
nohup cellphonedb method statistical_analysis Intermediate_CRC_CRC4_stromal_meta.txt Intermediate_fit_CRC4.h5ad --database ${cpdb_database} --counts-data=gene_name --threads=24 &