# mapped, cleaned bam files
wget https://www.encodeproject.org/files/ENCFF652XIT/@@download/ENCFF652XIT.bam c2c12_myoblast_1.bam 
wget https://www.encodeproject.org/files/ENCFF408MUF/@@download/ENCFF408MUF.bam c2c12_myoblast_2.bam
wget https://www.encodeproject.org/files/ENCFF202MCY/@@download/ENCFF202MCY.bam c2c12_myotube_1.bam
wget https://www.encodeproject.org/files/ENCFF396UFT/@@download/ENCFF396UFT.bam c2c12_myotube_2.bam

# reference
mkdir ref
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz ref/
gunzip ref/gencode.vM21.annotation.gtf.gz

# make talon config file
touch config.csv
printf "PB154,C2C12 myoblast,SequelII,c2c12_myoblast_1.bam\n" >> talon_config.csv
printf "PB155,C2C12 myoblast,SequelII,c2c12_myoblast_2.bam\n" >> talon_config.csv
printf "PB213,C2C12 myotube,SequelII,c2c12_myotube_1.bam\n" >> talon_config.csv
printf "PB214,C2C12 myotube,SequelII,c2c12_myotube_2.bam\n" >> talon_config.csv

# run talon
talon_initialize_database \
    --f ref/gencode.vM21.annotation.gtf \
    --g mm10 \
    --a vM21 \
    --o c2c12
    
talon \
    --f config.csv \
    --db c2c12.db \
    --o c2c12

talon_filter_transcripts \
    --db c2c12.db \
    -a vM21 \
    --maxFracA 0.5 \
    --minDatasets 2 \
    --minCount 5 \
    --o all_pass_list.csv
    
talon_create_GTF \
    --db c2c12.db \
    -b mm10 \
    -a vM21 \
    --whitelist all_pass_list.csv \
    --o all
    
talon_abundance \
    --db c2c12.db \
    -a vM21 \
    --whitelist all_pass_list.csv \
    -b mm10 \
    --o all
    