# This script allows to download the fastq files from GEO
for link in $(cat fastq_files.txt);
do
	wget -c $link
done
