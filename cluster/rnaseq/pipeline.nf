#!/usr/bin/env nextflow

/* Directories */
cwd = "/grid/mccandlish/home_norepl/martigo/projects/22_drugs/cluster/rnaseq"
ref_dir = "$cwd/reference"
bam_dir = "$cwd/bam"
fastq_dir = "$cwd/fastq"
rmats_dir = "$cwd/rmats"
counts_dir = "$cwd/counts"
qc_dir = "$cwd/qc"

/* Files */
data = "data.csv"
design = "$cwd/design.csv"
design_pred = "$cwd/design_pred.csv"
psi = "psi"

/* Other settings */
read_length = 300
overhang = 8
bam_suffix = "Aligned.sortedByCoord.out.bam"

/* Channels */
Channel
    .fromPath(data)
    .splitCsv(header:true)
    .map{ row-> row.Run}
    .into { samples_to_map ; samples_fastqc }


/* Reference preparation */

process download_genome {

   output:
   val "$ref_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" into fasta_gz_ch
   
   """
   wget -c https://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P "$ref_dir"
   """
}


process unzip_genome {

   conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/samtools'

   input:
   val fagz from fasta_gz_ch

   output:
   val "$ref_dir/genome.fa" into fastafile_ch
   
   """
   gunzip -c $fagz >"$ref_dir/genome.fa" 
   samtools faidx "$ref_dir/genome.fa"
   """
}


fastafile_ch.into{ fastafile_ch1 ; fastafile_ch2}


process download_annotation {

   output:
   val "$ref_dir/Homo_sapiens.GRCh38.103.gtf.gz" into gtf_gz_ch
   
   """
   wget -c https://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz -P "$ref_dir"
   """
}

process unzip_annotation {

   input: 
   val gtf_gz from gtf_gz_ch
   
   output:
   val "$ref_dir/annotation.gtf" into gtf_ch
   
   """
   gunzip -c $gtf_gz > "$ref_dir/annotation.gtf" 
   """
}

gtf_ch.into{ gtf_ch1; gtf_ch2 ; gtf_ch3}


process gff_to_bed{

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/gffread'

    input:
    val gtf from gtf_ch1

    output:
    val "$ref_dir/annotation.bed" into bed_ch

    """
    gffread $gtf --bed -M -K -Y --no-pseudo -U -l 500 > "$ref_dir/annotation.bed"
    """
}


process build_star_index{

    cpus 2
    memory '42 GB'
    executor 'sge'

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/star'

    input:
    path fasta from fastafile_ch2
    path gtf from gtf_ch2

    output:
    val "$ref_dir/index" into index_ch

    """
    mkdir star_index
    STAR --runThreadN 2 --runMode genomeGenerate --genomeDir $ref_dir/index\
         --genomeFastaFiles $fasta --sjdbGTFfile $gtf\
         --sjdbOverhang $read_length\
    """
}





/* Running analysis */

process star_map{

    cpus 4
    memory '42 GB'

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/star'

    input:
    val sample from samples_to_map
    val index from index_ch

    output:
    val "${bam_dir}/${sample}" into bams_ch

    """
    STAR --runThreadN 4 --genomeDir $index --readFilesCommand zcat\
         --readFilesIn "${fastq_dir}/${sample}_1.fastq.gz" "${fastq_dir}/${sample}_2.fastq.gz" \
         --alignSJoverhangMin $overhang --alignSJDBoverhangMin $overhang \
         --outFileNamePrefix "${bam_dir}/${sample}" --outSAMtype BAM SortedByCoordinate
    """
}

process index_bam{

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/samtools'

    input:
    val bam_prefix from bams_ch

    output:
    val "${bam_prefix}" into indexed_bam

    """
    samtools index "${bam_prefix}${bam_suffix}"
    """
}

indexed_bam.into{ rmats_bams ; frag_size_input ; bam_stat_input ; 
                  gene_body_coverage_input; jc_saturation_input ;
                  dups_input ; read_distribution_input;
                  feat_counts_input }

rmats_input = rmats_bams.collect()


process get_control_bamnames{

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/as_quant'

    output:
    val "$cwd/control.txt" into control_bams

    """
    python "$cwd/get_bam_names.py" "control"
    """
}


process get_treated_bamnames{

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/as_quant'

    output:
    val "$cwd/treated.txt" into treated_bams

    """
    python "$cwd/get_bam_names.py" "treated"
    """
}


process run_rmats {

    memory '32 GB'

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/rmats2'

    input:
    val bams from rmats_input
    val control from control_bams
    val treated from treated_bams
    val gtf from gtf_ch3

    output:
    val "${rmats_dir}/SE.MATS.JC.txt" into rmats_output

    """
    rmats.py --gtf $gtf --b1 $control --b2 $treated --novelSS\
             -t paired --libType fr-unstranded --od $rmats_dir\
             --readLength $read_length\
             --statoff --variable-read-length --tmp tmp
    """
}

process get_sample_names {

    output:
    path "sample_names" into sample_names_ch

    """
    cut -f 1 -d',' $design | grep -v "sample" > "sample_names"
    """

}


process parse_exon_skipping {

    memory '4 GB'

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/as_quant'

    input:
    val rmats_out from rmats_output
    path sample_names from sample_names_ch

    output:
    val "${counts_dir}/exon_skipping" into exon_skipping_counts_ch

    """
    parse_counts $rmats_out -p "rmats" -s "$sample_names" -e "exon_skipping" -o "${counts_dir}/exon_skipping"
    """
}


process filter_exon_counts {

    memory '4 GB'

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/as_quant'

    input:
    val counts from exon_skipping_counts_ch

    output:
    val "${counts}.filtered" into filtered_exon_counts_ch

    """
    filter_exon_counts -i $counts -T 10 -o "${counts}.filtered" --inc_suffix "inclusion" --skp_suffix "skipping" -s 1
    """
}

filtered_exon_counts_ch.into{ filtered_exon_counts_ch1 ; filtered_exon_counts_ch2; filtered_exon_counts_ch3}

process infer_exon_psis {

    cpus 4
    memory '8 GB'

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/as_quant'

    input:
    val counts from filtered_exon_counts_ch1

    """   
    fit_psi_regression "${counts}.inclusion.csv" "${counts}.total.csv" -d $design --design_pred $design_pred --log_bias "${counts}.log_bias.csv" -n 4 -o "$cwd/exons"
    """ 
}


process get_5ss_sequences {

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/as_quant'

    input:
    path fasta from fastafile_ch1
    val counts from filtered_exon_counts_ch3
    
    output:
    path "exons.5ss.csv" into exon_5ss_seqs_ch
    
    """
    python "$cwd/get_5ss_seqs.py" $fasta "${counts}.total.csv" "$exons.5ss.csv"
    """

}

process infer_allelic_manifolds{

    cpus 4
    memory '8 GB'

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/as_quant'

    input:
    val counts from filtered_exon_counts_ch3
    val exon_5ss_seqs from exon_5ss_seqs_ch
    
    """
    python "$cwd/fit_5ss_model.py" $counts $exon_5ss_seqs $design "$cwd/allelic_manifolds.csv"
    """
}


/* QC metrics */

process fastqc{

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/fastqc'

    input:
    val sample from samples_fastqc

    """
    echo "redo"
    fastqc "${fastq_dir}/${sample}_1.fastq.gz" "${fastq_dir}/${sample}_2.fastq.gz" -o "${qc_dir}"
    """
}

bed_ch.into {bed1 ; bed2 ; bed3 ; bed4;}

process bam_stat {

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/rseqc'

    input:
    val bam_prefix from bam_stat_input

    """
    bam_stat.py  -i "${bam_prefix}${bam_suffix}" > "${bam_prefix}.bam_stats.txt" 
    """
}


process fragment_size {

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/rseqc'

    input:
    val bam_prefix from frag_size_input
    val bed from bed1

    """
    inner_distance.py  -i "${bam_prefix}${bam_suffix}" -r $bed -o "${bam_prefix}.frag_size" -s 1 -l -200 -u 500
    """
}


process gene_body_coverage {

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/rseqc'

    input:
    val bam_prefix from gene_body_coverage_input
    val bed from bed2

    """
    geneBody_coverage.py -i "${bam_prefix}${bam_suffix}" -r $bed -o "${bam_prefix}.gene_body_coverage"
    """
}

process junction_saturation {

    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/rseqc'

    input:
    val bam_prefix from jc_saturation_input
    val bed from bed3

    """
    junction_saturation.py -i "${bam_prefix}${bam_suffix}" -r $bed -o "${bam_prefix}.junction_saturation"
    """
}

process read_duplication {
    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/rseqc'

    input:
    val bam_prefix from dups_input

    """
    read_duplication.py -i "${bam_prefix}${bam_suffix}" -o "${bam_prefix}.read_duplication"
    """
}

process read_distribution {
    conda '/grid/mccandlish/home_norepl/martigo/miniconda3/envs/rseqc'

    input:
    val bam_prefix from read_distribution_input
    val bed from bed4

    """
    read_distribution.py -i "${bam_prefix}${bam_suffix}" -r $bed > "${bam_prefix}.read_distribution.txt"
    """
}



