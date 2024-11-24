params.mc_fasta = false
params.adapters = false
params.megahit_tmp =  false /* enter <path> to specify a folder to store megahit temp data */
params.checkv_scripts = false
params.derep_split_query = false
params.db_dir = false
params.assembly_groups=["assembly", "co-assembly"]
params.download_only = false
params.ictv_vmr = false
mc_index_name = "mc-genomes"

if(!params.db_dir){
    db_dir = params.data_out + "dbs/"
} else {
    db_dir = params.db_dir
}


log.info """\
    ===================================
    MC Spike pipeline
    ===================================
    In
    -----------------------------------
    sample sheet      : ${params.samples}
    reads folder      : ${params.reads}
    reference genomes : ${params.mc_fasta ? params.mc_fasta : "-" }
    adapters          : ${params.adapters ? params.adapters : "-" }
    CheckV scripts    : ${params.checkv_scripts ? params.checkv_scripts : "-" }
    ICTV VMR file     : ${params.ictv_vmr ? params.ictv_vmr : "-" }

    Out
    -----------------------------------
    data folder       : ${params.data_out}
    database folder   : ${db_dir}
    MEGAHIT temp dir  : ${params.megahit_tmp ? params.megahit_tmp : "-" }
    ===================================
    """
    .stripIndent()


//==================== DOWNLOAD DATABASES =============================
process DL_GENOMAD{
    storeDir "${db_dir}"
    cpus 8
    
    input:
    val(db_dir)

    output:
    path("genomad_db")

    script:
    """
    genomad download-database ./
    """
}

process DL_CHECKV{
    storeDir "${db_dir}"
    cpus 8

    input:
    val(db_dir)

    output:
    path("checkv_db/checkv-db-v*")

    script:
    """
    checkv download_database ./checkv_db
    """
}

process DL_ICTV_VMR{
    storeDir "${params.data_out}/data/taxonomy"
    cpus 1

    input:
    val(db_dir)

    output:
    path("VMR_MSL38_v2.xlsx")

    script:
    """
    wget -O "VMR_MSL38_v2.xlsx" https://ictv.global/sites/default/files/VMR/VMR_MSL38_v2.xlsx 
    """
}

process COPY_METADATA{
    publishDir(
        path: { "${params.data_out}/data/" },
        pattern: "*.csv",
        mode: 'copy'
    )
    cpus 1

    input:
    tuple path(sample_meta), path(vlp_counts), path(mc_refnames)

    output:
    tuple path("sample_sheet.csv"), path("mock_community_particles_added.csv"), path("mock_community_ref_to_virus.csv")

    script:
    """
    if [ ! -f "sample_sheet.csv" ]; then
        mv ${sample_meta} "sample_sheet.csv"
    fi

    if [ ! -f "mock_community_particles_added.csv" ]; then
        mv ${vlp_counts} "mock_community_particles_added.csv"
    fi

    if [ ! -f "mock_community_particles_added.csv" ]; then
        mv ${mc_refnames} "mock_community_ref_to_virus.csv"
    fi
    """
}

//==================== READ QUALITY CONTROL ===========================
process READ_QC {
	cpus 8
	tag "${sample_name}"
	
	publishDir(
        path: { "${params.data_out}/data/read_qc/" },
        pattern: "fastqc/*/*.{html,zip}",
        mode: 'copy'
    )
    publishDir(
        path: { "${params.data_out}/data/read_qc/" },
        pattern: "fastp/*.{html,json}",
        mode: 'copy'
    )
    publishDir(
        path: { "${params.data_out}/reads/" },
        pattern: "cleaned/*.fastq.gz",
        mode: 'copy'
    )
    
	input:
	tuple val(sample_name), path(reads_fw), path(reads_rv), val(readname_fw), val(readname_rv)
	
	output:
	tuple val("${sample_name}"), path("cleaned/${sample_name}_R?.fastq.gz"), emit: reads
    path "fastqc/pre/pre-${sample_name}_R*.{html,zip}", emit: fastqc_pre
    path "fastqc/post/post-${sample_name}_R*.{html,zip}", emit: fastqc_post
	path "fastp/${sample_name}_fastp.{html,json}", emit: fastp
	
	script:
	"""
    mkdir -p fastqc/pre
    mkdir fastqc/post
    mkdir fastp
    mkdir cleaned
    mkdir tmp
    
    # Modify permissions to prevent errors in Singularity 
    chmod 777 fastqc
    chmod 777 fastqc/*

    fastqc --dir ./tmp -o ./fastqc/pre/ -f fastq -q ${reads_fw} ${reads_rv} -t ${task.cpus}
    for f in fastqc/pre/*.{html,zip}; do mv \$f fastqc/pre/pre-\$(basename \$f); done;
    
	fastp \
	-i ${reads_fw} \
	-I ${reads_rv} \
	-o cleaned/${readname_fw} \
	-O cleaned/${readname_rv} \
	--detect_adapter_for_pe \
	--correction \
	--cut_tail \
	-h fastp/${sample_name}_fastp.html \
	-j fastp/${sample_name}_fastp.json \
	--dont_overwrite ${params.adapters ? "--adapter_fasta " + params.adapters : "" } \\
	--thread ${task.cpus}

    fastqc -o fastqc/post -f fastq -q cleaned/${sample_name}_R1.fastq.gz cleaned/${sample_name}_R2.fastq.gz -t ${task.cpus}
	for f in fastqc/post/*.{html,zip}; do mv \$f fastqc/post/post-\$(basename \$f); done;
    """
}

process MULTIQC {
    publishDir "${params.data_out}/data/read_qc/multiqc/", mode:'copy'

    input:
    path '*'

    output:
	path 'multiqc_data/*'
	path "multiqc_report.html"

    script:
    """
	mkdir -p fastqc/pre
	mkdir fastqc/post
    mkdir fastp

	mv pre-*.{html,zip} fastqc/pre/
	mv post-*.{html,zip} fastqc/post/
	
    multiqc -d -dd 1 .
    """
}

process READ_STATS {
    publishDir "${params.data_out}/data/read_qc/", mode:'copy'

    input:
    path("reads_pre/*")
    path("reads_post/*")

    output:
    path("readstats.txt")   

    script:
    """
    seqkit stats -abT ./reads_pre/*.gz | awk '{FS=OFS="\t"}NR==1{ print "filter_step", \$0 } NR>1{print "pre-qc", \$0}' > readstats.txt
    seqkit stats -abT ./reads_post/*.gz | awk '{FS=OFS="\t"}NR>1{print "post-qc", \$0}' >> readstats.txt
    """
}

process READ_DEDUP {
    publishDir(
        path: { "${params.data_out}/data/read_qc/duplication_rate" },
        pattern: "${sample_name}_duplication-rate.json",
        mode: 'copy'
    )
    publishDir(
        path: { "${params.data_out}/reads/deduplicated" },
        pattern: "${sample_name}_dedup_R?.fastq.gz",
        mode: 'copy'
    )
    tag "${sample_name}"
    cpus 1

    input:
    tuple val(sample_name), path(reads_fw), path(reads_rv)

    output:
    tuple val(sample_name), path("${sample_name}_dedup_R1.fastq.gz"), path("${sample_name}_dedup_R2.fastq.gz"), emit: reads
    tuple val(sample_name), path("${sample_name}_duplication-rate.json"), emit: dup_rate

    script:
    """
    hts_SuperDeduper \
	-L ${sample_name}_duplication-rate.json \
	-1 ${reads_fw} \
	-2 ${reads_rv} \
	-f ${sample_name}_dedup ${sample_name}_dedup \
	-l 100 \
	-q 1 \
	-a 1 \
	--log_freq 25000
    """
}

process GET_DUP_RATE{
    publishDir(
        path: { "${params.data_out}/data/read_qc/duplication_rate" },
        pattern: "${sample_name}_duplication-rate.txt",
        mode: 'copy'
    )
    tag "${sample_name}"
    cpus 1
    
    input:
    tuple val(sample_name), path(dup_json)

    output:
    tuple val(sample_name), path("${sample_name}_duplication-rate.txt")

    script:
    """
    get_duplication_data.py ${sample_name} ${dup_json}
    """
}

process CONCAT_DUP_RATE{
    publishDir(
        path: { "${params.data_out}/data/read_qc/" },
        pattern: "duplication-rate.txt",
        mode: 'copy'
    )

    input:
    path("*")

    output:
    path("duplication-rate.txt")

    script:
    """
    DUP_FILES=( *_duplication-rate.txt )

    head -n 1 "\${DUP_FILES[0]}" > duplication-rate.txt
    for DUP_FILE in "\${DUP_FILES[@]}"; do 
        tail -n +2 \$DUP_FILE >> duplication-rate.txt
    done;
    """
}


//==================== MC READ MAPPING ================================
process MC_GENOME_STATS {
    publishDir "${params.data_out}/data/sequences", mode:'copy'
    cpus 1

    input:
    path(mc_fasta)

    output: 
    path("*.sequence_stats.txt")

    script:
    """
    fasta_file=${mc_fasta}
    seqkit fx2tab -nlg \${fasta_file} > \${fasta_file%.*}.sequence_stats.txt
    """
}

process BUILD_MC_IDX {
    publishDir "${params.data_out}/data/bowtie_idx/${mc_index_name}", mode:'copy'
    cpus 8

    input:
    path(mc_fasta)

    output: 
    path("${mc_index_name}_index*.bt2")

    script:
    """
    bowtie2-build --threads ${task.cpus} ${mc_fasta} ${mc_index_name}_index
    """
}

process MAP_MC_READS {
    publishDir(
        path: { "${params.data_out}/data/bowtie/${mc_index_name}" },
        pattern: "*.{bam,bai}",
        mode: 'symlink'
    )
    cpus 8
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(read_files), path(mc_index)

    output:
    tuple val(sample_name), path("*.bam"), path("*.bai"), emit: readmap


    script:
    """
    BAM_FILE=${sample_name}_${mc_index_name}.bam
    bowtie2 -p ${task.cpus} -x ${mc_index_name}_index -1 ${read_files[0]} -2 ${read_files[1]} | \
    samtools view -bS | \
    samtools sort - > \$BAM_FILE
    samtools index \$BAM_FILE
    """
}

process MC_READMAP_STATS {
    publishDir(
        path: { "${params.data_out}/data/bowtie/${mc_index_name}" },
        pattern: "*.{depth,idxstats,gc}.txt",
        mode: 'symlink'
    )
    cpus 8
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam_file), path(bam_idx), val(base_sample), val(spike), val(rt_method), val(replicate)

    output:
    path("*.idxstats.txt"), emit: idxstats
    path("*.depth.txt"), emit: depth
    path("*.gc.txt"), emit: gc

    script:
    """
    BAM_FILE=${bam_file}
    DEP_FILE=\${BAM_FILE/.bam/.depth.txt}
    IDX_FILE=\${BAM_FILE/.bam/.idxstats.txt}
    GC_FILE=\${BAM_FILE/.bam/.gc.txt}
    
    samtools depth -aa -o \$DEP_FILE --threads ${task.cpus} \$BAM_FILE
    samtools idxstats \$BAM_FILE --threads ${task.cpus} > \$IDX_FILE                                                 
    samtools view -F 4 \$BAM_FILE | cut -f 1,3,10 | awk '{ OFS = FS = "\t" }{ l=length(\$3); gcr=(gsub(/G/,"",\$3) + gsub(/C/,"",\$3)); cum_len[\$2]+=l; gctot[\$2] += gcr} END{ for (ref in cum_len) { print ref, gctot[ref]/cum_len[ref] } }' > \$GC_FILE

    for f in *.{depth,idxstats,gc}.txt; do
        mv \$f \$f.old
        awk '{ OFS = FS = "\t" }{ print "${sample_name}", "${base_sample}","${spike}","${rt_method}","${replicate}", \$0} ' \$f.old > \$f
        rm \$f.old
    done;
    """
}

process CONCAT_STATS {
    publishDir "${params.data_out}/data/bowtie/${mc_index_name}/", mode: "copy"
    input:
    path("*")
    path("*")
    path("*")

    output:
    path("mc_idxstats.txt"), emit: idxstats
    path("mc_depth.txt"), emit: depth
    path("mc_gcstats.txt"), emit: gc

    script:
    """
    echo -e "sample\tbase_sample\tspike\trt_method\treplicate\tref\tbp\tdepth" > mc_depth.txt
    for DEPTH in *.depth.txt; do
        awk '\$3!="NO"{ print \$0 }' \$DEPTH >> mc_depth.txt
    done; 

    echo -e "sample\tbase_sample\tspike\trt_method\treplicate\tref\tlength\tmapped_reads\tunmapped_reads" > mc_idxstats.txt
    for IDXSTAT in *.idxstats.txt; do 
        cat \$IDXSTAT >> mc_idxstats.txt
    done;

    echo -e "sample\tbase_sample\tspike\trt_method\treplicate\tref\tavg_gc" > mc_gcstats.txt
    for GCSTAT in *.gc.txt; do
        cat \$GCSTAT >> mc_gcstats.txt
    done;
    """
}


//==================== ASSEMBLY AND DEREPLICATION ======================
process ASSEMBLER {
    publishDir "${params.data_out}/data/assemblies/${assembly_group}", mode: "copy", pattern: "*.fasta"
    publishDir "${params.data_out}/data/assemblies/${assembly_group}/logs", mode: "copy", pattern: "*.log"
    cpus 8
    tag "${assembly_name}"

    input:
    tuple val(assembly_name), val(assembly_group), path(readfiles_fw), path(readfiles_rv)

    output:
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}_contigs.fasta"), emit: assembly
    path "${assembly_name}_readfiles.log", emit: log

    script:
    def r1_input = readfiles_fw.join(',')
	def r2_input = readfiles_rv.join(',')
	"""
	megahit -t ${task.cpus} -1 ${r1_input} -2 ${r2_input} -o ${assembly_name} ${params.megahit_tmp ? "--tmp "+params.megahit_tmp : ""}  # Add temp-dir due to error on linux with CIFS
    cp ${assembly_name}/final.contigs.fa ${assembly_name}_contigs.fasta
	echo -e "${readfiles_fw}\n${readfiles_rv}" > ${assembly_name}_readfiles.log
	"""	
}

process CDHIT {
	tag "${assembly_name}"
	cpus 8
    memory 32
	
	input:
	tuple val(assembly_name), val(assembly_group), path(assembly)

	output:
	tuple val(assembly_name), val(assembly_group), path("${assembly_name}.cdhit.fasta")

	script:
	def cdhit_contigs = "${assembly_name}.cdhit.fasta"
	"""
	cd-hit-est -i ${assembly} -o ${cdhit_contigs} -c 0.95 -n 10 -d 0 -M ${task.memory.toMega()} -T ${task.cpus} -aS 0.8
	"""	
}

process BLAST_DB {
	tag "${assembly_name}"
	cpus 1

	input:
	tuple val(assembly_name), val(assembly_group), path(cdhit_contigs)

	output:
	tuple val(assembly_name), val(assembly_group), path("blastdb/${assembly_name}*")

	script:
	def blast_db = "blastdb/${assembly_name}"
	"""
	makeblastdb -in  ${cdhit_contigs} -dbtype nucl -out ${blast_db}
	"""
}

process BLAST_CONTIGS {
	tag "${assembly_name}"
	cpus 8
	
	input:
	tuple val(assembly_name), val(assembly_group), path(cdhit_contigs), path(blast_db)

	output:
	tuple val(assembly_name), val(assembly_group), path("${assembly_name}.blast.tsv")

	script:
	def blast_table = "${assembly_name}.blast.tsv"
	"""
	blastn -query ${cdhit_contigs} -db ${assembly_name} -outfmt '6 std qlen slen' -max_target_seqs 10000 -out ${blast_table} -num_threads ${task.cpus}
	"""
}

process ANI_CALC {
	tag "${assembly_name}"
    cpus 1

	input:
	tuple val(assembly_name), val(assembly_group), path(blast_table)
    path(anicalc)

	output:
	tuple val(assembly_name), val(assembly_group), path("${assembly_name}_ani.tsv")

	script:
	"""
	python ${anicalc} -i ${blast_table} -o ${assembly_name}_ani.tsv
	"""
}

process ANI_FILTER {
	publishDir "${params.data_out}/data/assemblies/${assembly_group}", mode: "copy"
	tag "${assembly_name}"
    cpus 1

	input:
	tuple val(assembly_name), val(assembly_group), path(cdhit_contigs), path(ani_table)
    path(aniclust)

	output:
	tuple val(assembly_name), val(assembly_group), path("${assembly_name}_contigs.derep.fasta")

	script:
	def checkv_clust = "${assembly_name}_clusters.tsv"
	def checkv_list = "${assembly_name}_clusters.list"
	def checkv_contigs = "${assembly_name}_contigs.derep.fasta"
	"""
	python ${aniclust} --fna ${cdhit_contigs} --ani ${ani_table} --out ${checkv_clust} --min_ani 95 --min_tcov 85 --min_qcov 0
	
	cat ${checkv_clust} | cut -f1 > ${checkv_list}
	seqkit grep -f ${checkv_list} ${cdhit_contigs} > ${checkv_contigs}
	"""
}


//==================== ASSEMBLY READ MAPPING ============================
process BUILD_ASSEMBLY_INDEX {
    publishDir "${params.data_out}/data/bowtie_idx/assemblies/${assembly_group}", mode:'symlink'
    tag "${assembly_name}"
    cpus 8

    input:
    tuple val(assembly_name), val(assembly_group), path(assembly)

    output:
    tuple path("${assembly_name}-index*.bt2"), val("${assembly_name}"), val("${assembly_group}")
    

    script:
    """
    bowtie2-build --threads ${task.cpus} ${assembly} ${assembly_name}-index
    """
}

process MAP_ASSEMBLY_READS {
    publishDir "${params.data_out}/data/bowtie/assemblies/${assembly_group}/${assembly_name}", mode:'symlink'
    tag "${sample_name} / ${assembly_name}"
    cpus 8

    input:
    tuple val(assembly_name), val(assembly_group), val(sample_name), path(read_fw), path(read_rv), val(sample_meta), path(index)

    output:
    tuple val(assembly_name), val(assembly_group), val(sample_name), val(sample_meta), path("${sample_name}_vs_${assembly_name}.bam"), path("${sample_name}_vs_${assembly_name}.bam.bai")

    script:
    """
    BAM_FILE=${sample_name}_vs_${assembly_name}.bam
    bowtie2 -p ${task.cpus} -x ${assembly_name}-index -1 ${read_fw} -2 ${read_rv} | \
    samtools view -bS | \
    samtools sort - > \$BAM_FILE
    samtools index \$BAM_FILE
    """
}

process ASSEMBLY_READMAP_STATS{
    publishDir "${params.data_out}/data/bowtie/assemblies/${assembly_group}/${assembly_name}", mode:'symlink'
    tag "${sample_name} / ${assembly_name}"
    cpus 1

    input:
    tuple val(assembly_name), val(assembly_group), val(sample_name), val(sample_meta), path(bam_file), path(bam_idx)

    output:
    tuple val(assembly_name), val(assembly_group), path("${sample_name}_vs_${assembly_name}.idxstats.txt"), emit: idxstats
    tuple val(assembly_name), val(assembly_group), path("${sample_name}_vs_${assembly_name}.coverage.txt"), emit: coverage
    
    script:
    """
    STAT_FILENAME="${sample_name}_vs_${assembly_name}"
    samtools idxstats ${bam_file} > \$STAT_FILENAME.idxstats.txt
    samtools coverage ${bam_file} > \$STAT_FILENAME.coverage.txt

    mv \$STAT_FILENAME.idxstats.txt \$STAT_FILENAME.idxstats.txt.old
    echo -e "assembly\tsample\tbase_sample\tspike\trt_method\treplicate\tref\tlength\tmapped_reads\tunmapped_reads" > \$STAT_FILENAME.idxstats.txt
    awk '{ OFS = FS = "\t" }{ print "${assembly_name}", "${sample_name}", "${sample_meta["base_sample"]}","${sample_meta["spike"]}","${sample_meta["rt_method"]}","${sample_meta["replicate"]}", \$0} ' \$STAT_FILENAME.idxstats.txt.old >> \$STAT_FILENAME.idxstats.txt
    rm \$STAT_FILENAME.idxstats.txt.old
    
    mv \$STAT_FILENAME.coverage.txt \$STAT_FILENAME.coverage.txt.old
    awk '{ OFS = FS = "\t" }NR==1{ print "assembly\tsample\tbase_sample\tspike\trt_method\treplicate", \$0} NR>1{ print "${assembly_name}", "${sample_name}", "${sample_meta["base_sample"]}","${sample_meta["spike"]}","${sample_meta["rt_method"]}","${sample_meta["replicate"]}", \$0} '  \$STAT_FILENAME.coverage.txt.old > \$STAT_FILENAME.coverage.txt
    rm \$STAT_FILENAME.coverage.txt.old
    """
}

process CONCAT_ASSEMBLY_READMAP_STATS {
    publishDir "${params.data_out}/data/bowtie/assemblies/", mode: "copy"
    cpus 1
    tag "${assembly_group}"

    input:
    tuple val(assembly_group), path("*")
    tuple val(assembly_group), path("*")

    output:
    path("${assembly_group}_collated_idxstats.txt"), emit: idxstats
    path("${assembly_group}_collated_coverage.txt"), emit: coverage

    script:
    """
    IDXSTAT_FILES=( *.idxstats.txt )
    COVSTAT_FILES=( *.coverage.txt )

    head -n 1 "\${IDXSTAT_FILES[0]}" > ${assembly_group}_collated_idxstats.txt
    for IDXSTAT in "\${IDXSTAT_FILES[@]}"; do 
        tail -n +2 \$IDXSTAT >> ${assembly_group}_collated_idxstats.txt
    done;

    head -n 1 "\${COVSTAT_FILES[0]}" > ${assembly_group}_collated_coverage.txt
    for COVSTAT in "\${COVSTAT_FILES[@]}"; do 
        tail -n +2 \$COVSTAT >> ${assembly_group}_collated_coverage.txt
    done;
    """
}


//==================== ASSEMBLY QUALITY CONTROL =========================
process ASSEMBLY_STATS {
    publishDir "${params.data_out}/data/assemblies", mode: "copy"
    tag "${assembly_group}"

    input:
    tuple val(assembly_group), path("*")

    output:
    path("assembly-stats_${assembly_group}.txt")

    script:
    """
    seqkit stats -aT *_contigs.fasta | awk '{FS=OFS="\t"} NR==1{ print "stage", \$0} NR>1{ print "raw", \$0}' > assembly_stats_raw.txt
    seqkit stats -aT *_contigs.derep.fasta | awk '{FS=OFS="\t"} NR==1{ print "stage", \$0} NR>1{ print "dereplicated", \$0}' > assembly_stats_derep.txt
    seqkit stats -aT *_contigs.vir.fasta | awk '{FS=OFS="\t"} NR==1{ print "stage", \$0} NR>1{ print "viral", \$0}' > assembly_stats_vir.txt
    
    head -n 1 assembly_stats_raw.txt > "assembly-stats_${assembly_group}.txt"

    tail -q -n +2 assembly_stats_*.txt | awk '{FS=OFS="\t"}{gsub("_contigs.*.fasta", "", \$2); print }' >> "assembly-stats_${assembly_group}.txt"
    """
}

process CONTIG_STATS{
    tag "${assembly_name}"

    input:
    tuple val(assembly_name), val(assembly_group), path(assembly)

    output:
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}_contigs-stats-single.txt")

    script:
    """
    # sort contigs by length, long to short
    seqkit sort -lr ${assembly} | \
    # get contig id, length and gc content
    seqkit fx2tab --length --gc --name --only-id | \
    # print sample name and seqkit output
    awk '{FS=OFS="\t"} {print "${assembly_name}", \$0}' > ${assembly_name}_contigs-stats-single.txt
    """
}

process CONCAT_CONTIG_STATS{
    publishDir "${params.data_out}/data/assemblies", mode: "copy"
    tag "${assembly_group}"

    input:
    tuple val(assembly_group), path(assembly)

    output:
    tuple val(assembly_group), path("contigs-stats_${assembly_group}.txt")

    script:
    """
    echo -e "assembly_name\tcontig\tlength\tgc" > contigs-stats_${assembly_group}.txt
    cat *_contigs-stats-single.txt >> contigs-stats_${assembly_group}.txt
    """
}

process QUAST {
	publishDir "${params.data_out}/data/quast/${assembly_group}", mode:'copy'
	tag "${assembly_group}"
    cpus 8

	input:
    tuple val(assembly_group), path(assemblies)
    path(mc_fasta)
	
	output:
	tuple val(assembly_group), path("quast_results/latest/*"), emit: quast_total
    tuple val(assembly_group), path("quast_results/latest/combined_reference/transposed_report.tsv"), emit: quast_stats
		
	script:
	"""
    seqkit split -i ${mc_fasta}; 
    
    mc_fasta_name=\$(basename ${mc_fasta})
    mc_fasta_prefix=\${mc_fasta_name%.*}.part_
    
    function join_by {
    local d=\${1-} f=\${2-}
    if shift 2; then
        printf %s "\$f" "\${@/#/\$d}"
    fi
    }

    for f in \${mc_fasta_name}.split/*.fasta; do 
        fname=\$(basename \$f);  
        mv \$f \${mc_fasta_name}.split/\${fname#\$mc_fasta_prefix};
    done

    ref_files=( \${mc_fasta_name}.split/*.fasta )
    refs=\$( join_by , \${ref_files[@]} )

	metaquast.py ${assemblies} -r \${refs} -t ${task.cpus}
	"""
}

//==================== VIRUS SEQUENCE IDENTIFICATION ====================
process GENOMAD {
    publishDir "${params.data_out}/data/virus_classification/${assembly_group}/genomad", mode: 'copy'
    cpus 8
    tag "${assembly_name}"

    input:
    tuple val(assembly_name), val(assembly_group), path(assembly)
    path(genomad_db)


    output:
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}_virus_summary.txt"), emit: genomad_summary
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}/*"), emit: genomad_total

    script:
    """
    genomad end-to-end --min-score 0.8 --threads ${task.cpus} \
	    ${assembly} \
	    ${assembly_name} \
	    ${genomad_db}
    
    assembly_filename=\$(basename ${assembly})
    assembly_filename=\${assembly_filename%.*}
    
    summary_file=${assembly_name}/\${assembly_filename}_summary/\${assembly_filename}_virus_summary.tsv
    
    awk '{FS=OFS="\t"}NR==1{ print "assembly_name", \$0} NR>1{ print "${assembly_name}", \$0}' \$summary_file > "${assembly_name}_virus_summary.txt"
    """
}

process CONCAT_GENOMAD{
    publishDir "${params.data_out}/data/virus_classification", mode: 'copy'
    cpus 1
    tag "${assembly_group}"

    input:
    tuple val(assembly_group), path(genomad_summaries)

    output:
    tuple val(assembly_group), path("genomad_summary_${assembly_group}.txt")

    script:
    """
    FILES=( *_virus_summary.txt )

    head -n 1 "\${FILES[0]}" > genomad_summary_${assembly_group}.txt
    for FILE in "\${FILES[@]}"; do 
        tail -n +2 \$FILE >> genomad_summary_${assembly_group}.txt
    done;
    """
}

process MC_BLAST {
    publishDir "${params.data_out}/data/alignments/mc-genomes/${assembly_group}", pattern: "*_blast_mc-genomes.txt", mode: 'copy'
    publishDir "${params.data_out}/data/virus_classification/${assembly_group}/mc-genomes", pattern: "*_ani_mc-genomes.txt", mode: 'copy'
    cpus 1
    tag "${assembly_name}"

    input:
    tuple val(assembly_name), val(assembly_group), path(assembly)
    path(mc_fasta)
    path(anicalc)

    output:
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}_blast_mc-genomes.txt"), emit: blast
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}_ani_mc-genomes.txt"), emit: ani

    script:
    def blast_table="${assembly_name}_blast_mc-genomes.txt"
    def ani_table="${assembly_name}_ani_mc-genomes.txt"
    """
    blastn \
	-query ${assembly} \
	-subject ${mc_fasta} \
	-outfmt "6 std qlen slen" \
	-max_target_seqs 10000 \
	-out  ${blast_table}

    if [[ -s ${blast_table } ]]; then
        python ${anicalc} -i ${blast_table} -o ${ani_table}.tmp
    else
        touch ${ani_table}.tmp
    fi

    
    awk '
    BEGIN{
        FS = OFS = "\t"; 
        print "assembly_name\tcontig\tref\tnum_alns\tani\tqcov\ttcov"
        } 
    {
    if(NR > 1 && \$4 >=98 && \$5>=95){ 
        print "${assembly_name}", \$0 
        }
    }
    ' ${ani_table}.tmp > ${ani_table}

    """
}

process CONCAT_MC_BLAST {
    publishDir "${params.data_out}/data/virus_classification/", mode: 'copy'
    cpus 1
    tag "${assembly_group}"

    input:
    tuple val(assembly_group), path(ani_table)

    output:
    tuple val(assembly_group), path("mc-genome_ani_${assembly_group}.txt")

    script:
    """
    ANI_FILES=( *_ani_mc-genomes.txt )

    head -n 1 "\${ANI_FILES[0]}" > mc-genome_ani_${assembly_group}.txt
    for FILE in "\${ANI_FILES[@]}"; do 
        tail -n +2 \$FILE >> mc-genome_ani_${assembly_group}.txt
    done;    
    """

}

process FILTER_VIRAL{
    publishDir "${params.data_out}/data/virus_classification/${assembly_group}", pattern: "*.txt", mode: 'copy'
    publishDir "${params.data_out}/data/assemblies/${assembly_group}", pattern: "*.fasta", mode: 'copy'
    cpus 1
    tag "${assembly_name}"

    input:
    tuple val(assembly_name), val(assembly_group), path(genomad), path(mc_ani), path(assembly)

    output:
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}_contigs.vir.fasta"), emit: assembly
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}_vir-contig-list.txt"), emit: list

    script:
    """
    cat <(awk 'NR>1' ${genomad} | cut -f2) <(awk 'NR>1' ${mc_ani} | cut -f2) | sort | uniq > ${assembly_name}_vir-contig-list.txt
    seqkit grep -f ${assembly_name}_vir-contig-list.txt ${assembly} > ${assembly_name}_contigs.vir.fasta
    """
}

process CHECKV {
    publishDir "${params.data_out}/data/virus_classification/${assembly_group}/checkv/${assembly_name}", mode: 'copy'
    cpus 8
    tag "${assembly_name}"

    input:
    tuple val(assembly_name), val(assembly_group), path(assembly)
    path(checkv_db)


    output:
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}_checkv_summary.txt"), emit: checkv_summary
    tuple val(assembly_name), val(assembly_group), path("${assembly_name}/*"), emit: checkv_total

    script:
    """
    checkv end_to_end ${assembly} ${assembly_name} -t ${task.cpus} -d ${checkv_db}

    assembly_filename=\$(basename ${assembly})
    assembly_filename=\${assembly_filename%.*}
    
    summary_file=${assembly_name}/quality_summary.tsv
    
    awk '{FS=OFS="\t"}NR==1{ print "assembly_name", \$0} NR>1{ print "${assembly_name}", \$0}' \$summary_file > "${assembly_name}_checkv_summary.txt"

    """
}

process CONCAT_CHECKV{
    publishDir "${params.data_out}/data/virus_classification", mode: 'copy'
    cpus 1
    tag "${assembly_group}"

    input:
    tuple val(assembly_group), path(checkv_summaries)

    output:
    tuple val(assembly_group), path("checkv_summary_${assembly_group}.txt")

    script:
    """
    FILES=( *_checkv_summary.txt )

    head -n 1 "\${FILES[0]}" > checkv_summary_${assembly_group}.txt
    for FILE in "\${FILES[@]}"; do 
        tail -n +2 \$FILE >> checkv_summary_${assembly_group}.txt
    done;
    """
}

//==================== VIRUS TAXONOMY ====================================
process CREATE_TAX_TABLE{
    publishDir "${params.data_out}/data/taxonomy", mode: "copy"
    cpus 1

    input:
    path(ictv_vmr)

    output:
    path("ICTV-Clade-to-Rank.txt")

    script:
    """
    translate_ictv_clade_to_rank.py ${ictv_vmr} ICTV-Clade-to-Rank.txt
    """
}

process ADD_TAX_LINEAGE{
    publishDir "${params.data_out}/data/taxonomy", mode: "copy"
    tag "${assembly_group}"
    cpus 1

    input:
    tuple val(assembly_group), path(genomad_summary)
    path(ictv_clade_ranks)

    output:
    tuple val(assembly_group), path("genomad_summary_taxonomy_${assembly_group}.txt")

    script:
    """
    get_taxonomic_lineage.py ${genomad_summary} ${ictv_clade_ranks} genomad_summary_taxonomy_${assembly_group}.txt
    """
}

workflow download_dbs{
    take:
    db_dir

    main:
    DL_GENOMAD(db_dir)
    DL_CHECKV(db_dir)
    
    if(!params.ictv_vmr){
        ictv_vmr = DL_ICTV_VMR(db_dir)
    } else{
        ictv_vmr = channel.fromPath(params.ictv_vmr)
    }
    
    emit:
    genomad = DL_GENOMAD.out
    checkv = DL_CHECKV.out
    ictv_vmr = ictv_vmr
}

workflow quality_filter  {
    take:
    sample_metadata

    main:
    read_qc_data = READ_QC(
        sample_metadata.map( it -> [it["sample_name"], params.reads + "/" + it["fw_readfile"], params.reads + "/" + it["rv_readfile"], it["fw_readfile"], it["rv_readfile"]] )
        )

    MULTIQC(read_qc_data.fastp.concat(read_qc_data.fastqc_pre.concat(read_qc_data.fastqc_post)).collect())
    
    READ_STATS(
        sample_metadata
        .map( it -> [params.reads + "/" + it["fw_readfile"], params.reads + "/" + it["rv_readfile"]])
        .collect(),
        read_qc_data.reads
        .map( it -> it[1])
        .collect()
    )

    dedup = READ_DEDUP(sample_metadata.map( it -> [it["sample_name"], params.reads + "/" + it["fw_readfile"], params.reads + "/" + it["rv_readfile"]] ))
    dup_rate_files = GET_DUP_RATE(dedup.dup_rate)
    dup_rate = CONCAT_DUP_RATE(dup_rate_files.map{ it[1] }.collect())


    emit:
    reads = read_qc_data.reads
    dup_rate = dup_rate
}

workflow mc_readmap {
    take:
    sample_metadata
    reads

    main:
        mc_genome_stats = MC_GENOME_STATS(params.mc_fasta)
        mc_index = BUILD_MC_IDX(params.mc_fasta)
        mc_readmaps = MAP_MC_READS(reads.combine(mc_index.toList()))
        
        sample_metadata
        .map( it -> [it["sample_name"], it["base_sample"], it["spike"], it["rt_method"], it["replicate"]])
        .set { sample_metadata }
        
        mc_readmap_statfiles = MC_READMAP_STATS(mc_readmaps.readmap.join(sample_metadata))
        mc_readmap_stats = CONCAT_STATS(mc_readmap_statfiles.idxstats.collect(), mc_readmap_statfiles.depth.collect(), mc_readmap_statfiles.gc.collect())

    emit:
    mc_readmap_stats.idxstats
    mc_readmap_stats.depth
    mc_readmap_stats.gc
}

workflow coassembly {
    take:
    sample_metadata
    reads

    main:
       
    assemblies = ASSEMBLER(reads) // out: assembly_name, assembly_group, assembly_fasta
    
    
    cdhit_contigs = CDHIT(assemblies.assembly)
	
	blast_dbs = BLAST_DB(cdhit_contigs)
    blast_tables = BLAST_CONTIGS(cdhit_contigs.join(blast_dbs, by: [0, 1]))
    
    ani_tables = ANI_CALC(blast_tables, params.checkv_scripts + "/anicalc.py")
	derep = ANI_FILTER(cdhit_contigs.join(ani_tables, by: [0,1]), params.checkv_scripts + "/aniclust.py")

	contig_collection = assemblies.assembly
	.concat(
		derep
		.map( it -> it[2])
		)
	.collect()

    emit:
    raw = assemblies.assembly 
    derep = derep

}

workflow assembly_readmap{
    take:
    reads_to_map
    assemblies

    main:
    assembly_idxs = BUILD_ASSEMBLY_INDEX(assemblies)
    
    readmaps = MAP_ASSEMBLY_READS(reads_to_map.combine(assembly_idxs, by: [1,2]))

    readmap_stats = ASSEMBLY_READMAP_STATS(readmaps)
    
    idxstats_to_collate = readmap_stats.idxstats
    .map( it -> it[1..-1])
    .groupTuple()

    coverage_to_collate = readmap_stats.coverage
    .map( it -> it[1..-1])
    .groupTuple()

    collated_readstats = CONCAT_ASSEMBLY_READMAP_STATS(idxstats_to_collate, coverage_to_collate)
    
    emit:
    idxstats=collated_readstats.idxstats
    coverage=collated_readstats.coverage
}

workflow virus_id {
    take:
    assemblies
    genomad_db
    checkv_db

    main:
    genomad = GENOMAD(assemblies, genomad_db)

    blast = MC_BLAST(assemblies, params.mc_fasta, params.checkv_scripts + "/anicalc.py")

    vir_list = genomad.genomad_summary
    .join(blast.ani, by: [0, 1])
    .join(assemblies, by: [0, 1])

    vir_contigs = FILTER_VIRAL(vir_list)

    genomad_summary_single = genomad.genomad_summary
    .map{ [it[1], it[2]] }
    .groupTuple()

    genomad_summary = CONCAT_GENOMAD(genomad_summary_single)
    
    blast_ani_single = blast.ani
    .map{ [it[1], it[2]] }
    .groupTuple()

    blast_ani = CONCAT_MC_BLAST(blast_ani_single)

    checkv = CHECKV(vir_contigs.assembly, checkv_db)

    checkv_summary_single = checkv.checkv_summary
    .map{ [it[1], it[2]] }
    .groupTuple()
    
    checkv_summary = CONCAT_CHECKV(checkv_summary_single)

    emit:
    genomad = genomad_summary
    blast = blast_ani
    checkv = checkv_summary
    contigs = vir_contigs.assembly
}

workflow assembly_statistics{
    take:
    assemblies
    dereplicated_assemblies
    virus_contigs

    main:
    // COLLECT ASSEMBLY STATISTICS
    assembly_collection = assemblies
        .concat(dereplicated_assemblies)
        .concat(virus_contigs)
        .map( it -> [it[1], it[2]])
        .groupTuple()
    
    assembly_stats = ASSEMBLY_STATS(assembly_collection)

    contig_stats_single = CONTIG_STATS(assemblies)
    .map { [it[1], it[2]]}
    .groupTuple()
    
    contig_stats = CONCAT_CONTIG_STATS(contig_stats_single)
    
    quast = QUAST(assembly_collection, params.mc_fasta)

    emit:
    assembly_stats
    contig_stats
    quast.quast_stats
}

workflow virus_taxonomy{
    take:
    ictv_vmr
    genomad_summary

    main:
    clade_to_ranks = CREATE_TAX_TABLE(ictv_vmr)

    ADD_TAX_LINEAGE(genomad_summary, clade_to_ranks)
}

workflow {
    sample_metadata = Channel.fromPath( "${params.samples}" )
	    .splitCsv(header: true)

    // DOWNLOAD DATABASES
    download_dbs(db_dir)    

    if(!params.download_only){
        // READ QUALITY CONTROL
        clean_reads = quality_filter(sample_metadata).reads
        
        // MAP READS TO MC REFERENCES
        mc_readmap_stats = mc_readmap(sample_metadata, clean_reads)
        
        // ORGANISE READS BY ASSEMBLY GROUP(S)
        assembly_groups = params.assembly_groups

        read_groups = sample_metadata
        .map( it -> [it["sample_name"],  it[assembly_groups[0]], assembly_groups[0]])
        
        if(assembly_groups.size > 1){
            for(assembly_group in assembly_groups[1..-1]){
                read_groups = sample_metadata
                .map( it -> [it["sample_name"], it[assembly_group], assembly_group])
                .filter { it[1] != "" && it[1] != "-"}
                .concat(read_groups)
            }
        }

        read_groups
        .combine(clean_reads.map( it -> [it[0], it[1][0], it[1][1]]), by:0)
        .set { read_groups }    

        reads_by_assembly = read_groups
        .map(it -> it[1..-1])
        .groupTuple(by: [0,1])
        
        // ASSEMBLE READS
        assemblies = coassembly(sample_metadata, reads_by_assembly)
        
        // MAP READS TO ASSEMBLY
        reads_to_map = read_groups
        .combine(
            sample_metadata
            .map( it -> it.findAll{!(assembly_groups+["fw_readfile", "rv_readfile"]).contains(it.key)} )
            .map( it -> [it["sample_name"], it]),
            by: 0
        )

        assembly_readmap_stats = assembly_readmap(reads_to_map, assemblies.derep)
        
        // IDENTIFY VIRAL SEQUENCES
        virus_contigs = virus_id(assemblies.derep, download_dbs.out.genomad, download_dbs.out.checkv)
        
        // COLLECT ASSEMBLY STATISTICS
        assembly_stats = assembly_statistics(assemblies.raw, assemblies.derep, virus_contigs.contigs)

        // ADD VIRUS TAXONOMY
        
        virus_taxonomy(download_dbs.out.ictv_vmr, virus_id.out.genomad)

        COPY_METADATA([params.samples, params.vlp_counts, params.mc_refnames])
    }
}