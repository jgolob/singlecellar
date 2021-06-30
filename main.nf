#!/usr/bin/env nextflow

/*
  singlecellar: A pipeline to wrap around and run scRNA analysis.

*/

// Using DSL-2
nextflow.enable.dsl=2

// containers
container__cellranger = 'golob/cellranger:6.0.2A'
container__trimgalore = 'quay.io/biocontainers/trim-galore:0.6.6--0'

// Default parameters
params.manifest = false
params.help = false
params.cr_reference = false
params.output = './output'

workflow SingleCellRNA {
    take:
        read_mainfest_ch
        cr_ref_tgz

    main:

    // trimgalore to clean up / trim

    Cellranger_count(
        read_mainfest_ch.map{
            [it.specimen, file(it.R1), file(it.R2)]
        },
        cr_ref_tgz
    )

    CombineCRout(
        Cellranger_count.out
            .collect{ it[1] }
    )


}

process Cellranger_count {
    container "${container__cellranger}"
    label 'multithread'
    errorStrategy 'finish'

    input:
    tuple val(specimen), path(R1), path(R2)
    path cr_ref

    output:
    tuple val(specimen), path("${specimen}_result.tgz")
    """
    mkdir reads/
    mv ${R1} reads/${specimen}_S1_L001_R1_001.fastq.gz
    mv ${R2} reads/${specimen}_S1_L001_R2_001.fastq.gz

    mkdir ref/
    tar xzvf ${cr_ref} -C ref/ --strip-components=1
    ls -l ref/

    cellranger count \
    --id=${specimen} \
    --disable-ui \
    --localcores ${task.cpus} \
    --localmem ${task.memory.toGiga()} \
    --fastqs=reads/ \
    --transcriptome=ref/

    tar czf ${specimen}_result.tgz ${specimen}/outs/* 
    rm -r ${specimen}
    rm -r ref/*
    """
}

process CombineCRout {
    container "${container__cellranger}"
    label 'io_limited'
    errorStrategy 'finish'
    publishDir "${params.output}", mode: 'copy'

    input:
    path(out_tgz)

    output:
    path('cellranger/*')

    """
    mkdir cellranger
    for f in *.tgz; do tar -xzvf  "\$f" -C cellranger/ ; done
    """
}


def helpMessage() {
    log.info"""
    Workflow to process scRNA reads
    Usage:
    --manifest      Path to a read-pairs to be analyzed (REQUIRED)
                    'specimen'                  This can be any string, which is a sequence of alpha-numeric characters, underscores, 
                                                or dashes and no spaces, that is less than 64 characters.
                                                MUST be unique for this data set. 
                    'R1'                        First read, in fasta format and gzipped
                    'R2'                        Second / reverse read, in fasta format and gzipped  
    --cr_ref        Cellranger reference in tar.gz format
    --output        Path where to place the output files
    
    Parameters:
    """.stripIndent()
}

def ReadManifest(manifest_file){
    manifest_file.splitCsv(
        header: true, 
        sep: ","
    ).branch{
        //valid_paired_indexed:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty()) && (it.R2 != null ) && (it.R2 != "") && (!file(it.R2).isEmpty()) && (it.I1 != null ) && (it.I1 != "" ) && (!file(it.I1).isEmpty()) && (it.I2 != null ) && (it.I2 != "") && (!file(it.I2).isEmpty())
        valid_paired:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty()) && (it.R2 != null ) && (it.R2 != "") && (!file(it.R2).isEmpty())
        //valid_unpaired:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty())
        other: true
    }
}

workflow {
    if ((!params.manifest | params.help | !params.cr_ref)) {
        helpMessage()
        exit 0
    }

    // Load manifest!
    manifest = ReadManifest(
        Channel.from(
            file(params.manifest)
        )
    )

    SingleCellRNA(
        manifest.valid_paired,
        file(params.cr_ref)
    )

}
