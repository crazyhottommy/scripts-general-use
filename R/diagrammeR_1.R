library(DiagrammeR)

grViz("
digraph boxes_andcircles {

  node [shape = box]
  ChIPseq_fastq; bam; rm_dup_sorted_bam;bigwig; visualization;peaks; downstream_analysis

  edge [arrowhead = diamond]
  ChIPseq_fastq-> bam [label = ' fastqc bowtie',
              fontname = Helvetica];
  bam-> rm_dup_sorted_bam[label = ' samtools',
                          fontname = Helvetica]
  rm_dup_sorted_bam-> bigwig[label = ' deeptools'
                            fontname = Helvetica]
  bam-> peaks [label = 'MACS2',
              fontname = Helvetica]
  peaks-> downstream_analysis
  rm_dup_sorted_bam-> downstream_analysis
  bigwig-> visualization[label =' IGV',
              fontname = Helvetica]

  node [shape = box,
        color = blue]
  WGS_bam; remapped_bam; discordants_bam; splitter_bam; SVs; dysfunctional_enhancer; annotation; prioritization
  filtered_SVs
  edge [arrowhead = diamond,
        color =blue]
  WGS_bam-> remapped_bam[label = ' speedseq realign',
                        fontname = Helvetica]
  WGS_bam-> discordants_bam
  WGS_bam-> splitter_bam
  remapped_bam-> SVs
  discordants_bam-> SVs[label =' speedseq sv',
                        fontname = Helvetica]
  splitter_bam-> SVs
  peaks-> dysfunctional_enhancer[label = ' bedtools',
                              fontname = Helvetica,
                              color = red]
  filtered_SVs-> dysfunctional_enhancer[color = red]
  SVs-> filtered_SVs[label = ' vawk pybamview',
                              fontname = Helvetica]
  dysfunctional_enhancer-> annotation[label =' VEP',
                                      fontname = Helvetica]
  annotation-> prioritization[label = ' Gemini',
                              fontname = Helvetica]


}
")

