params.input = "Aurelia_Transcriptomes.input.csv"
params.output = "Aurelia_Transcriptomes"

input_ch = Channel.fromPath(params.input).splitCsv(header: true).map { row -> file(row.path) }.into { input_ch_1; input_ch_2; input_ch_3 }

process blast {
  input:
    set val(a), val(b), file("left.fa"), file("right.fa") from input_ch_1.combine(input_ch_2).map { tuple(it[0].getBaseName(0), it[1].getBaseName(0), it[0], it[1]) }
  output:
    file("${a}_vs_${b}.txt.gz") into blast_ch

  module 'DIAMOND/2.0.4-GCC-9.3.0'

  """
  diamond blastp -p ${task.cpus} -q left.fa -d right.fa | gzip -c > ${a}_vs_${b}.txt.gz
  """
}

process orthofinder {
  input:
    file(fasta) from input_ch_3.collect()
    file(blast) from blast_ch.collect()

  output:
    file(params.output)

  module 'Biopython/1.75-foss-2020a-Python-3.8.2'
  module 'OrthoFinder/2.5.2'
  publishDir ".", overwrite: true, mode: 'link'

  """
set -o errexit
set -o pipefail

cat << EOF > fasta.txt
${fasta.join('\n')}
EOF

cat << EOF > blast.txt
${blast.join('\n')}
EOF

prepare_orthofinder_input.py fasta.txt blast.txt

orthofinder -a ${task.cpus} -n Nextflow -b .

mv OrthoFinder/Results_Nextflow ${params.output}
  """
}
