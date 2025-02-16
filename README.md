Consensus maker

This tool makes consensuses from an input FASTA file and a VCF file containing variant calls.
Currently does not support  symbolic structural variant alleles other than <DEL>.

usage: consensus_maker.py [-h] [-c CONS_COUNT] [-l CONS_LENGTH] [-m MIN_FREQ] [-s] [-e CONS_EXTRA_BASES] [--version] input_fasta vcf output_fasta

This tool makes consensuses from an input FASTA file and a VCF file containing variant calls.

positional arguments:
  input_fasta           Input FASTA file path.
  vcf                   Input VCF file path.
  output_fasta          Output FASTA path.

optional arguments:
  -h, --help            show this help message and exit
  -c CONS_COUNT, --cons_count CONS_COUNT
                        Number of consensuses (default: 1000)
  -l CONS_LENGTH, --cons_length CONS_LENGTH
                        Length of consensuses (default: 160)
  -m MIN_FREQ, --min_freq MIN_FREQ
                        Minimum allele frequency to include a variant in the consensus (default: 0.0)
  -s, --sorted          Sort consensuses by coordinate (default: False)
  -e CONS_EXTRA_BASES, --cons_extra_bases CONS_EXTRA_BASES
                        (Advanced) Extra bases (max) to be generated in each consensus to compensate for deletions (default: 200)
  --version             show program's version number and exit
