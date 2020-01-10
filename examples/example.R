library(outsider)
# import
pasta <- module_import('pasta', repo = 'dombennett/om..pasta')
# help
pasta('-h')


# Real Example
# download
wd <- file.path(tempdir(), 'example_job')
if (!dir.exists(wd)) {
  dir.create(wd)
}
# example DNA, based on MAFFT example
seq_file <- file.path(wd, 'example_seq.fasta')
url <- 'https://raw.githubusercontent.com/DomBennett/om..pasta/master/example_seq.fasta'
download.file(url = url, destfile = seq_file)
# run DNA alignment
pasta(arglist = c('-i', seq_file, '-d', 'dna',
                  '--alignment-suffix=alignment.fasta',
                  '--job=example'), outdir = wd)
