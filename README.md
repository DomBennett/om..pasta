
<!--
The README should be used to describe the program. It acts like the homepage of
your module.

Edit README.Rmd not README.md. The .Rmd file can be knitted to parse real-code
examples and show their output in the .md file.

To knit, use devtools::build_readme() or outsider.devtools::build()

Edit the template to describe your program: how to install, import and run;
run exemplary, small demonstrations; present key arguments; provide links and
references to the program that the module wraps.

Learn more about markdown and Rmarkdown:
https://daringfireball.net/projects/markdown/syntax
https://rmarkdown.rstudio.com/
-->

# Run [`pasta`](https://github.com/smirarab/pasta) with `outsider` in R

[![Build
Status](https://travis-ci.org/dombennett/om..pasta.svg?branch=master)](https://travis-ci.org/dombennett/om..pasta)

> pasta (Practical Alignment using SATe and Transitivity): Multiple
> sequence alignment of biological sequences.

<!-- Install information -->

## Install and look up help

``` r
library(outsider)
#> ----------------
#> outsider v 0.1.0
#> ----------------
#> - Security notice: be sure of which modules you install
module_install(repo = "dombennett/om..pasta")
#> -----------------------------------------------------
#> Warning: You are about to install an outsider module!
#> -----------------------------------------------------
#> Outsider modules install and run external programs
#> via Docker <https://www.docker.com>. These external
#> programs may communicate with the internet and could
#> potentially be malicious.
#> 
#> Be sure to know the module you are about to install:
#> Is it from a trusted developer? Are colleagues using
#> it? Is it supposed to download lots of data? Is it
#> well used (e.g. check number of stars on GitHub)?
#> -----------------------------------------------------
#>  Module information
#> -----------------------------------------------------
#> program: pasta
#> details: PASTA (Practical Alignment using Sate and TrAnsitivity) algorithm.
#> docker: dombennett
#> github: dombennett
#> url: https://github.com/dombennett/om..pasta
#> image: dombennett/om_pasta
#> container: om_pasta
#> package: om..pasta
#> Travis CI: Passing
#> -----------------------------------------------------
#> Enter any key to continue or press Esc to quit
#module_help(repo = "dombennett/om..pasta")
```

<!-- Detailed examples -->

## Aligning a small set of DNA sequences

> The below example is adapted from the [PASTA
> tutorial](https://github.com/smirarab/pasta/blob/master/pasta-doc/pasta-tutorial.md#using-pasta).

``` r
library(outsider)
# import
pasta <- module_import('pasta', repo = 'dombennett/om..pasta')
# get help
pasta('-h')
#> Usage: run_pasta.py [options] <settings_file1> <settings_file2> ...
#> 
#> 
#> PASTA performs iterative realignment and tree inference, similar to SATe, but
#> uses a very different merge algorithm which improves running time, memory
#> usage, and accuracy. The current code is heavily based on SATe, with lots of
#> modifications, many related to algorithmic differences between PASTA and SATe,
#> but also many scalability improvements (parallelization, tree parsing,
#> defaults, etc.)
#> 
#> Minimally you must provide a sequence file (with the '--input' option); a
#> starting tree is optional. By default, important algorithmic parameters are
#> set based on automatic rules.
#> 
#> The command line allows you to alter the behavior of the algorithm
#> (termination criteria, when the algorithm switches to "Blind" acceptance of
#> new alignments, how the tree is decomposed to find subproblems to be used, and
#> the external tools to use).
#> 
#> Options can also be passed in as configuration files.
#> 
#> With the format:
#> ####################################################
#> [commandline]
#> option-name = value
#> 
#> [sate]
#> option-name = value
#> ####################################################
#> 
#> With every run, PASTA saves the configuration file for that run as a temporary
#> file called [jobname]_temp_pasta_config.txt in your output directory.
#> 
#> Configuration files are read in the order they occur as arguments (with values
#> in later files replacing previously read values). Options specified in the
#> command line are read last. Thus these values "overwrite" any settings from
#> the configuration files. Note that the use of --auto option can overwrite some
#> of the other options provided by commandline or through configuration files.
#> 
#> 
#> Options:
#>   --version             show program's version number and exit
#>   -h, --help            show this help message and exit
#> 
#>   commandline options:
#>     -a, --aligned       If used, then the input file be will treated as
#>                         aligned for the purposes of the first round of tree
#>                         inference (the algorithm will start with tree
#>                         searching on the input before re-aligning). This
#>                         option only applies if a starting tree is NOT given.
#>     --alignment-suffix=ALIGNMENT-SUFFIX
#>                         suffix for alignment name (default: .marker001.[input
#>                         name].aln)
#>     --auto              This option is mostly for backward compatibility. If
#>                         used, then automatically identified default values for
#>                         the max_subproblem_size, number of cpus, tools,
#>                         breaking strategy, masking criteria, and stopping
#>                         criteria will be used. This is just like using the
#>                         default options. However, [WARNING] when auto option
#>                         is used PASTA overrides the value of these options
#>                         even if you have supplied them; we recommend that you
#>                         run this option with --exportconfig to see the exact
#>                         set of options that will be used in your analysis.
#>     -d DATATYPE, --datatype=DATATYPE
#>                         Specify DNA, RNA, or Protein to indicate what type of
#>                         data is specified. Note that this option is NOT
#>                         automatically determined [default: dna]
#>     --exportconfig=EXPORTCONFIG
#>                         Export the configuration to the specified file and
#>                         exit. This is useful if you want to combine several
#>                         configurations and command line settings into a single
#>                         configuration file to be used in other analyses.
#>     -i INPUT, --input=INPUT
#>                         input sequence file
#>     -j JOB, --job=JOB   job name [pastajob]
#>     --keepalignmenttemps
#>                         Keep even the realignment temporary running files
#>                         (this only has an effect if keeptemp is also
#>                         selected).
#>     -k, --keeptemp      Keep temporary running files? [default: disabled]
#>     --missing=MISSING   How to deal with missing data symbols. Specify either
#>                         "Ambiguous" or "Absent" if the input data contains
#>                         ?-symbols
#>     -m, --multilocus    Analyze multi-locus data? NOT SUPPORTED IN CURRENT
#>                         PASTA version.
#>     --raxml-search-after
#>                         If used, the completion of the PASTA algorithm will be
#>                         followed by a tree search using RAxML on the masked
#>                         alignment. This can be useful if a very fast and
#>                         approximate tree estimator is used during the PASTA
#>                         algorithm. [default: disabled]
#>     --temporaries=TEMPORARIES
#>                         directory that will be the parent for this job's
#>                         temporary file [default in PASTA home]
#>     --timesfile=TIMESFILE
#>                         optional file that will store the times of events
#>                         during the PASTA run. If the file exists, new lines
#>                         will be
#>     -t TREEFILE, --treefile=TREEFILE
#>                         starting tree file
#>     --two-phase         If used, then the program will not perform the PASTA
#>                         algorithm. Instead it will simply call the sequence
#>                         aligner to align the entire dataset then will call the
#>                         tree estimator to obtain the tree.
#>     --untrusted         If used, then the data in the input file will be
#>                         parsed using a more careful procedure. This will
#>                         generate more helpful error messages, but will use
#>                         more memory and be much slower for large inputs. If
#>                         this option is omitted, the error messages resulting
#>                         from invalid input data will be more cryptic.
#> 
#>   SATe acceptance options:
#>     --blind-after-iter-without-imp=#
#>                         Maximum number of iterations without an improvement in
#>                         likelihood score that PASTA will run before switching
#>                         to blind mode. [default: disabled]
#>     --blind-after-time-without-imp=#.#
#>                         Maximum time (in seconds) that PASTA will run without
#>                         an improvement in likelihood score before switching to
#>                         blind mode. [default: disabled]
#>     --blind-after-total-iter=#
#>                         Maximum number of iterations that PASTA will run
#>                         before switching to blind mode. [default: 0]
#>     --blind-after-total-time=#.#
#>                         Maximum time (in seconds) that PASTA will run before
#>                         switching to blind mode. [default: disabled]
#>     --no-blind-mode-is-final
#>                         When the blind mode is final, then PASTA will never
#>                         leave blind mode once it is has entered blind mode.
#>     --move-to-blind-on-worse-score
#>                         If True then PASTA will move to the blind mode as soon
#>                         it encounters a tree/alignment pair with a worse
#>                         score. This is essentially the same as running in
#>                         blind mode from the beginning, but it does allow one
#>                         to terminate a run at an interval from the first time
#>                         the algorithm fails to improve the score.
#> 
#>   SATe decomposition options:
#>     --break-strategy=BREAK_STRATEGY
#>                         The method for choosing an edge when bisecting the
#>                         tree during decomposition [default: mincluster]
#>     --max-subproblem-frac=#.#
#>                         The maximum size (number of leaves) of subproblems
#>                         specified in terms as a proportion of the total number
#>                         of leaves.  When a subproblem contains this number of
#>                         leaves (or fewer), then it will not be decomposed
#>                         further. [default: automatically picked based on
#>                         alignment size]
#>     --max-subproblem-size=#
#>                         The maximum size (number of leaves) of subproblems.
#>                         When a subproblem contains this number of leaves (or
#>                         fewer), then it will not be decomposed further.
#>                         [default: automatically picked based on alignment
#>                         size]
#>     --max-subtree-diameter=#.#
#>                         The maximum diameter of each subtree. [default: 2.5]
#>     --min-subproblem-size=#
#>                         The minimum size (number of leaves) of subproblems.
#>                         [default: 0]
#> 
#>   SATe filtering options:
#>     --treeshrink-filter
#>                         If used, then the inferred FastTree will be filtered
#>                         by TreeShrink for long branch outliers.
#> 
#>   SATe output options:
#>     -o OUTPUT_DIRECTORY, --output-directory=OUTPUT_DIRECTORY
#>                         directory for output files (defaults to input file
#>                         directory)
#>     --no-return-final-tree-and-alignment
#>                         Return the best likelihood tree and alignment pair
#>                         instead of those from the last iteration; this is
#>                         discouraged with masking option enabled.
#> 
#>   SATe platform options:
#>     --max-mem-mb=#      The maximum memory available to OPAL (for the Java
#>                         heap size when running Java tools).
#>     --num-cpus=#        The number of processing cores that you would like to
#>                         assign to PASTA.  This number should not exceed the
#>                         number of cores on your machine. [default: number of
#>                         cores available on the machine]
#> 
#>   SATe searching options:
#>     --mask-gappy-sites=#
#>                         The minimum number of non-gap characters required in
#>                         each column passed to the tree estimation step.
#>                         Columns with fewer non-gap characters than the given
#>                         threshold will be masked out before passing the
#>                         alignment into the tree estimation module. These
#>                         columns will be present in the final alignment.
#>                         [default: 0.1% of alignment size]
#>     --start-tree-search-from-current
#>                         If selected that the tree from the previous iteration
#>                         will be given to the tree searching tool as a starting
#>                         tree.
#> 
#>   SATe spaning-tree options:
#>     --build-MST         Construct the spanning tree using minimum spanning
#>                         tree algorithm [default: False]
#> 
#>   SATe termination options:
#>     --after-blind-iter-term-limit=#
#>                         The maximum number of iteration that the PASTA
#>                         algorithm will run after PASTA has entered blind mode.
#>                         If the number is less than 1, then no iteration limit
#>                         will be used. [default: disabled]
#>     --after-blind-iter-without-imp-limit=#
#>                         The maximum number of iterations without an
#>                         improvement in score that the PASTA algorithm will run
#>                         after entering BLIND mode.  If the number is less than
#>                         1, then no iteration limit will be used. [default:
#>                         disabled]
#>     --after-blind-time-term-limit=#.#
#>                         Maximum time (in seconds) that PASTA will continue
#>                         starting new iterations of realigning and tree
#>                         searching after PASTA has entered blind mode. If the
#>                         number is less than 0, then no time limit will be
#>                         used. [default: disabled]
#>     --after-blind-time-without-imp-limit=#.#
#>                         Maximum time (in seconds) since the last improvement
#>                         in score that PASTA will continue starting new
#>                         iterations of realigning and tree searching after
#>                         entering BLIND mode. If the number is less than 0,
#>                         then no time limit will be used. [default: disabled]
#>     --iter-limit=#      The maximum number of iteration that the PASTA
#>                         algorithm will run.  If the number is less than 1,
#>                         then no iteration limit will be used. [default: 3]
#>     --iter-without-imp-limit=#
#>                         The maximum number of iterations without an
#>                         improvement in score that the PASTA algorithm will
#>                         run.  If the number is less than 1, then no iteration
#>                         limit will be used. [default: disabled]
#>     --time-limit=#.#    Maximum time (in seconds) that PASTA will continue
#>                         starting new iterations of realigning and tree
#>                         searching. If the number is less than 0, then no time
#>                         limit will be used. [default: disabled]
#>     --time-without-imp-limit=#.#
#>                         Maximum time (in seconds) since the last improvement
#>                         in score that PASTA will continue starting new
#>                         iterations of realigning and tree searching. If the
#>                         number is less than 0, then no time limit will be
#>                         used. [default: disabled]
#> 
#>   SATe tools options:
#>     --aligner=ALIGNER   The name of the alignment program to use for
#>                         subproblems. [default: mafft]
#>     --merger=MERGER     The name of the alignment program to use to merge
#>                         subproblems. [default: OPAL]
#>     --tree-estimator=TREE_ESTIMATOR
#>                         The name of the tree inference program to use to find
#>                         trees on fixed alignments. [default: fasttree]

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
#> PASTA INFO: Reading input sequences from 'example_seq.fasta'...
#> PASTA INFO: Configuration written to "/working_dir/example_temp_pasta_config.txt".
#> 
#> PASTA INFO: Directory for temporary files created at /root/.pasta/example/tempab7o4zwo
#> PASTA INFO: Name translation information saved to /working_dir/example_temp_name_translation.txt as safe name, original name, blank line format.
#> PASTA INFO: Creating a starting tree for the PASTA algorithm...
#> PASTA INFO: Performing initial alignment of the entire data matrix...
#> PASTA INFO: Performing initial tree search to get starting tree...
#> PASTA INFO: Starting PASTA algorithm on initial tree...
#> PASTA INFO: Max subproblem set to 5
#> PASTA INFO: Step 0. Realigning with decomposition strategy set to mincluster
#> PASTA INFO: Step 0. Alignment obtained. Tree inference beginning...
#> PASTA INFO: realignment accepted and score improved.
#> PASTA INFO: current score: -21536.292, best score: -21536.292
#> PASTA INFO: TreeShrink option has been turned off!
#> PASTA INFO: Step 1. Realigning with decomposition strategy set to mincluster
#> PASTA INFO: Step 1. Alignment obtained. Tree inference beginning...
#> PASTA INFO: realignment accepted despite the score not improving.
#> PASTA INFO: current score: -21478.217, best score: -21478.217
#> PASTA INFO: TreeShrink option has been turned off!
#> PASTA INFO: Step 2. Realigning with decomposition strategy set to mincluster
#> PASTA INFO: Step 2. Alignment obtained. Tree inference beginning...
#> PASTA INFO: realignment accepted despite the score not improving.
#> PASTA INFO: current score: -21478.217, best score: -21478.217
#> PASTA INFO: TreeShrink option has been turned off!
#> PASTA INFO: Writing resulting alignment to /working_dir/example.alignment.fasta
#> PASTA INFO: Writing resulting tree to /working_dir/example.tre
#> PASTA INFO: Writing resulting likelihood score to /working_dir/example.score.txt
#> PASTA INFO: The resulting alignment (with the names in a "safe" form) was first written as the file "/working_dir/example_temp_iteration_2_seq_alignment.txt"
#> PASTA INFO: The resulting tree (with the names in a "safe" form) was first written as the file "/working_dir/example_temp_iteration_2_tree.tre"
#> PASTA INFO: Total time spent: 20.288654565811157s

# Verify results
list.files(wd)
#>  [1] "example_seq.fasta"                                             
#>  [2] "example_temp_iteration_0_seq_alignment.txt"                    
#>  [3] "example_temp_iteration_0_seq_unmasked_alignment.gz"            
#>  [4] "example_temp_iteration_0_tree.tre"                             
#>  [5] "example_temp_iteration_1_seq_alignment.txt"                    
#>  [6] "example_temp_iteration_1_seq_unmasked_alignment.gz"            
#>  [7] "example_temp_iteration_1_tree.tre"                             
#>  [8] "example_temp_iteration_2_seq_alignment.txt"                    
#>  [9] "example_temp_iteration_2_seq_unmasked_alignment.gz"            
#> [10] "example_temp_iteration_2_tree.tre"                             
#> [11] "example_temp_iteration_initialsearch_seq_alignment.txt"        
#> [12] "example_temp_iteration_initialsearch_seq_unmasked_alignment.gz"
#> [13] "example_temp_iteration_initialsearch_tree.tre"                 
#> [14] "example_temp_name_translation.txt"                             
#> [15] "example_temp_pasta_config.txt"                                 
#> [16] "example.alignment.fasta"                                       
#> [17] "example.err.txt"                                               
#> [18] "example.out.txt"                                               
#> [19] "example.score.txt"                                             
#> [20] "example.tre"
```

The final alignment is contained within the “example.alignment.fasta”
file.

<!-- Remove module after running above example -->

### Key arguments

All PASTA arguments can be provided to the `arglist` as a character
vector. Below is a table describing some of the most common. All
arguments can be displayed using the `-h`
flag.

| Argument         | Usage                    | Description                                         |
| ---------------- | ------------------------ | --------------------------------------------------- |
| i                | \-i `file`               | Specify input sequence file                         |
| t                | \-t `file`               | Specify starting tree file                          |
| d                | \-d `type`               | Specify sequence type (protein, RNA, DNA)           |
| alignment suffix | –alignment-suffix=`name` | Specify the word ending for the resulting alignment |
| job              | –job=`name`              | Specify a name for the job                          |

In addition to the arguments provided in the `arglist` a user can
provide an `outdir` specifying the file location of where all resulting
PASTA analysis files should be returned.

## Links

Find out more by visiting [PASTA’s
homepage](https://github.com/smirarab/pasta#options).

## Please cite

  - Mirarab S, Nguyen N, Warnow T. PASTA: ultra-large multiple sequence
    alignment. Sharan R, ed. Res Comput Mol Biol. 2014:177-191.
  - Mirarab S, Nguyen N, Guo S, Wang L-S, Kim J, Warnow T. PASTA:
    Ultra-Large Multiple Sequence Alignment for Nucleotide and
    Amino-Acid Sequences. J Comput Biol. 2015;22(5):377-386.
    <doi:10.1089/cmb.2014.0156>.
  - Bennett et al. (2020). outsider: Install and run programs, outside
    of R, inside of R. *Journal of Open Source Software*, In
review

## <!-- Footer -->

<img align="left" width="120" height="125" src="https://raw.githubusercontent.com/ropensci/outsider/master/logo.png">

**An `outsider` module**

Learn more at [outsider
website](https://docs.ropensci.org/outsider/). Want to build your
own module? Check out [`outsider.devtools`
website](https://docs.ropensci.org/outsider.devtools/).
