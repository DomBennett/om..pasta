#' @name pasta
#' @title Practical Alignment using SATe and Transitivity
#' @description Run pasta in R.
#' @param arglist Arguments for pasta
#' @param outdir Filepath to output files are placed
#' @example /examples/example.R
#' @export
pasta <- function(arglist = arglist_get(...), outdir = getwd()) {
  files_to_send <- filestosend_get(arglist = arglist)
  parsed_arglist <- arglist_parse(arglist = arglist, normalise_paths = TRUE)
  otsdr <- outsider_init(pkgnm = 'om..pasta', cmd = 'run_pasta.py',
                         arglist = parsed_arglist, wd = outdir,
                         files_to_send = files_to_send)
  run(otsdr)
}
