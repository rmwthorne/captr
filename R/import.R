#' Converts comma-containing strings to numeric values
#'
#' @param x A character object
#' @return The numeric value of \code{x}.
#' @examples
#' strip_comma(c("1,000", "2,100"))
#' @export

strip_comma <- function (x) {
  as.numeric(gsub(",","",x))
}


#' Removes _# if at the end of a string
#'
#' @param x A character object
#' @return The cleaned version of \code{x}.
#' @examples
#' strip_suffix("MYC_1")
#' @export

strip_suffix <- function(x) gsub(pattern = "_1|_2|_3", replacement = "", x = x)


#' Splits a UCSC position into three valid values
#'
#' @param x A character object in the format chr:start-stop
#' @return A 3-column data frame form the values in \code{x}.
#' @examples
#' split_position
#'
split_position <- function(chr) {
  # Regex-split a Chr:Start-Stop string
  if (!grepl(":", chr) | !grepl("-", chr)) {
    # Error handling: make sure chr is formatted correctly
    stop("chr must be formatted as: \"chr:start-stop\"", call. = F)
  }
  # Return a list of chromosome, start and stop positions
  data.frame("to.split" = chr) %>%
    separate(to.split, c("split.chr","split.start","split.stop"), sep = ":|-") %>%
    mutate(split.start = strip_comma(split.start),
           split.stop = strip_comma(split.stop))
}



tidy_gfc <- function(gene.expn, window, smooth = 10, path = "data") {
  message("Extracting data for ", gene.expn, "...")

  # Load and merge the dataframes
  gene.files <- list.files(path = path, pattern = gene.expn, full.names = T)
  # message("Files matching pattern: \n - ", paste0(gene.files, collapse = "\n - "))
  r_an.df <- read.table(gene.files[grepl(pattern = "R_analysis", x = gene.files)])
  r_an.df <- r_an.df %>% mutate("Fragment" = rownames(r_an.df))
  deseq.df <- read.table(gene.files[grepl(pattern = "DESeq_analysis", x = gene.files)])
  deseq.df <- deseq.df %>% mutate("Fragment" = rownames(deseq.df))
  df <- full_join(x = r_an.df, y = deseq.df, by = c("Fragment"))
  df <- df %>% separate(Fragment, c("Chr","Start","Stop"), sep=":|-") %>%
    mutate(Start = strip_comma(Start), Stop = strip_comma(Stop)) %>%
    mutate(Position = (Start+Stop)/2 )

  # Filter for the datapoints in the window
  window <- split_position(window) %>% as.list()
  df <- df %>% filter_(~Chr == window[[1]],
                       ~Position > window[[2]],
                       ~Position < window[[3]])

  df <- df %>%
    mutate_(Locus = ~gene.expn) %>%
    mutate_(S.condA_mean = ~rollmean(condA_mean, smooth, fill = "extend")) %>%
    mutate_(S.condB_mean = ~rollmean(condB_mean, smooth, fill = "extend")) %>%
    mutate_(S.diff_mean = ~rollmean(dif_mean, smooth, fill = "extend")) %>%
    mutate_(S.condA_SD = ~rollmean(condA_SD, smooth, fill = "extend")) %>%
    mutate_(S.condB_SD = ~rollmean(condB_SD, smooth, fill = "extend"))

  # Tidy
  df.tidy <- df %>% gather_(key_col = "Condition", value_col = "Normalised.Mean",
                            gather_cols = c("condA_mean", "condB_mean"))
  df.tidy <- df %>% gather_(key_col = "S.Condition", value_col = "S.Normalised.Mean",
                            gather_cols = c("S.condA_mean", "S.condB_mean"))
  df.tidy
}

