suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(dplyr)
  library(stringr)
})

input_file  <- "tb_immune_genesets_symbols.txt"
output_file <- "tb_immune_genesets.txt"

lines <- readLines(input_file)

convert_line <- function(line) {
  if (str_starts(line, "#") || nchar(trimws(line)) == 0) {
    return(NA)
  }

  tokens <- str_split(line, "\\s+")[[1]]
  geneset <- tokens[1]
  symbols <- tokens[-1]

  entrez <- mapIds(
    org.Hs.eg.db,
    keys = symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  entrez <- entrez[!is.na(entrez)]
  if (length(entrez) == 0) return(NA)

  paste(c(geneset, unique(entrez)), collapse = " ")
}

converted <- sapply(lines, convert_line)
converted <- converted[!is.na(converted)]

writeLines(converted, output_file)

cat("✓ Converted gene sets written to:", output_file, "\n")
