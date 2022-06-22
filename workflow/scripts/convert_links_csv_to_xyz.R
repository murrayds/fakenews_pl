# Convert the nodes from a CSV format into an XYZ format for touchdesigner
library(readr)
library(dplyr)

df <- read_csv(snakemake@input[[1]]) %>%
  select(x, y, z, dirvec_x, dirvec_y, dirvec_z, dist)

write_delim(df,
            snakemake@output[[1]],
            delim = "\t",
            col_names = FALSE)
