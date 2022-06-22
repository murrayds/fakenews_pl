# Convert the nodes from a CSV format into an XYZ format for touchdesigner
library(readr)
library(dplyr)

df <- read_csv(snakemake@input[[1]])


df.xyz <- df %>%
  mutate(
    A = numfake / max(numfake)
  ) %>%
  select(x, y, z, A) %>%
  arrange(x) %>%
  slice(1:17280)


write_delim(df.xyz %>% select(x, y, z, A),
            snakemake@output[[1]],
            delim = "\t",
            col_names = FALSE)
