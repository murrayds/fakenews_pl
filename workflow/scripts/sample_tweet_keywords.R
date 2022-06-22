# Convert the nodes from a CSV format into an XYZ format for touchdesigner
library(readr)
library(dplyr)
library(vroom)
library(stringr)

tweets <- vroom(snakemake@input[[1]],
                delim = "\t",
                col_types = c(id_str = "c", userid = "c", created_at = "_", text_keyword = "c"))

nodes <- read_csv(snakemake@input[[2]])

# Random seed, ensure that we have the data to work with
set.seed(1)

# Sample people, weighing by their numfake
node_sample <- nodes %>%
  sample_n(1200, weight = numfake)

# This is where the logic is actually get the sample and format the sentence
df <- tweets %>%
  inner_join(node_sample, by = c('userid')) %>%
  group_by(userid) %>%
  filter(row_number() < 5) %>%
  summarize(
    # Remove the pipe characters and marge keywords across group
    label = gsub("|", "", paste0(text_keyword, collapse = " "), fixed = T),
    # Then, select unique terms, then take the first 6 words of the aggregate list
    label = word(paste0(unique(unlist(str_split(label, " "))), collapse = " "), 1, 6),
    # get other data for each group
    x = first(x),
    y = first(y),
    z = first(z),
    numfake = first(numfake)
  ) %>%
  rowwise() %>%
  mutate(
    dist = sqrt(x^2 + y^2 + z^2),
    dirvec_x = x / dist,
    dirvec_y = y / dist,
    dirvec_z = z / dist,
  ) %>%
  select(x, y, z, dist, dirvec_x, dirvec_y, dirvec_z, label, numfake) %>%
  filter(!is.na(label))

# Write the output
write_csv(df, snakemake@output[[1]])
