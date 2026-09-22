# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(tidyr)

# Read the CSV files
t1_haplotypes <- read_csv("results/T1_S1/T1_S1_haplotypes.csv")
t2_haplotypes <- read_csv("results/T2_S2/T2_S2_haplotypes.csv")
t3_haplotypes <- read_csv("/Users/semiquant/Bioinformatics/Projects/BDQarmy/cloneArmyInterrogate/files/results/aln.sorted_haplotypes.csv")

# Combine the dataframes
combined_haplotypes <- bind_rows(
  t1_haplotypes %>% mutate(sample = "T1"),
  t2_haplotypes %>% mutate(sample = "T2"),
  t3_haplotypes %>% mutate(sample = "T3")
)



# combined_haplotypes = combined_haplotypes %>%
#   filter(sample == "T3")



combined_haplotypes %>%
  filter(snp_count != 1) %>%
  filter(mutations != 1) %>%
  group_by(haplotype, sample) %>%  # Group by sample to count unique occurrences
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count)) %>%
  filter(total_count > 100) %>%
  nrow()

combined_haplotypes %>%
  filter(indel_count == 1) %>%
  group_by(haplotype, sample) %>%  # Group by sample to count unique occurrences
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count)) %>%
  filter(total_count > 100) %>%
  nrow()






# Filter for snp_count == 1 and summarize by haplotype
result <- combined_haplotypes %>%
  filter(snp_count == 1) %>%
  filter(mutations == 1) %>%
  group_by(haplotype, sample) %>%  # Group by sample to count unique occurrences
  summarise(total_count = sum(count), .groups = 'drop') %>%
  arrange(desc(total_count))

result %>%
  filter(total_count > 100) %>%
  nrow()/1495*100


result %>%
  filter(total_count > 1000) %>%
  nrow()



# Count lowercase at each position
mutation_counts <- result %>%
  filter(total_count > 100) %>%
  mutate(
    chars = strsplit(haplotype, "")
  ) %>%
  unnest(chars) %>%
  group_by(haplotype) %>%
  mutate(position = row_number()) %>%
  ungroup() %>%
  filter(str_detect(chars, "[a-z]")) %>%  # Keep only lowercase positions
  group_by(position) %>%
  summarise(
    lowercase_count = n_distinct(sample)  # Count unique samples at each position
  ) %>%
  filter(lowercase_count > 0)

# Print the counts to check
print("\nMutation counts (should be max 3):")
print(mutation_counts %>% arrange(position))

# Create the plot
ggplot(mutation_counts, aes(x = position, y = lowercase_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(
    x = "Position",
    y = "Number of samples with mutation"
  ) +
  scale_y_continuous(breaks = 1:3)  # Force y-axis to show 1,2,3


