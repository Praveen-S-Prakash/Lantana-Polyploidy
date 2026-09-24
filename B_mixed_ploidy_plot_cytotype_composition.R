
## Plot cytotype composition across locations
## Includes relative stacked bar plot and overall pie chart

library(tidyverse)

setwd("/Users/praveenp/Desktop/Chapter_2_analysis/27_10_25/")

# Load ploidy data
df <- read.csv("Ploidy_table.csv")

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = -Location,
    names_to = "Ploidy",
    values_to = "Count"
  )

# Define colours for ploidies
ploidy_colors <- c(
  "Diploids" = "#8863ce",
  "Triploids" = "#edf41b",
  "Tetraploids" = "#33d5db",
  "Hexaploids" = "#ec501b"
)

# Stacked bar plot showing relative proportions
p_bar <- ggplot(
  df_long,
  aes(x = Location, y = Count, fill = Ploidy)
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  scale_fill_manual(
    values = ploidy_colors
  ) +
  labs(
    x = "Location",
    y = "Proportion of individuals",
    title = "Cytotype composition across locations"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  )

print(p_bar)

# Save stacked bar plot
ggsave(
  "Cytotype_composition_barplot.pdf",
  plot = p_bar,
  width = 8,
  height = 5
)

# Sum counts across all locations for each ploidy
total_counts <- df_long %>%
  group_by(Ploidy) %>%
  summarise(
    Total = sum(Count, na.rm = TRUE)
  )

# Overall cytotype composition pie chart
p_pie <- ggplot(
  total_counts,
  aes(x = "", y = Total, fill = Ploidy)
) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = NA
  ) +
  coord_polar("y") +
  scale_fill_manual(
    values = ploidy_colors
  ) +
  theme_void(base_size = 14) +
  labs(
    title = "Overall cytotype composition"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

print(p_pie)

# Save pie chart
ggsave(
  "Cytotype_composition_pie_total.pdf",
  plot = p_pie,
  width = 6,
  height = 6
)
