# mass-accuracy-error-tracking

Quick check of the mass accuracy error using FragPipe output

The libraries required are the following:
```r
library(tidyverse) # all the data wrangling
library(here) # we need this library to make our lives easier to deal with the path to the files
library(janitor) # we need the clean_names() function
library(ggtext) # to make the customized annotation
```

All we need to do is load the `psm.tsv` files of interest, as in the following example.
We added a new column containing the mass accuracy error in ppm using the ratio of ∆m/z (observed m/z minus theoretical m/z) over theoretical m/z multiplied by 1e6 (making ppm). We are taking advantage of the fact that FragPipe gives us the **observed_m_z** and **calculated_m_z** columns.

```r
HeLa_DDA <- read_tsv("QC/hela01/psm.tsv") %>% 
  clean_names() %>% 
  dplyr::mutate(sample = "HeLa 01") %>%
  rbind(
    read_tsv("QC/hela02/psm.tsv") %>%
      clean_names() %>%
      dplyr::mutate(sample = "HeLa 02")
    ) %>% 
      rbind(
        read_tsv("QC/hela03/psm.tsv") %>%
          clean_names() %>%
          dplyr::mutate(sample = "HeLa 03")
        ) %>% 
          rbind(
            read_tsv("QC/hela04/psm.tsv") %>%
              clean_names() %>%
              dplyr::mutate(sample = "HeLa 04")
            ) %>% 
  dplyr::mutate(delta_mass_ppm = (observed_m_z-calculated_m_z)/calculated_m_z*1e6) # this is where the magic happens
```

We can add usefull annotations, like the proportion of ions out of the calibrated mass. Let's put the treshhold at 10 ppm.

```r
error_dataframe <- HeLa_DDA %>% 
  dplyr::filter(abs(delta_mass_ppm) > 10) %>% 
  group_by(sample) %>%
  summarise(n = n()) %>%
  left_join(HeLa_DDA %>% 
              group_by(sample) %>%
              summarise(total = n())
            ) %>% 
  dplyr::mutate(percent = round((n/total)*100, 2)
                )
```

Finally, we can plot the precision error over retention time for each sample. Note that we are filtering out ions with error >100 ppm. This will make the difference visible at scale. Our plot contains valuable piece of information like:
- The proportion of ions out of calibration.
- A [Generalized Additive Model](https://medium.com/data-science-in-your-pocket/understanding-generalized-additive-models-gams-for-regression-with-mathematics-a19bb9915a73) (GAM) fitted to the data (should be close to zero).
- The general profile of error over the elution time for each sample.

```r
error_plot <- HeLa_DDA %>% 
  dplyr::filter(abs(delta_mass_ppm) < 100) %>% 
  ggplot(aes(x = retention/60,
             y = delta_mass_ppm)
         ) +
  geom_point(alpha = 0.1, color = "black", size = 1) +
  geom_smooth(method = "gam", se = FALSE, color = "blue", linewidth = 0.5) +
  facet_wrap(~sample, ncol = 2, scales = "free_y") +
  geom_text(data = error_dataframe,
            aes(label = paste0(percent, "% ≥ 10 ppm"), x = 22, y = -30),
            size = 5, color = "red") +
  geom_hline(yintercept = c(10, 0, -10), color = "red", linetype = "dashed", linewidth = 0.2) +
  labs(title = "Mass error in ppm by sample",
       caption = "Data from 4 different HeLa injections analyzed by LC-MS/MS in DDA mode on Orbitrap Exploris 480.<br>
       <span style='color:blue;'>**Blue line**</span> - Generalized Additive Model smooth line.<br>
       <span style='color:red;'>**Red dashed lines**</span> - 10 ppm threshold.<br>
       <span style='color:red;'>**Red text**</span> - percentage of ions with mass error ≥ 10 ppm.",
      x = "Retention time (min)",
       y = "Mass error (ppm)") +
  theme_minimal() +
  theme(text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_markdown(),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(face = "bold"),
        axis.ticks = element_line(color = "black"),
        line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "none"
        )

# save the plot in your working directory
ggsave(
  "QC_error_tracking.png",
  plot = error_plot,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)
```

<p align="center">
<img src="https://github.com/41ison/mass-accuracy-error-tracking/blob/main/HeLa_QC.png" width="500">
</p>

