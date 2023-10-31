################################################################################
### Plotting alpha and beta estimates of TL relationship with complex traits ###
### Author: Samuel Moix                                                      ###
### Date: 17.03.2023                                                         ###
################################################################################

################################################
### Libraries ##################################
library(ggplot2)
library(dplyr)
library(readr)
require(ggforce)
require(gtable_filter)
require(ggExtra)
library(ggforestplot)

################################################
### Working directories ########################

### Input files
LM_file <- "/SET/PATH/TO/DIRECTORY/lm_res_scaled_no_outliers_no_bc.tsv"
metadata_file <- "/SET/PATH/TO/DIRECTORY/TL_metadata.csv"
MR_res_forward_file <- "/SET/PATH/TO/DIRECTORY/TELOMERE_on_TRAIT_mr_results.txt"
MR_res_reverse_file <-"/SET/PATH/TO/DIRECTORY/TRAIT_on_TELOMERE_mr_results.txt"
indep_trait <- 141

################################################
### Load data ##################################

### Load metadata
metadata <- read.csv(metadata_file, sep = ',', header = TRUE)
# Add a star to diseases
metadata <- metadata %>%
  mutate(description = ifelse(Group == "Diseases", paste0("*", description), description))


### Load linear regression file
lm_results <- read_delim(LM_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)

### Load bi-directional MR files
MR_F <- read.table(MR_res_forward_file, header = T, sep = "\t")
MR_R <- read.table(MR_res_reverse_file, header = T, sep = "\t")

################################################
### Prepare plotting dataframe #################

### Select and rename columns of interest
lm_results <- lm_results %>%
  select(c("pheno","rcoef","se","p_value"))
colnames(lm_results) <- c("pheno","b","se","pval")
lm_results$analysis <- "Linear regression"

MR_F <- MR_F %>% 
  dplyr::filter(method == "Inverse variance weighted") %>%
  select(c("outcome","b","se","pval"))
colnames(MR_F) <- c("MR_name","b","se","pval")
MR_F$analysis <- "Forward MR"

MR_R <- MR_R %>% 
  dplyr::filter(method == "Inverse variance weighted") %>%
  select(c("exposure","b","se","pval"))
colnames(MR_R) <- c("MR_name","b","se","pval")
MR_R$analysis <- "Reverse MR"

### Add description to dataframes
lm_results <- merge(lm_results, metadata[,c("pheno","description","Category")], by = "pheno")
lm_results <- lm_results[,-1] # Remove pheno columns

MR_F <- merge(MR_F, metadata[,c("MR_name","description","Category")], by = "MR_name")
MR_F <- MR_F[,-1] # Remove pheno columns

MR_R <- merge(MR_R, metadata[,c("MR_name","description","Category")], by = "MR_name")
MR_R <- MR_R[,-1] # Remove pheno columns

### Bind dataframes together
plot_df <- rbind(lm_results,MR_F,MR_R)

### Keep only trait at least significant in one of the tested methods
plot_df <-plot_df %>%
  group_by(description) %>%
  filter(any(pval < 0.05/indep_trait)) %>%
  ungroup()

### Only MR significant results
# mr_filter <- plot_df %>%
#   dplyr::filter(analysis %in% c("Forward MR", "Reverse MR")) %>%
#   group_by(description) %>%
#   dplyr::filter(any(pval < 0.05/141)) %>%
#   ungroup() %>%
#   distinct(description)
# 
# plot_df <- plot_df %>%
#   dplyr::filter(description %in% mr_filter$description)

### Reorder dataframe by lm estimate
order_df <- plot_df %>%
  dplyr::filter(analysis == "Linear regression") %>%
  arrange(b)

plot_df <- plot_df %>%
  group_by(description) %>%  
  arrange(factor(description, levels = order_df$description)) %>%
  ungroup()

################################################
### Plotting ###################################

################################################################################
################################################################################
## get strip and axis of a given panel
## Assumptions:
## - axis are adjacent to the panel, that is exactly +1/-1 positions to the t/b/l/r ...
## - ... unless there is a strip then it is +2/-2 
get_whole_panel <- function(panel_name,
                            table_layout) {
  target <- table_layout$layout %>%
    dplyr::filter(name == panel_name) %>%
    dplyr::select(row = t, col = l)
  stopifnot(NROW(target) == 1)
  pos <- unlist(target)
  dirs <- list(t = c(-1, 0),
               b = c(1, 0),
               l = c(0, -1),
               r = c(0, 1))
  filter_elems <- function(dir, 
                           type = c("axis", "strip")) {
    type <- match.arg(type)
    new_pos <- pos + dir
    res <- table_layout$layout %>%
      dplyr::filter(grepl(type, name),
                    l == new_pos["col"],
                    t == new_pos["row"]) %>%
      dplyr::pull(name)
    if (length(res)) res else NA
  }
  strip <- purrr::map_chr(dirs, filter_elems, type = "strip")
  strip <- strip[!is.na(strip)]
  dirs[[names(strip)]] <- 2 * dirs[[names(strip)]]
  axes  <- purrr::map_chr(dirs, filter_elems, type = "axis")
  gtable::gtable_filter(table_layout, paste(c(panel_name, axes, strip), collapse = "|"))
}


facet_multi_col <- function(facets, layout, scales = "fixed", space = "fixed",
                            shrink = TRUE, labeller = "label_value",
                            drop = TRUE, strip.position = "top", 
                            min_prop = ifelse(strip.position %in% c("top", "bottom"), 
                                              0.12, 0.1)) {
  space <- match.arg(space, c("free", "fixed"))
  if (space == "free") {
    ## if we ask for free space we need scales everywhere, so make sure they are included
    scales <- "free"
  }
  facet <- facet_wrap(facets, ncol = 1, scales = scales, shrink = shrink, 
                      labeller = labeller, drop = drop, strip.position = strip.position)
  params <- facet$params
  params$space_free <- space == "free"
  params$layout <- layout
  params$parent <- facet
  params$min_prop <- min_prop
  ggproto(NULL, FacetMultiCol, shrink = shrink, params = params)
}



render <- function(self, panels, layout, 
                   x_scales, y_scales, ranges, 
                   coord, data, theme, params) {
  combined <- ggproto_parent(FacetWrap, self)$draw_panels(panels, layout, 
                                                          x_scales, y_scales, ranges, 
                                                          coord, data, theme, params)
  if (params$space_free) {
    panel_names <- combined$layout$name
    panels <- lapply(panel_names[grepl("panel", panel_names)],
                     get_whole_panel,
                     table_layout = combined)
    
    ## remove zeroGrob panels
    zG <- sapply(panels, function(tg) all(sapply(tg$grobs, ggplot2:::is.zero)))
    panels <- panels[!zG]
    ## calculate height for each panel
    heights <- matrix(NA, NROW(params$layout), NCOL(params$layout))
    ## store the rounded range in the matrix cell corresponding to its position
    ## allow for a minimum space in dependence of the overall number of rows to
    ## render small panels well
    
    heights[as.matrix(layout[, c("ROW", "COL")])] <- vapply(ranges, function(r) 
      round(diff(r$y.range), 0), numeric(1))
    
    ## 12% should be the minimum height used by any panel if strip is on top otherwise 10%
    ## these values are empirical and can be changed
    min_height <- round(params$min_prop * max(colSums(heights, TRUE)), 0)
    heights[heights < min_height] <- min_height
    idx <- c(heights)
    idx[!is.na(idx)] <- seq_along(idx[!is.na(idx)])
    max_height_sum <- max(colSums(heights, na.rm = TRUE))
    i <- 0
    layout_matrix <- apply(heights, 2, function(col) {
      res <- unlist(lapply(col, function(n) {
        i <<- i + 1
        mark <- idx[i]
        if (is.na(n)) {
          NA
        } else {
          rep(mark, n)
        }
      }))
      len <- length(res)
      if (len < max_height_sum) {
        res <- c(res, rep(NA, max_height_sum - len))
      }
      res
    })
    
    ## set width of left axis to maximum width to align plots
    max_width <- max(do.call(grid::unit.c, lapply(panels, function(gt) gt$widths[1]))) + grid::unit(1, "pt") ## Add margin <--------
    panels <- lapply(panels, function(p) { #                                                                      between the two
      p$widths[1] <- max_width
      p
    })
    
    combined <- gridExtra::arrangeGrob(grobs = panels,
                                       layout_matrix = layout_matrix,
                                       as.table = FALSE)
    ## add name, such that find_panel can find the plotting area
    combined$layout$name <- paste("panel_", layout$LAB)
  }
  combined
}

layout <- function(data, params) {
  parent_layout <- params$parent$compute_layout(data, params)
  msg <- paste0("invalid ",
                sQuote("layout"),
                ". Falling back to ",
                sQuote("facet_wrap"),
                " layout")
  if (is.null(params$layout) ||
      !is.matrix(params$layout)) {
    warning(msg)
    parent_layout
  } else {
    ## smash layout into vector and remove NAs all done by sort
    layout <- params$layout
    panel_numbers <- sort(layout)
    if (!isTRUE(all.equal(sort(as.numeric(as.character(parent_layout$PANEL))),
                          panel_numbers))) {
      warning(msg)
      parent_layout
    } else {
      ## all good
      indices <- cbind(ROW = c(row(layout)),
                       COL = c(col(layout)),
                       PANEL = c(layout))
      indices <- indices[!is.na(indices[, "PANEL"]), ]
      ## delete row and col number from parent layout
      parent_layout$ROW <- parent_layout$COL <- NULL
      new_layout <- merge(parent_layout, 
                          indices,
                          by = "PANEL") %>%
        dplyr::arrange(PANEL)
      new_layout$PANEL <- factor(new_layout$PANEL)
      labs <- new_layout %>%
        dplyr::select(-PANEL,
                      -SCALE_X,
                      -SCALE_Y,
                      -ROW,
                      -COL) %>%
        dplyr::mutate(sep = "_") %>%
        do.call(paste, .)
      new_layout$LAB <- labs
      new_layout
      
      
    }
  }
}

FacetMultiCol <- ggproto("FacetMultiCol", FacetWrap,
                         compute_layout = layout,
                         draw_panels    = render)

################################################################################
################################################################################

### Layout code
# 01 "Anthropometric" 9 (1.65)
# 02 "Cancer" 3
# 03 "Cardiovascular" 8 (1.45)
# 04 "Endocrine & Metabolic" 7
# 05 "Hematological" 12 (2.20)
# 06 "Hepatic" 7
# 07 "Life history" 7
# 08 "Lifestyle" 14 (2.60)
# 09 "Musculoskeletal & Connective" 6
# 10 "Neuropsychiatric" 5
# 11 "Other diseases" 3
# 12 "Pulmonary" 5
# 13 "Renal"  4
# 14 "Reproductive" 7
# 15 "Serum Lipids" 6
# I have 11 elements of 1.25cm, 1 of 1.45cm, 1 of 1.65cm, 1 of 2.20cm and 1 of 2.60cm
# 3  3  4  5  5  6  6  7  7  7  7  8  9 12 14

# my_layout1 <- matrix(c(7, 6, 14, 4, 2, 12, 3, 5, 10, 1, 15, 8, 9, 13, 11, NA), nrow = 8, ncol = 2, byrow = FALSE) # unordered

my_layout1 <- matrix(c(7, 8, 1, 14, 15, 13, 5, 6, 3, 12, 2, 9, 10, 4, 11, NA), nrow = 8, ncol = 2, byrow = FALSE)

# my_layout1 <- matrix(c(1:14), nrow = 7, ncol = 2, byrow = FALSE) # Only causal evidence


my_colors <- c("black","#CC3311","#6699CC")
plot_df$analysis <- factor(plot_df$analysis, levels = c("Linear regression", "Forward MR", "Reverse MR"))


# Create the forest plot
global_tile <- forestplot(
  df = plot_df,
  name = description,
  estimate = b,
  pvalue = pval,
  psignif = 0.05/indep_trait,
  xlab = expression(""*alpha*" or "*beta*" (95% CI)"),
  colour = analysis,
  se = se
) +
  coord_cartesian(xlim = c(-0.3, 0.3)) +
  labs(color=NULL) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 17),
        strip.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        legend.position = "bottom",
        legend.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(size = 3.5))) +
  scale_color_manual(values = my_colors, # Use custom colors
                     labels=c(expression("Linear regression ("*beta*")"), 
                              expression("Telomere on trait MR ("*alpha*")"), 
                              expression("Trait on telomere MR ("*alpha*")")))
  


pdf(file = "/SET/PATH/TO/DIRECTORY/plot1.pdf",
    height = 20,
    width = 17)

plot_2_col <- global_tile + facet_multi_col("Category", my_layout1, scales = "free_y", 
                              space = "free", strip.position = "top")
dev.off()

ggsave(file.path("/SET/PATH/TO/DIRECTORY/", "plot2.svg"),
       plot_2_col, width=17, height=20, units="in")

######

pdf(file = "/SET/PATH/TO/DIRECTORY/plot3.pdf",
    height = 40,
    width = 8)

forestplot(
  df = plot_df,
  name = description,
  estimate = b,
  pvalue = pval,
  psignif = 0.05/indep_trait,
  xlab = expression(""*alpha*" or "*beta*" (95% CI)"),
  colour = analysis,
  #shape = analysis,
  se = se,
) + 
  coord_cartesian(xlim = c(-0.3, 0.3)) +
  labs(color=NULL) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  facet_col(~Category, scales = "free", space = "free") +
  scale_color_manual(values = my_colors, # Use custom colors
                     labels=c(expression("Linear regression ("*beta*")"), 
                              expression("Telomere on trait MR ("*alpha*")"), 
                              expression("Trait on telomere MR ("*alpha*")")))

dev.off()

