# First we plot the time-scaled maximum clade credibility (MCC) tree as Figure 2b (under Figure 2a, bottom, output width: 10, height: 8); then we plot individual clades (serveral subplots) under the ML tree with bootstrap support as supplementary figures.

# Load required libraries
library(tidyverse)
library(treeio)
library(ggtree)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(phytools)
library(patchwork)

# Supplementary Figures

## Read wards metadata
metadata_hai <- readxl::read_xlsx("data/hospital_data/HAI-28MAY2022-16AUG2022.xlsx")
metadata_cases <- readxl::read_xlsx("data/hospital_data/HK_nosocomial_metadata_2025-09-13_GS_031025.xlsx")
metadata_cases$new_name <- sapply(metadata_cases$FASTA, function(x) {
  this_name <- strsplit(x, "/", fixed = TRUE)[[1]][1:3] %>% paste(collapse = "_")
  this_name
})
df_Hospital_anonymized <- read_csv("data/hospital_data/Hospital_anonymized.csv")
metadata_cases <- left_join(metadata_cases, df_Hospital_anonymized, by = c("HOSPITAL" = "Hospital"))
metadata_hai <- left_join(metadata_hai, df_Hospital_anonymized, by = c("HOSPITAL" = "Hospital")) 
write_csv(metadata_hai, "data/hospital_data/metadata_hai.csv")

unique_hospital_wards <- paste0(metadata_hai$Hospital_anonymized, "-", metadata_hai$`WARD/CLUSTERS`)
unique_hospital_wards_official <- unique_hospital_wards[metadata_hai$`Confirmed cluster`]
unique_hospitals <- metadata_hai$Hospital_anonymized %>% unique()
unique_hospitals_official <- metadata_hai$Hospital_anonymized[metadata_hai$`Confirmed cluster`] %>% unique()

(tmp <- metadata_hai %>% filter(`Confirmed cluster`) %>% select(Hospital_anonymized, `WARD/CLUSTERS`, Lineages) %>% arrange(Lineages, Hospital_anonymized))
length(unique(tmp$Hospital_anonymized))

source("scripts/helper/ward_colors.R")

## Read ML tree for BA.2.2
tree_ml_ba22 <- treeio::read.tree("results/trees/ba22_gtr+g.treefile")
tree_ml_ba56 <- treeio::read.tree("results/trees/ba56_gtr+g.treefile")
tree_ml_ba2121 <- treeio::read.tree("results/trees/ba2121_gtr+g.treefile")

fig_ml_trees_hospitals <- lapply(unique_hospitals, function(this_hospital) {
  # this_hospital = unique_hospitals[10]
  these_wards <- metadata_hai$`WARD/CLUSTERS`[metadata_hai$Hospital_anonymized == this_hospital]
  print(this_hospital)
  print(these_wards)
  if (this_hospital %in% unique_hospitals_official) {
    this_plot_title <- this_hospital
  } else {
    this_plot_title <- paste0(this_hospital, "*") # add asterisk for non-officially reported nosocimial infection hospitals
  }
  these_lineages <- metadata_cases %>%
    filter(Hospital_anonymized == this_hospital) %>%
    .$LINEAGE %>%
    unique()

  lapply(these_lineages, function(this_lineage) {
    # this_lineage <- these_lineages[1]
    print(this_lineage)
    lapply(c(TRUE, FALSE), function(check_only_official_clusters) {
      # check_only_official_clusters = TRUE
      print(paste0("Check only official clusters: ", check_only_official_clusters))
      this_tree <- switch(this_lineage,
        "BA.2.2" = tree_ml_ba22,
        "BA.5.6" = tree_ml_ba56,
        "BA.2.12.1" = tree_ml_ba2121
      )
      metadata_this_hospital_lineage <- metadata_cases %>% filter(Hospital_anonymized == this_hospital, LINEAGE == this_lineage)
      if (check_only_official_clusters) {
        metadata_this_hospital_lineage <- metadata_this_hospital_lineage %>% filter(paste0(Hospital_anonymized, "-", `WARD/CLUSTERS`) %in% unique_hospital_wards_official)
      }
      if (nrow(metadata_this_hospital_lineage) == 0) {
        print(paste0("No officially reported data for ", this_hospital, " in ", this_lineage))
        return(list(NULL))
      }
      these_tip_names <- metadata_this_hospital_lineage$FASTA %>% list()
      names(these_tip_names) <- this_hospital
      if(length(unlist(these_tip_names)) > 1){
        this_MRCA_node <- ape::getMRCA(this_tree, unlist(these_tip_names))
      } else {
        # output the parent node of the tip
        tip_index <- which(this_tree$tip.label == unlist(these_tip_names))
        this_MRCA_node <- this_tree$edge[this_tree$edge[,2] == tip_index, 1]
      }
      this_clade <- ape::extract.clade(this_tree, this_MRCA_node)
      this_clade <- groupOTU(this_clade, these_tip_names)

      tips_sequenced_by_us <- this_clade$tip.label[this_clade$tip.label %in% metadata_cases$FASTA]

      # then check the descendants from the MRCA node, if specific clade doesnot contain any of the tips, then output the clade node and collapse later
      # Find all internal nodes within this extracted clade's phylo structure
      # Exclude the root node of this extracted clade itself
      internal_nodes_in_clade <- (Ntip(this_clade) + 2):(Ntip(this_clade) + Nnode(this_clade))

      nodes_to_collapse <- c()
      # Define the target tips for the current clade
      target_tips_in_clade <- unlist(these_tip_names)

      if (Nnode(this_clade) > 1) { # Only proceed if there are internal nodes other than the root
        processed_nodes <- c() # Keep track of nodes already processed or part of a collapsed clade

        for (node_idx in internal_nodes_in_clade) {
          # Skip if this node is already a descendant of a node marked for collapse
          is_descendant_of_collapsed <- FALSE
          for (collapsed_node in nodes_to_collapse) {
            # Check if node_idx is a descendant of collapsed_node
            # Need error handling in case getDescendants fails
            descendants_of_collapsed <- tryCatch(
              phytools::getDescendants(this_clade, collapsed_node),
              error = function(e) NULL
            )
            if (!is.null(descendants_of_collapsed) && node_idx %in% descendants_of_collapsed) {
              is_descendant_of_collapsed <- TRUE
              break # No need to check other collapsed nodes
            }
          }
          if (is_descendant_of_collapsed) {
            next # Skip to the next node_idx
          }

          # Get tips descending from this internal node
          subclade_tips <- tryCatch(
            {
              phytools::getDescendants(this_clade, node_idx) %>% # Get all descendant nodes
                .[. <= Ntip(this_clade)] %>% # Filter for tip nodes
                this_clade$tip.label[.] # Get tip labels
            },
            error = function(e) {
              warning(paste("Could not get tips for node", node_idx, "in clade", x, ":", e$message))
              NULL
            }
          )

          if (!is.null(subclade_tips) && length(subclade_tips) > 0) {
            # Check if any of the tips in this subclade are among the target tips for this ward/cluster
            has_target_tips <- any(subclade_tips %in% target_tips_in_clade)

            # If the subclade contains NO target tips, mark its root node for collapsing
            if (!has_target_tips) {
              # Node numbers in this_clade should correspond to node numbers in the ggtree plot
              nodes_to_collapse <- c(nodes_to_collapse, node_idx)
              # Add all descendants of this node to processed_nodes to avoid checking them later
              # This might be redundant given the check at the beginning of the loop, but can be kept for clarity
              # all_descendants <- tryCatch(phytools::getDescendants(this_clade, node_idx), error = function(e) NULL)
              # if (!is.null(all_descendants)) {
              #    processed_nodes <- unique(c(processed_nodes, all_descendants))
              # }
            }
          }
          # Mark this node as processed
          # processed_nodes <- c(processed_nodes, node_idx) # Redundant if the initial check works correctly
        }
      }

      # Check which nodes contains highest number of tips
      # Calculate the number of tips for each node to be collapsed
      if (length(nodes_to_collapse) > 0) {
        num_tips_per_node <- sapply(nodes_to_collapse, function(node_id) {
          # Get all descendants (tips and internal nodes)
          descendants <- phytools::getDescendants(this_clade, node_id)
          # Count how many descendants are tips (indices <= Ntip)
          num_tips <- sum(descendants <= Ntip(this_clade))
          return(num_tips)
        })
        # Optional: Name the vector elements with the node IDs for clarity
        names(num_tips_per_node) <- nodes_to_collapse
      } else {
        print(paste("No nodes identified for collapsing in clade", this_hospital))
        num_tips_per_node <- integer(0) # Create an empty integer vector if no nodes
      }

      # Filter out nodes with fewer than 20 tips
      num_tips_per_node_gt_treshd <- num_tips_per_node[num_tips_per_node >= 20]
      # reorder the nodes_to_collapse based on the number of tips
      nodes_to_collapse <- names(num_tips_per_node_gt_treshd) %>% as.numeric()

      tree_tmp <- ggtree(this_clade, alpha = 0.8, size=0.2) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) # Increase upper expansion factor
      if (length(nodes_to_collapse) > 0) {
        for (i in 1:length(nodes_to_collapse)) {
          this_node_to_collapse <- nodes_to_collapse[i]
          tree_tmp <- ggtree::collapse(tree_tmp, this_node_to_collapse)
        }
      }
      tree_tmp <- tree_tmp +
        geom_point2(aes(subset = (node %in% nodes_to_collapse)), shape = 23, size = 5, fill = "black", alpha = 0.6) +
        # highlight tips in these_tip_names
        NULL

      # --- Node Labeling Logic --- # no need to label all the internal nodes, only the ones subtending >= 2 target tips
      phylo_obj <- ape::as.phylo(this_clade) # Use as.phylo to be safe

      # Get target tip labels for the current hospital/group
      # 'these_tip_names' is a list with one element named after the hospital
      target_tip_labels <- unlist(these_tip_names)

      # Get node numbers for target tips within the phylo object
      # Node numbers for tips are 1 to Ntip
      target_tip_nodes <- match(target_tip_labels, phylo_obj$tip.label)
      target_tip_nodes <- target_tip_nodes[!is.na(target_tip_nodes)] # Remove NAs if any tips weren't found in the extracted clade

      # Get all internal node numbers from the phylo object
      # Node numbers for internal nodes are Ntip+1 to Ntip+Nnode
      Ntip_phylo <- ape::Ntip(phylo_obj)
      Nnode_phylo <- ape::Nnode(phylo_obj)
      # Check if there are any internal nodes
      if (Nnode_phylo > 0) {
        internal_nodes <- (Ntip_phylo + 1):(Ntip_phylo + Nnode_phylo)
      } else {
        internal_nodes <- integer(0) # No internal nodes
      }

      # --- Step 1: Find internal nodes subtending >= 2 target tips AND are parent/grandparent ---
      candidate_nodes_list <- list()
      # Only proceed if there are at least 2 target tips and some internal nodes
      if (length(target_tip_nodes) >= 2 && length(internal_nodes) > 0) {
        edge_matrix <- phylo_obj$edge

        for (node_id in internal_nodes) {
          # Get all descendant node numbers using phytools::getDescendants
          descendant_nodes <- tryCatch(
            {
              phytools::getDescendants(phylo_obj, node_id)
            },
            error = function(e) {
              # warning(paste("Could not get descendants for node", node_id, ":", e$message))
              integer(0) # Return empty if error
            }
          )

          # Filter for descendant tip nodes (node numbers <= Ntip)
          descendant_tip_nodes <- descendant_nodes[descendant_nodes <= Ntip_phylo]

          # Find which descendant tips are in our target group
          target_descendants_in_clade <- intersect(descendant_tip_nodes, target_tip_nodes)

          # Keep node if it subtends at least 2 target tips
          if (length(target_descendants_in_clade) >= 2) {
            # --- New Check: Is node parent or grandparent of any target tip? ---
            is_parent_or_grandparent <- FALSE

            # Find direct children
            children <- edge_matrix[edge_matrix[, 1] == node_id, 2]
            children_tips <- children[children <= Ntip_phylo]

            # Check if any direct children are target tips
            if (any(children_tips %in% target_descendants_in_clade)) {
              is_parent_or_grandparent <- TRUE
            } else {
              # Find grandchildren (only check children that are internal nodes)
              internal_children <- children[children > Ntip_phylo]
              if (length(internal_children) > 0) {
                grandchildren <- edge_matrix[edge_matrix[, 1] %in% internal_children, 2]
                grandchildren_tips <- grandchildren[grandchildren <= Ntip_phylo]
                # Check if any grandchildren are target tips
                if (any(grandchildren_tips %in% target_descendants_in_clade)) {
                  is_parent_or_grandparent <- TRUE
                }
              }
            }
            # --- End New Check ---

            # Store the node if it meets both criteria
            if (is_parent_or_grandparent) {
              # Store the node and its sorted target descendant tip nodes
              candidate_nodes_list[[as.character(node_id)]] <- sort(target_descendants_in_clade)
            }
          }
        }
      }

      # --- Step 2: Filter out ancestral nodes with identical target tip sets ---
      nodes_to_keep <- numeric(0) # Initialize as empty numeric vector
      candidate_nodes_list <- candidate_nodes_list[order(as.numeric(names(candidate_nodes_list)), decreasing = TRUE)]
      nodes_to_keep <- names(candidate_nodes_list)[!duplicated(candidate_nodes_list)]

      # --- Step 3: Create node_label_data for plotting ---
      # Filter the data from the ggtree object (tree_tmp) for the nodes we decided to keep
      # tree_tmp is the ggtree plot object created earlier, its $data contains node info
      if (length(nodes_to_keep) > 0) {
        node_label_data <- tree_tmp$data %>%
          dplyr::filter(node %in% nodes_to_keep)
      } else {
        # Create an empty tibble with the expected columns if no nodes are kept
        # Adjust column names/types if necessary based on what geom_label_repel expects
        # Need at least node, x, y, potentially label, isTip etc. Get structure from tree_tmp$data
        node_label_data <- tree_tmp$data[0, ] # Create empty df with same columns
        # Or more robustly:
        # node_label_data <- dplyr::slice(tree_tmp$data, 0)
        warning("No internal nodes identified for labeling for hospital: ", this_hospital)
      }
      # --- End Node Labeling Logic ---

      tip_label_data <- tree_tmp$data %>% filter((node %in% nodes_to_collapse) | (group == this_hospital & isTip))
      tip_label_data$new_label <- ""
      tip_label_data$new_label[tip_label_data$isTip] <- sapply(tip_label_data$label[tip_label_data$isTip], function(x) {
        this_name <- strsplit(x, "/", fixed = TRUE)[[1]][c(1,3)] %>% paste(collapse = "_")
        this_ward <- metadata_cases$`WARD/CLUSTERS`[metadata_cases$FASTA == x]
        if(!paste0(this_hospital, "-", this_ward) %in% unique_hospital_wards_official){
          this_ward <- paste0(this_ward, "*") # add asterisk for non-officially reported nosocimial infection hospitals
        }
        paste0(this_name, " (", this_ward, "-")
      })
      tip_label_data$new_label[!tip_label_data$isTip] <- paste0("Collpased (N=", num_tips_per_node_gt_treshd, ")")

      tip_label_data <- left_join(tip_label_data, metadata_cases[, c("FASTA", "case_type_updated")], by = c("label" = "FASTA"))
      tip_label_data$new_label[tip_label_data$isTip] <- paste0(tip_label_data$new_label[tip_label_data$isTip], tip_label_data$case_type_updated[tip_label_data$isTip], ")")

      tip_label_data$size <- ifelse(tip_label_data$isTip, 3, 2)
      tip_label_data$color <- sapply(seq_len(nrow(tip_label_data)), function(i) {
        if(!tip_label_data$isTip[i]){
          "#000000"
        } else {
          this_ward <- metadata_cases$`WARD/CLUSTERS`[metadata_cases$FASTA == tip_label_data$label[i]]
          return(ward_colors[paste0(this_hospital, "-", this_ward)])
        }
      })

      p_out <- tree_tmp +
        geom_tippoint(aes(color = color), size = 2, show.legend = FALSE, data = tip_label_data %>% filter(group == this_hospital), alpha = 0.8) +
        geom_label_repel(aes(label = new_label, color = color, size = size), nudge_x = 0.00001, hjust = "left", force = 0.5, direction="y", alpha = 0.8, data = tip_label_data, arrow = arrow(length = unit(0.01, "npc")), show.legend = FALSE, max.overlaps = Inf) +
        geom_label_repel(aes(label = label), nudge_x = -0.00002, hjust = "right", force = 0.5, direction="y", color = "darkblue", alpha = 0.8, data = node_label_data, arrow = arrow(length = unit(0.01, "npc"), type = "closed"), show.legend = FALSE, max.overlaps = Inf) +
        scale_color_identity() +
        scale_size_identity() +
        ggtitle(paste0(this_plot_title, " (", this_lineage, ")")) +
        NULL
      if(max(tree_tmp$data$x, na.rm=TRUE)>0){
        p_out <- p_out +
          geom_treescale(x = max(tree_tmp$data$x, na.rm = TRUE) * 0.8) +
          NULL
      }
      return(list(file_name = paste0(this_hospital, "_", this_lineage, "_", c("All", "Official")[check_only_official_clusters + 1]), plot = p_out))
    })
  })
})

# save the plots using ggsave, adjust the width and height per the actual plot
dir.create("results/figs_ml_trees_hospitals", showWarnings = FALSE)
sapply(1:length(fig_ml_trees_hospitals), function(i) {
  this_hospital_plot_list <- fig_ml_trees_hospitals[[i]]
  sapply(this_hospital_plot_list, function(this_lineage_plot_list){
    # this_lineage_plot_list <- this_hospital_plot_list[[1]]
    sapply(this_lineage_plot_list, function(this_plot_list) {
      # this_plot_list <- this_lineage_plot_list[[1]]
      if(is.null(unlist(this_plot_list))) return(NULL) # Skip if NULL
      this_plot_title <- this_plot_list$file_name
      print(this_plot_title)
      p <- this_plot_list$plot

      # Calculate data range directly from p$data
      x_range <- range(p$data$x, na.rm = TRUE)
      y_range <- range(p$data$y, na.rm = TRUE)
      x_diff <- diff(x_range)
      y_diff <- diff(y_range)

      # Handle cases where range might be zero or invalid
      if (!is.finite(x_diff) || x_diff <= 0) {
        x_diff <- 1 # Assign a default non-zero difference if needed
      }
      if (!is.finite(y_diff) || y_diff <= 0) {
        y_diff <- 1 # Assign a default non-zero difference if needed
      }

      # Define base size and scaling factors (adjust these as needed)
      base_width_cm <- 18
      base_height_cm <- 18
      scale_factor_width <- 5000 # Increase width based on branch length range
      scale_factor_height <- 0.5 # Increase height based on y-range (roughly tip count)
      max_width_cm <- 30
      max_height_cm <- 50

      # Calculate dynamic width and height
      plot_width <- max(base_width_cm, min(max_width_cm, base_width_cm + x_diff * scale_factor_width))
      plot_height <- max(base_height_cm, min(max_height_cm, base_height_cm + y_diff * scale_factor_height))

      # Handle potential NaN/Inf if ranges are zero or invalid
      if (!is.finite(plot_width)) plot_width <- base_width_cm
      if (!is.finite(plot_height)) plot_height <- base_height_cm

      # --- Adjust Label Sizes Based on Plot Height ---
      # Define reference size and dimension
      ref_height_cm <- 50 # Reference height in cm for base label sizes
      ref_tip_label_base_size <- 3 # Base size for tip labels at ref height
      ref_node_label_base_size <- 3 # Base size for internal node labels at ref height
      ref_collapsed_label_base_size <- 2.5 # Base size for collapsed node labels at ref height
      min_scale_factor <- 0.7 # Minimum scaling factor
      max_scale_factor <- 1.5 # Maximum scaling factor

      # Calculate scaling factor based on plot height relative to reference height
      scale_factor <- plot_height / ref_height_cm
      # Clamp the scaling factor to avoid excessively small or large labels
      scale_factor <- max(min_scale_factor, min(max_scale_factor, scale_factor))

      # Calculate scaled sizes
      scaled_tip_size <- ref_tip_label_base_size * scale_factor
      scaled_node_size <- ref_node_label_base_size * scale_factor
      scaled_collapsed_size <- ref_collapsed_label_base_size * scale_factor

      # Modify the plot object 'p' to update label sizes
      # Find the geom_label_repel layers dynamically
      layer_indices <- which(sapply(p$layers, function(layer) inherits(layer$geom, "GeomLabelRepel")))

      if (length(layer_indices) >= 2) {
        tip_label_layer_idx <- layer_indices[1] # Assumes first is tips/collapsed
        node_label_layer_idx <- layer_indices[2] # Assumes second is internal nodes

        # Modify the size data for the first geom_label_repel (tips and collapsed)
        # This layer uses scale_size_identity via aes(size=size)
        original_tip_data <- p$layers[[tip_label_layer_idx]]$data
        if ("isTip" %in% names(original_tip_data)) {
          # Recalculate size based on scaled values
          new_sizes <- ifelse(original_tip_data$isTip, scaled_tip_size, scaled_collapsed_size)
          p$layers[[tip_label_layer_idx]]$data$size <- new_sizes
        } else {
          warning("Column 'isTip' not found in data for the first GeomLabelRepel layer for tree: ", this_plot_title)
          # Fallback: apply an average scaled size or keep original
          p$layers[[tip_label_layer_idx]]$data$size <- (scaled_tip_size + scaled_collapsed_size) / 2
        }


        # Modify the size parameter for the second geom_label_repel (internal nodes)
        # This layer seems to use a default size, so we set aes_params$size
        # Initialize aes_params if it's NULL
        if (is.null(p$layers[[node_label_layer_idx]]$aes_params)) {
          p$layers[[node_label_layer_idx]]$aes_params <- list()
        }
        p$layers[[node_label_layer_idx]]$aes_params$size <- scaled_node_size
      } else {
        warning("Could not find two GeomLabelRepel layers to adjust size for tree: ", this_plot_title)
      }
      # --- End Label Size Adjustment ---


      # Construct filename
      filename <- paste0("results/figs_ml_trees_hospitals/fig_ml_tree_", gsub("[^A-Za-z0-9_]", "_", this_plot_title), ".pdf") # Sanitize filename

      # Save the plot
      ggsave(
        filename = filename,
        plot = p, # Save the modified plot object with adjusted label sizes
        width = plot_width,
        height = plot_height,
        units = "cm",
        limitsize = FALSE,
        device = cairo_pdf # Use cairo_pdf for better rendering, especially with transparency/alpha
      )
    })
  })
})

# Figure 2b
# Read BEAST MCC tree
tree <- read.beast("results/trees/ba22_phylogeo_edit_mcc.tre")
# Drop tip(s) as they are no longer classified as nosocomial infection
tip_labels <- ape::as.phylo(tree)$tip.label
to_drop <- grep("^22TM347527($|/)", tip_labels, value = TRUE)
if (length(to_drop) > 0) {
  tree <- treeio::drop.tip(tree, to_drop)
}
to_drop <- grep("^22MB358559($|/)", tip_labels, value = TRUE)
if (length(to_drop) > 0) {
  tree <- treeio::drop.tip(tree, to_drop)
}

data_collection_date <- read_tsv("results/trees/date_ba22.tsv")
ward_colors <- c(ward_colors[unique_hospital_wards_official], "Excluded clusters" = "#808080")

p <- ggtree(tree, mrsd = max(data_collection_date$date, na.rm = TRUE), alpha = 0.8, size=0.01)
time_range <- c(min(p$data$x), max(p$data$x))

# Convert the ggtree object to a phylo object so we can get descendants
phylo_obj <- as.phylo(tree)
Ntip_phylo <- ape::Ntip(phylo_obj)
Nnode_phylo <- ape::Nnode(phylo_obj)

# Identify all internal nodes (for a phylo object, tips are 1..Ntips, internal nodes are (Ntip+1)..(Ntip+Nnode))
internal_nodes <- (Ntip_phylo + 1):(Ntip_phylo + Nnode_phylo)

# Sort nodes from root to tips (descending node IDs) to prioritize collapsing higher-level nodes
internal_nodes <- sort(internal_nodes, decreasing = TRUE)

# We'll store nodes to collapse
nodes_to_collapse <- c()
# Track all descendants of nodes that will be collapsed
all_descendants_to_skip <- c()

# We'll also track node-to-tip counts so we can label collapsed clades with geom_label_repel
collapsed_clades_info <- data.frame(node = numeric(0), n_tips = numeric(0))

for (node_id in internal_nodes) {
  print(node_id)
  
  # Skip if this node is already a descendant of a node marked for collapsing
  if (node_id %in% all_descendants_to_skip) {
    next
  }
  
  all_descendants <- tryCatch(
    {
      phytools::getDescendants(phylo_obj, node_id)
    },
    error = function(e) {
      NULL
    }
  )
  if (is.null(all_descendants)) next
  
  descendant_tips <- all_descendants[all_descendants <= Ntip_phylo]
  tip_labels <- phylo_obj$tip.label[descendant_tips]
  desc_groups <- p$data %>%
    dplyr::filter(label %in% tip_labels) %>%
    dplyr::pull(group)

  # If none are "hospital" and number of tips >= 20, record for collapse
  if (!any(desc_groups == "hospital") && length(descendant_tips) >= 20) {
    nodes_to_collapse <- c(nodes_to_collapse, node_id)
    collapsed_clades_info <- rbind(
      collapsed_clades_info,
      data.frame(node = node_id, n_tips = length(descendant_tips))
    )
    
    # Add all descendants of this node to the skip list
    all_descendants_to_skip <- c(all_descendants_to_skip, all_descendants)
  }
}

# Collapse each node
if (length(nodes_to_collapse) > 0) {
  for (nid in nodes_to_collapse) {
    p <- ggtree::collapse(p, nid)
  }
}

# Join collapsed_clades_info with p$data to get node positions for labeling
collapsed_labels_df <- p$data %>%
  dplyr::inner_join(collapsed_clades_info, by = c("node" = "node"))

tree_data <- as_tibble(p$data)
tree_data <- left_join(tree_data, data_collection_date, by = c("label" = "strain"))
tree_data <- left_join(tree_data, metadata_cases[, c("FASTA", "Hospital_anonymized", "WARD/CLUSTERS")], by = c("label" = "FASTA")) %>%
  mutate(
    hospital_ward = paste0(Hospital_anonymized, "-", `WARD/CLUSTERS`)
)
tree_data$hospital_ward[tree_data$hospital_ward=="NA-NA"] <- NA # Remove NA 
tree_data$hospital_ward[!tree_data$hospital_ward %in% unique_hospital_wards_official] <- "Excluded clusters"
tree_data$hospital_ward <- factor(
  tree_data$hospital_ward,
  levels = c(sort(unique_hospital_wards_official), "Excluded clusters"),
  labels = c(sort(unique_hospital_wards_official), "Excluded clusters")
)

date_breaks <- lubridate::decimal_date(lubridate::ymd(
  c(
    "2021-10-01", "2022-01-01", "2022-04-01",
    "2022-07-01", "2022-10-01", "2023-01-01"
  )
))

p_2c <- p + 
  geom_tree(aes(color = group), linewidth = 0.15, alpha = 0.8) +
  ggrepel::geom_label_repel(
    data = collapsed_labels_df,
    aes(x = x, y = y, label = paste0("N=", n_tips)),
    color = "black",
    size = 1.5,
    label.padding = unit(0.15, "lines"),
    alpha = 0.8
  ) +
  geom_tippoint(
    aes(shape = "Cases under study", fill = hospital_ward),
    alpha = 0.8, stroke = 0.3,
    data = tree_data %>% filter(!is.na(Hospital_anonymized)),
    size = 2,
    show.legend = TRUE
  ) +
  geom_point2(
    aes(subset = (node %in% nodes_to_collapse), shape = "Collapsed community clades"),
    size = 2,
    alpha = 0.6,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    name = "Hospital-clusters",
    values = ward_colors,
    breaks = c(sort(unique_hospital_wards_official), "Excluded clusters"),
    guide = guide_legend(override.aes = list(shape = 21))
  ) +
  scale_shape_manual(
    name = "Node type",
    values = c("Cases under study" = 21, "Collapsed community clades" = 23)
  ) +
  scale_color_manual(
    name = "Branch states",
    values = c(
      "community" = "#606e67",
      "community+hospital" = "#2E8B57",
      "hospital" = "#8B0000",
      "hospital+community" = "#A0522D"
    ),
    guide = guide_legend(override.aes = list(shape = NA, linewidth = 1))
  ) +
  scale_x_continuous(
    name = "Date",
    labels = function(x) format(lubridate::date_decimal(x), "%Y-%m"),
    breaks = date_breaks,
    expand = c(0, 0)
  ) +
  theme_tree2() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  scale_y_continuous(expand = c(0.02, 0))

ggsave("results/fig_2c.pdf", p_2c, width = 20, height = 20, units = "cm", limitsize = FALSE)

# Figure 2a: upper panel, will be the plots showing population dynamics and markov jumps; x-axis: date, y-axis: eff estimates and MJ events.

## Plotting the effective population size estimates
data_skygrid <- read_tsv("results/ruixuan/Hospital infection analysis data_3/ba22_skygrid_reconstruction.tsv", skip = 1)
p_a1 <- ggplot(data_skygrid)+
  geom_line(aes(x = time, y = median), color = "#2e568b") +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper), fill = "#2e568b", alpha = 0.3) +
  theme_minimal() +
  scale_y_log10(
    name = "Effective population size",
    labels = scales::comma_format()
  ) +
  geom_vline(xintercept = as.Date("2022-05-28"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  geom_vline(xintercept = as.Date("2022-08-18"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  scale_x_continuous(
    name = "Date",
    labels = function(x) format(lubridate::date_decimal(x), "%Y-%m"),
    expand = c(0, 0),
    limits = time_range,
    breaks = date_breaks
  ) +
  NULL

## MJ
data_markov_jump <- read_csv("results/ruixuan/Hospital infection analysis data_2/ba22_mj_timing.csv")
# Filter transitions and normalize by 501 trees
community_to_hospital <- data_markov_jump %>% 
  filter(startLocation == "community", endLocation == "hospital") %>% 
  mutate(transition = "community_to_hospital")
hospital_to_community <- data_markov_jump %>% 
  filter(startLocation == "hospital", endLocation == "community") %>% 
  mutate(transition = "hospital_to_community")

# Combine and normalize
combined <- bind_rows(community_to_hospital, hospital_to_community) %>% 
  mutate(rate = 1/501)  # Each jump contributes 1/501
combined$date <- decimal2Date(combined$time)

# Create daily bins for 2022
all_days <- seq(decimal2Date(time_range)[1], decimal2Date(time_range)[2], by = "day")

# Create daily rates
daily_rates <- combined %>% 
  mutate(day = floor_date(date, "day")) %>% 
  group_by(day, transition) %>% 
  summarize(daily_rate = sum(rate), .groups = "drop") %>% 
  complete(day = all_days, transition, fill = list(daily_rate = 0)) %>% 
  group_by(transition) %>% 
  ungroup()

p_a2 <- ggplot(daily_rates, aes(x = day, y = daily_rate, fill = transition)) +
  geom_col(
    position = position_dodge(width = 1),  # Matches bar width
    width = 0.95, alpha = 0.8) +
  geom_vline(xintercept = as.Date("2022-05-28"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  geom_vline(xintercept = as.Date("2022-08-18"), linetype = "dotted", color = "black", size = 1, alpha = 0.8) +
  scale_fill_manual(
    values = c("community_to_hospital" = "#d95f02", 
               "hospital_to_community" = "#1b9e77"),
    labels = c("Community → Hospital", "Hospital → Community")
  ) +
  scale_x_date(
    breaks = decimal2Date(date_breaks),
    date_labels = "%Y-%m",
    limits = decimal2Date(time_range),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Date",
    y = "Markov Jumps (mean)",
    fill = "Transition Direction"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(3, "pt"),  # Slightly longer ticks
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)

fig_2 <- p_a2 + p_2c +
  plot_layout(ncol = 1, heights = c(1.5, 7)) +
  plot_annotation(tag_levels = "a") & 
  theme(legend.position='right')
ggsave("results/fig_2.pdf", fig_2, width = 8, height = 10, limitsize = FALSE)
