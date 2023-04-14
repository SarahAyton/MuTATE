#' PlotTree: Plot Decision Tree
#'
#' Plot a multi-target decision tree.
#'
#' @param tree A decision tree object.
#'
#' @return A plot of the multi-target decision tree.
#' @import igraph
#' @import ggraph
#' @import grDevices
#' @import graphics
#' @export

PlotTree <- function(tree) {
  # Extract parent and child node IDs
  tree[[1]]$parent[1] <- 0
  parent <- tree[[1]]$parent
  parent <- parent[complete.cases(parent)]
  child <- gsub(" *\\*", "", tree[[1]]$child)
  child <- child[complete.cases(child)]

  # Replace NA values with 0
  tree[[2]] <- lapply(tree[[2]], function(x) {
    x[is.na(x)] <- 0
    x
  })

  # create data frame with node data
  node_data <- bind_rows(lapply(unique(as.numeric(child)), function(i) {
    data.frame(parent = floor(i/2), node = i,
               Variable = ifelse(((i %in% unique(parent))), tree[[2]][[i]][["PartVar"]], paste0("Subtype")),
               Partition = tree[[2]][[i]][["SplitVar"]],
               Threshold = sub(paste0("^", tree[[2]][[i]][["Var"]], "\\s+"), "", tree[[2]][[i]][["SplitVar"]]),
               N = tree[[2]][[i]][["N"]],
               Pnode = tree[[2]][[i]][["Pnode"]], CP = tree[[2]][[i]][["CP"]],
               relerr = tree[[2]][[i]][["relerr"]], SXerror = tree[[2]][[i]][["SXerror"]])
  }))
  # Find the rows where Variable is "Subtype"
  subtype_rows <- rownames(node_data[which(node_data$Variable == "Subtype"), c(1:ncol(node_data))])
  # iterate over the unique subtypes and rename the corresponding rows
  for (i in 1:length(subtype_rows)) {
    node_data[subtype_rows[i],]$Variable <- paste0("Subtype ", i)
  }
  # Replace 0 with NA in CV and SXerror variables
  node_data$CP <- ifelse(node_data$CP == 0 & !(node_data$node %in% node_data$parent), NA, node_data$CP)
  node_data$SXerror <- ifelse(node_data$SXerror == 0, NA, node_data$SXerror)
  # Add nsplit and leaves columns to node data
  node_data <- node_data %>%
    mutate(nsplit = as.numeric(match(parent, unique(parent))-1),
           leaves = as.numeric(match(parent, unique(parent))))
  node_data$Pnode <- as.numeric(node_data$Pnode)
  node_data$N <- as.numeric(node_data$N)
  node_data$Label <- paste0(node_data$N, " (", round(node_data$Pnode * 100, 2), "%)")
  # round numeric values in continuous "Threshold" partitions to 2 decimal places
  suppressWarnings(node_data$Threshold <- ifelse(grepl("^[<>]=?\\s*[0-9.]+$", node_data$Threshold),
                                                 paste0(sub("^([<>]=?\\s*)([0-9.]+)$", "\\1", node_data$Threshold),
                                                        sprintf("%.2f", as.numeric(sub("^([<>]=?\\s*)([0-9.]+)$", "\\2", node_data$Threshold)))),
                                                 node_data$Threshold))

  # Summarize leaf target data
  # Extract target names
  leafids <- node_data[which(is.na(node_data$CP)), c(1:ncol(node_data))]$node
  leaf_data <- bind_rows(lapply(leafids, function(i) {
    targets <- names(tree[[2]][[i]]$Targets)
    output <- data.frame(Node = rep(i, length(targets)),
                         Target = targets,
                         Value = NA,
                         relerr = NA)
    for (t in 1:length(targets)) {
      if ("median" %in% names(tree[[2]][[i]]$Targets[[t]])) {
        val <- paste0("Median = ", round(tree[[2]][[i]]$Targets[[t]]$median, digits=2), " (IQR: ", paste0(round(tree[[2]][[i]]$Targets[[t]]$IQR25, digits=2), ", ", round(tree[[2]][[i]]$Targets[[t]]$IQR75, digits=2)), ")")
      }
      else if ("estimatedrate" %in% names(tree[[2]][[i]]$Targets[[t]])) {
        val <- paste0("Rate = ", round(tree[[2]][[i]]$Targets[[t]]$estimatedrate, digits=2), " (Events: ", paste0(tree[[2]][[i]]$Targets[[t]]$events_count, ")"))
      }
      else {
        val <- paste0("Predicted = ", tree[[2]][[i]]$Targets[[t]]$predicted, " (Freq: ", round((1-as.numeric(tree[[2]][[i]]$Targets[[t]]$expectedloss))*100), "%)")
      }
      output$Value[t] <- val
      output$relerr[t] <- tree[[2]][[i]]$Targets[[t]]$relerr
    }
    return(output)
  }))

  # Define a function to get the default color codes for a given number of categories
  get_default_colors <- function(n) {
    gg_colors <- scales::hue_pal()(n) # Get default ggplot2 colors
    rgb_cols <- sapply(gg_colors, function(x) {col2rgb(x)}) # Convert colors to RGB format
    rgb_cols <- t(rgb_cols) # Transpose to make RGB codes the rows
    return(apply(rgb_cols, 1, function(x) {paste0("#", paste0(as.hexmode(x), collapse = ""))})) # Convert RGB codes to hex format
  }
  leafcolors <- names(get_default_colors(nrow(node_data[which(is.na(node_data$CP)), c(1:ncol(node_data))])))
  # Create an igraph object from the tree data
  g <- graph.data.frame(node_data[2:nrow(node_data), 1:2], directed = TRUE) # V(g)$label
  # Generate colors based on node type:
  colrs <- c(rep("gray80",nrow(node_data[which(!is.na(node_data$CP)), c(1:ncol(node_data))])), leafcolors)
  V(g)$color <- colrs[V(g)]
  # Convert to a ggraph layout
  layout <- create_layout(g, "tree")
  # Increase the spacing between nodes
  #spacing_factor <- 5
  #layout[,1] <- layout[,1] * spacing_factor
  #layout[,2] <- layout[,2] * spacing_factor
  label_node <- sapply(as.numeric(layout$name), function(i) {
    if (i %in% leafids) {
      leafi <- leaf_data[leaf_data$Node == i, ]
      target_summary <- paste0("\n",('Outcomes'),"\n")
      for (t in 1:nrow(leafi)) {
        target_summary <- paste0(target_summary, leafi[t, "Target"], ": ", leafi[t, "Value"], "\n")
      }
      paste0("\n", node_data[node_data$node == i, "Variable"], "\n",
             node_data[node_data$node == i, "N"], " (",
             round(node_data[node_data$node == i, "Pnode"] * 100, 2), "%)\n", target_summary)
    }
    else {
      paste0(node_data[which(node_data$node==i), 1:ncol(node_data)]$Variable, "\n",
             node_data[which(node_data$node==i), 1:ncol(node_data)]$N, " (",
             round(node_data[which(node_data$node==i), 1:ncol(node_data)]$Pnode*100, 2), "%)")
    }
  })
  label_arrow <- sapply(2:nrow(node_data), function(i) {
    paste0(node_data[i, 1:ncol(node_data)]$Threshold)
  })
  # Function for plotting an elliptical node
  myellipse <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size <- 1/175 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    plotrix::draw.ellipse(x=coords[,1], y=coords[,2], a = vertex.size, b=vertex.size/2, col=vertex.color)
  }
  ## Register the shape with igraph
  add_shape("ellipse", clip=shapes("circle")$clip, plot=myellipse)
  # Generate shapes based on node type:
  nodeshapes <- c(rep("rectangle",nrow(node_data[which(!is.na(node_data$CP)), c(1:ncol(node_data))])),
                  rep("rectangle",nrow(node_data[which(is.na(node_data$CP)), c(1:ncol(node_data))])))
  V(g)$shape <- nodeshapes[V(g)]
  try <- plot(g, layout=layout.reingold.tilford)
  # Find the leaf nodes of the graph
  plot(g, layout = layout.reingold.tilford(g, root=1), rescale=FALSE, asp=1,
       ylim=c(min(layout$y)-0.5,max(layout$y)+0.5), xlim=c(min(layout$x)-0.5,max(layout$x)+0.5),
       vertex.size=((strwidth(label_node) + strwidth("oooo")) * 75),
       vertex.size2= strheight(label_node) * 2 * 75, #vertex.shape="ellipse",
       vertex.frame.color = "black", lineend = 'round',
       vertex.label=label_node, vertex.label.family='Times', vertex.label.cex=0.5, vertex.label.color= "black",
       edge.label=label_arrow, edge.label.family='Times', edge.label.cex=0.5, edge.label.color= "black",
       edge.arrow.size=0.15, arr.type="triangle", edge.arrow.color= "black",
       margin=c(0,0,0,0))
}
