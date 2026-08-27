# Interactive, result-grounded visual summaries shared by the Shiny Results
# Center and the self-contained Combined HTML Report.

visual_numeric_column <- function(data, candidates) {
  if (is.null(data) || !ncol(data)) return(NULL)
  normalized <- tolower(gsub("[^a-z0-9]", "", names(data)))
  for (candidate in candidates) {
    index <- match(tolower(gsub("[^a-z0-9]", "", candidate)), normalized)
    if (!is.na(index)) {
      values <- suppressWarnings(as.numeric(as.character(data[[index]])))
      if (any(is.finite(values))) return(names(data)[[index]])
    }
  }
  NULL
}

visual_clean_labels <- function(values, maximum = 72L) {
  values <- as.character(values)
  values[is.na(values)] <- "Unlabelled result"
  values <- trimws(gsub("[[:cntrl:]<>]+", " ", values))
  values[!nzchar(values)] <- "Unlabelled result"
  substr(values, 1L, maximum)
}

result_visual_summary <- function(data, agent_id, maximum_rows = 12L) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(data)) return(NULL)
  label_column <- result_label_column(data, agent_id)
  if (is.null(label_column) || !label_column %in% names(data)) return(NULL)
  labels <- visual_clean_labels(data[[label_column]])

  if (agent_id %in% c("gsva", "immune")) {
    numeric_columns <- setdiff(numeric_result_columns(data), label_column)
    if (!length(numeric_columns)) return(NULL)
    matrix_values <- vapply(numeric_columns, function(column) {
      suppressWarnings(as.numeric(as.character(data[[column]])))
    }, numeric(nrow(data)))
    if (is.null(dim(matrix_values))) matrix_values <- matrix(matrix_values, nrow = nrow(data))
    primary <- apply(matrix_values, 1L, stats::sd, na.rm = TRUE)
    secondary <- rowMeans(matrix_values, na.rm = TRUE)
    tertiary <- apply(matrix_values, 1L, function(values) diff(range(values, na.rm = TRUE)))
    metric <- "Across-sample standard deviation"
    secondary_label <- "Mean score"
    tertiary_label <- "Observed range"
  } else {
    p_column <- visual_numeric_column(data, c("p.adjust", "Adjusted.P.value", "qvalue", "FDR", "pvalue", "P.value"))
    primary_column <- switch(
      agent_id,
      chea = visual_numeric_column(data, c("Combined.Score", "Odds.Ratio")),
      string = visual_numeric_column(data, c("degree", "combined_score_prop", "combined_score")),
      drug = visual_numeric_column(data, c("Mean_Response", "Response", "AUC", "IC50", "sensitivity")),
      NULL
    )
    if (is.null(primary_column) && !is.null(p_column)) {
      probability <- suppressWarnings(as.numeric(as.character(data[[p_column]])))
      probability[!is.finite(probability) | probability <= 0] <- .Machine$double.xmin
      primary <- pmin(-log10(probability), 300)
      metric <- paste0("−log10(", p_column, ")")
    } else {
      primary_column <- primary_column %or_else% visual_numeric_column(
        data,
        c("FoldEnrichment", "RichFactor", "Count", "score", "Rank")
      )
      if (is.null(primary_column)) return(NULL)
      primary <- suppressWarnings(as.numeric(as.character(data[[primary_column]])))
      metric <- gsub("_", " ", primary_column)
    }

    secondary_column <- visual_numeric_column(data, c("Count", "Measurements", "degree", "FoldEnrichment", "Odds.Ratio"))
    tertiary_column <- visual_numeric_column(data, c("FoldEnrichment", "RichFactor", "Rank", "GeneRatio", "Combined.Score"))
    secondary <- if (is.null(secondary_column)) seq_len(nrow(data)) else suppressWarnings(as.numeric(as.character(data[[secondary_column]])))
    tertiary <- if (is.null(tertiary_column)) rank(primary, ties.method = "average", na.last = "keep") else suppressWarnings(as.numeric(as.character(data[[tertiary_column]])))
    secondary_label <- if (is.null(secondary_column)) "Result order" else gsub("_", " ", secondary_column)
    tertiary_label <- if (is.null(tertiary_column)) "Metric rank" else gsub("_", " ", tertiary_column)
  }

  keep <- is.finite(primary)
  if (!any(keep)) return(NULL)
  summary <- data.frame(
    label = labels[keep],
    primary = primary[keep],
    secondary = secondary[keep],
    tertiary = tertiary[keep],
    stringsAsFactors = FALSE
  )
  summary$secondary[!is.finite(summary$secondary)] <- 0
  summary$tertiary[!is.finite(summary$tertiary)] <- 0

  order_index <- if (identical(agent_id, "drug") && "Rank" %in% names(data)) {
    order(suppressWarnings(as.numeric(as.character(data$Rank[keep]))), na.last = TRUE)
  } else {
    order(summary$primary, decreasing = TRUE, na.last = TRUE)
  }
  summary <- summary[utils::head(order_index, maximum_rows), , drop = FALSE]
  summary <- summary[order(summary$primary, decreasing = FALSE), , drop = FALSE]

  list(
    data = summary,
    metric = metric,
    secondary_label = secondary_label,
    tertiary_label = tertiary_label,
    agent_id = agent_id
  )
}

professional_bar_plot <- function(summary, title) {
  if (is.null(summary) || !requireNamespace("plotly", quietly = TRUE)) return(NULL)
  data <- summary$data
  # R plotly passes named palettes through grDevices::col2rgb() when a
  # continuous colour variable is mapped. Plotly.js names such as "Tealgrn"
  # are therefore interpreted as literal R colours and fail at render time.
  palette <- if (identical(summary$agent_id, "drug")) {
    c("#443983", "#31688e", "#21918c", "#5ec962", "#fde725")
  } else {
    c("#e8f1eb", "#a9c7b1", "#5f9780", "#d98468", "#a94f3b")
  }
  colour_ramp <- grDevices::colorRampPalette(palette)(max(2L, nrow(data)))
  colour_rank <- rank(data$secondary, ties.method = "first", na.last = "keep")
  bar_colours <- colour_ramp[pmax(1L, pmin(length(colour_ramp), colour_rank))]
  plotly::plot_ly(
    data = data,
    x = ~primary,
    y = ~factor(label, levels = label),
    type = "bar",
    orientation = "h",
    hovertext = ~paste0(
      "<b>", label, "</b><br>", summary$metric, ": ", signif(primary, 4),
      "<br>", summary$secondary_label, ": ", signif(secondary, 4),
      "<br>", summary$tertiary_label, ": ", signif(tertiary, 4)
    ),
    hoverinfo = "text",
    marker = list(
      color = bar_colours,
      line = list(color = "rgba(36,51,45,.25)", width = 0.6)
    )
  ) |>
    plotly::layout(
      title = list(text = paste0("<b>", title, "</b><br><sup>Top result-wide features · hover for exact evidence</sup>"), x = 0.02),
      xaxis = list(title = summary$metric, zeroline = FALSE, gridcolor = "#e4e9e4"),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 225, r = 30, t = 74, b = 58),
      paper_bgcolor = "#fffdf8",
      plot_bgcolor = "#fffdf8",
      font = list(family = "-apple-system, BlinkMacSystemFont, Segoe UI, sans-serif", color = "#24332d"),
      showlegend = FALSE
    ) |>
    plotly::config(displaylogo = FALSE, responsive = TRUE, scrollZoom = TRUE)
}

rescale_visual_values <- function(values, to = c(8, 18)) {
  values <- suppressWarnings(as.numeric(values))
  values[!is.finite(values)] <- 0
  observed <- range(values)
  if (!diff(observed)) return(rep(mean(to), length(values)))
  to[[1L]] + (values - observed[[1L]]) / diff(observed) * diff(to)
}


string_network_3d_plot <- function(summary, interactions, title = "STRING") {
  required <- c("from_name", "to_name", "combined_score")
  if (
    is.null(interactions) || !nrow(interactions) ||
    !all(required %in% names(interactions)) ||
    !requireNamespace("igraph", quietly = TRUE)
  ) return(NULL)

  nodes <- summary$data[order(summary$data$primary, decreasing = TRUE), , drop = FALSE]
  node_names <- unique(nodes$label)
  edges <- interactions[
    interactions$from_name %in% node_names &
      interactions$to_name %in% node_names &
      interactions$from_name != interactions$to_name,
    required,
    drop = FALSE
  ]
  edges$combined_score <- suppressWarnings(as.numeric(edges$combined_score))
  edges <- edges[is.finite(edges$combined_score), , drop = FALSE]
  if (!nrow(edges)) return(NULL)
  edge_key <- vapply(seq_len(nrow(edges)), function(index) {
    paste(sort(c(edges$from_name[[index]], edges$to_name[[index]])), collapse = "::")
  }, character(1))
  edges <- edges[!duplicated(edge_key), , drop = FALSE]

  graph <- igraph::graph_from_data_frame(
    edges[, c("from_name", "to_name", "combined_score")],
    directed = FALSE,
    vertices = data.frame(name = node_names, stringsAsFactors = FALSE)
  )
  score <- suppressWarnings(as.numeric(unlist(
    igraph::edge_attr(graph, "combined_score"),
    use.names = FALSE
  )))

  # Preserve a maximum-confidence spanning tree so every displayed node remains
  # connected, then add the strongest remaining associations without drawing
  # the complete dense graph over itself.
  spanning <- igraph::mst(graph, weights = -score)
  spanning_keys <- apply(igraph::as_edgelist(spanning, names = TRUE), 1L, function(pair) {
    paste(sort(pair), collapse = "::")
  })
  graph_edges <- igraph::as_edgelist(graph, names = TRUE)
  graph_keys <- apply(graph_edges, 1L, function(pair) paste(sort(pair), collapse = "::"))
  strongest <- utils::head(order(score, decreasing = TRUE), 24L)
  keep_edges <- which(graph_keys %in% unique(c(spanning_keys, graph_keys[strongest])))

  set.seed(1701)
  coordinates <- igraph::layout_with_fr(graph, dim = 3L, weights = pmax(score, 1), niter = 700L)
  coordinates <- igraph::norm_coords(coordinates, xmin = -1, xmax = 1, ymin = -1, ymax = 1, zmin = -1, zmax = 1)
  rownames(coordinates) <- igraph::V(graph)$name

  line_x <- line_y <- line_z <- numeric()
  line_hover <- character()
  for (edge_index in keep_edges) {
    pair <- graph_edges[edge_index, ]
    edge_description <- paste0(
      "<b>", pair[[1L]], " ↔ ", pair[[2L]], "</b><br>",
      "STRING combined score: ", signif(score[[edge_index]], 5)
    )
    line_x <- c(line_x, coordinates[pair[[1L]], 1L], coordinates[pair[[2L]], 1L], NA_real_)
    line_y <- c(line_y, coordinates[pair[[1L]], 2L], coordinates[pair[[2L]], 2L], NA_real_)
    line_z <- c(line_z, coordinates[pair[[1L]], 3L], coordinates[pair[[2L]], 3L], NA_real_)
    line_hover <- c(line_hover, edge_description, edge_description, "")
  }

  node_degree <- nodes$primary[match(igraph::V(graph)$name, nodes$label)]
  node_hover <- paste0(
    "<b>", igraph::V(graph)$name, "</b><br>",
    summary$metric, ": ", signif(node_degree, 5),
    "<br>Displayed network: ", length(keep_edges), " high-confidence STRING associations"
  )

  plotly::plot_ly() |>
    plotly::add_trace(
      x = line_x, y = line_y, z = line_z,
      type = "scatter3d", mode = "lines",
      line = list(color = "rgba(94,124,109,0.72)", width = 5),
      hovertext = line_hover, hoverinfo = "text", name = "STRING association"
    ) |>
    plotly::add_trace(
      x = coordinates[, 1L], y = coordinates[, 2L], z = coordinates[, 3L],
      type = "scatter3d", mode = "markers+text",
      text = igraph::V(graph)$name, textposition = "top center",
      hovertext = node_hover, hoverinfo = "text",
      marker = list(
        size = rescale_visual_values(node_degree, c(12, 24)),
        color = node_degree,
        colorscale = list(c(0, "#bdd4c5"), c(0.55, "#d98468"), c(1, "#a94f3b")),
        line = list(color = "#fffdf8", width = 1.2),
        opacity = 0.98,
        showscale = TRUE,
        colorbar = list(title = summary$metric, thickness = 13)
      ),
      textfont = list(size = 11, color = "#24332d"),
      name = "Hub protein"
    ) |>
    plotly::layout(
      title = list(
        text = paste0(
          "<b>", title, " · connected 3D interaction network</b><br>",
          "<sup>Balls are ranked hub proteins · lines are retrieved STRING associations · hover for exact degree</sup>"
        ),
        x = 0.02
      ),
      scene = list(
        xaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
        yaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
        zaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE),
        aspectmode = "cube",
        camera = list(eye = list(x = 1.35, y = 1.45, z = 1.05)),
        bgcolor = "#fffdf8"
      ),
      margin = list(l = 0, r = 25, t = 88, b = 0),
      paper_bgcolor = "#fffdf8",
      font = list(family = "-apple-system, BlinkMacSystemFont, Segoe UI, sans-serif", color = "#24332d"),
      showlegend = FALSE
    ) |>
    plotly::config(displaylogo = FALSE, responsive = TRUE, scrollZoom = TRUE)
}


evidence_3d_plot <- function(summary, title, interactions = NULL) {
  if (is.null(summary) || nrow(summary$data) < 3L || !requireNamespace("plotly", quietly = TRUE)) return(NULL)
  if (!identical(summary$agent_id, "string")) return(NULL)
  string_network_3d_plot(summary, interactions, title)
}

string_network_guide_ui <- function(summary, interactions = NULL) {
  if (is.null(summary) || !identical(summary$agent_id, "string")) return(NULL)
  ranked <- summary$data[order(summary$data$primary, decreasing = TRUE), , drop = FALSE]
  top <- utils::head(ranked, 6L)
  eligible_edges <- 0L
  if (!is.null(interactions) && nrow(interactions) && all(c("from_name", "to_name") %in% names(interactions))) {
    node_names <- as.character(ranked$label)
    edges <- interactions[
      as.character(interactions$from_name) %in% node_names &
        as.character(interactions$to_name) %in% node_names &
        as.character(interactions$from_name) != as.character(interactions$to_name),
      , drop = FALSE
    ]
    if (nrow(edges)) {
      edge_key <- vapply(seq_len(nrow(edges)), function(index) {
        paste(sort(c(as.character(edges$from_name[[index]]), as.character(edges$to_name[[index]]))), collapse = "::")
      }, character(1))
      eligible_edges <- length(unique(edge_key))
    }
  }

  shiny::tags$aside(
    class = "string-app-network-guide",
    `aria-label` = "STRING 3D network guide",
    shiny::span(class = "network-guide-kicker", "NETWORK GUIDE"),
    shiny::h4("How to read this view"),
    shiny::div(
      class = "string-app-network-stats",
      shiny::div(shiny::strong(nrow(ranked)), shiny::span("ranked nodes considered")),
      shiny::div(shiny::strong(eligible_edges), shiny::span("eligible retrieved edges"))
    ),
    shiny::tags$ul(
      class = "string-app-network-legend",
      shiny::tags$li(shiny::i(class = "legend-node"), shiny::span(shiny::strong("Balls"), " are submitted proteins ranked by STRING interaction degree.")),
      shiny::tags$li(shiny::i(class = "legend-edge"), shiny::span(shiny::strong("Lines"), " are retrieved STRING associations, not directional regulatory arrows.")),
      shiny::tags$li(shiny::i(class = "legend-size"), shiny::span(shiny::strong("Size and colour"), " encode interaction degree within the retrieved network."))
    ),
    shiny::h5("Highest-degree hubs"),
    shiny::tags$ol(
      class = "string-app-hub-list",
      lapply(seq_len(nrow(top)), function(index) {
        shiny::tags$li(shiny::span(top$label[[index]]), shiny::strong(paste(signif(top$primary[[index]], 5L), "degree")))
      })
    ),
    shiny::h5("Mouse controls"),
    shiny::p(class = "string-app-network-controls", shiny::strong("Drag"), " to rotate · ", shiny::strong("scroll"), " to zoom · ", shiny::strong("hover"), " for exact values · ", shiny::strong("double-click"), " to reset."),
    shiny::p(class = "string-app-network-warning", shiny::strong("Interpret carefully: "), "connectivity prioritizes candidates; it does not prove activity, functional importance, direction, or causality.")
  )
}
