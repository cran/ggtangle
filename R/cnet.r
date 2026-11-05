#' @rdname cnetplot
#' @importFrom igraph 'V<-'
#' @importFrom ggplot2 scale_size
#' @importFrom ggrepel geom_text_repel
#' @method cnetplot list
#' @details For list input, `showCategory` accepts either a single integer for
#'   the first `n` categories, a numeric vector of category indices, or a
#'   character vector of category names.
#' @param hilight_alpha transparent value for not selected to be highlight
#' @param ... additional parameters. One important parameter is 'curvature'
#'   (default is 0), which can be used to curve the edges (e.g.,
#'   `curvature = 0.2`). For circular layouts, this uses a cubic Bezier arc
#'   layer with the same radial control-point geometry as the former
#'   `ggraph::geom_edge_arc()` implementation. In circular mode, any non-zero
#'   curvature enables this arc layer; the arc geometry is intentionally kept
#'   compatible with the former backend. Legacy `enrichplot` arguments are
#'   also accepted: `colorEdge = TRUE` colors edges by category, and
#'   `circular = TRUE` selects the circular layout. A legacy
#'   `categorySize = "pvalue"` style argument is accepted as well: a single
#'   column name is evaluated the same way as `categorySizeBy`, provided the
#'   gene-set list carries numeric values for it (enrichment results attach
#'   them automatically).
#' @export
cnetplot.list <- function(
  x,
  layout = "nicely",
  showCategory = 5,
  color_category = "#E5C494",
  size_category = 1,
  color_item = "#B3B3B3",
  size_item = 1,
  color_edge = "grey",
  size_edge = .5,
  categorySizeBy = ~itemNum,
  node_label = "all",
  node_label_size = NULL,
  foldChange = NULL,
  fc_threshold = NULL,
  hilight = "none",
  hilight_alpha = .3,
  ...
) {
    dots <- list(...)
    ## Legacy `enrichplot`-era arguments. `colorEdge` / `circular` have no
    ## formals, so they arrive via `...`. `categorySize` has no formal either,
    ## but R partial-matches it onto `categorySizeBy`, where it is handled
    ## below; the `...` route is kept for exact named calls (e.g. `do.call`).
    legacy_category_size <- dots$categorySize
    dots$categorySize <- NULL

    if (isTRUE(dots$circular)) {
        layout <- "layout_circular"
    }
    dots$circular <- NULL

    color_edge_legacy <- dots$colorEdge
    dots$colorEdge <- NULL
    if (!is.null(color_edge_legacy) && identical(color_edge, "grey")) {
        if (isTRUE(color_edge_legacy)) {
            color_edge <- "category"
        } else if (rlang::is_string(color_edge_legacy)) {
            color_edge <- color_edge_legacy
        }
    }

    edge_layer_params <- list(linewidth = size_edge)
    if (!is.null(dots$curvature)) {
        edge_layer_params$curvature <- dots$curvature
        edge_layer_params$geom <- ggfun::geom_bezier
        dots$curvature <- NULL
    }

    ## `categorySizeBy` usage demonstrations
    ## x <- list(A = letters[1:10], B = letters[5:12])
    ## attr(x, "p.adjust") <- c(A = 0.01, B = 0.2)
    ## cnetplot(x, node_label = "none", categorySizeBy = ~ -log10(p.adjust))

    x_names <- names(x)
    x_attr_names <- setdiff(names(attributes(x)), c("names", "row.names", "class"))
    categoryAttr <- list()
    for (attr_name in x_attr_names) {
        attr_value <- attr(x, attr_name)
        if (!is.numeric(attr_value)) {
            next
        }
        if (!is.null(names(attr_value))) {
            attr_value <- attr_value[x_names]
        } else if (length(attr_value) == length(x)) {
            names(attr_value) <- x_names
        } else {
            next
        }
        categoryAttr[[attr_name]] <- attr_value
    }

    x <- subset_cnet_list(x, showCategory)
    originalItemNum <- setNames(vapply(x, length, FUN.VALUE = numeric(1)), names(x))

    # node_label <- match.arg(node_label, c("category", "all", "none", "item", "gene", "exclusive", "share"))
    if (length(node_label) > 1) {
        if (getOption("cnetplot_subset", default = FALSE)) {
            x <- subset_cnet_list_item(x, node_label)
            node_label <- "all"
        }
    } else if (
        !node_label %in%
            c("category", "all", "none", "item", "gene", "exclusive", "share")
    ) {
        if (!grepl("[><=]", node_label)) {
            stop("wrong parameter for 'node_label'")
        } else if (is.null(foldChange)) {
            stop(
                "'foldChange' should not be NULL with the 'node_label' setting"
            )
        }
    }

    if (length(node_label) == 1 && node_label == "gene") {
        node_label <- "item"
    }

    # Filter genes by absolute fold change threshold
    if (!is.null(fc_threshold) && !is.null(foldChange)) {
        keep_genes <- names(foldChange)[abs(foldChange) > fc_threshold]
        foldChange <- foldChange[keep_genes]
        x <- lapply(x, function(gs) {
            intersect(gs, keep_genes)
        })
        # Remove empty sets
        x <- x[sapply(x, length) > 0]
    }

    categorySizeQuo <- rlang::enquo(categorySizeBy)
    categorySizeExpr <- rlang::quo_get_expr(categorySizeQuo)
    if (rlang::is_formula(categorySizeExpr)) {
        categorySizeExpr <- rlang::f_rhs(categorySizeExpr)
    }
    ## The generic has no `categorySize` formal, so R binds a legacy
    ## `categorySize = "pvalue"` onto `categorySizeBy` by partial matching
    ## and the value shows up here as a character literal. Translate it to
    ## the size column, mirroring `enrichplot::cnetplot()` for enrichment
    ## results (YuLab-SMU/enrichplot#346).
    if (is.null(legacy_category_size) && rlang::is_string(categorySizeExpr)) {
        legacy_category_size <- categorySizeExpr
    }
    if (!is.null(legacy_category_size)) {
        categorySizeExpr <- translate_cnet_category_size(legacy_category_size)
    }
    exprIsItemNum <- is.symbol(categorySizeExpr) && as.character(categorySizeExpr) == "itemNum"

    categoryNames <- names(x)
    categoryDf <- data.frame(
        category = categoryNames,
        itemNum = vapply(x, length, FUN.VALUE = numeric(1)),
        stringsAsFactors = FALSE
    )
    if (length(categoryAttr) > 0) {
        for (attr_name in names(categoryAttr)) {
            categoryDf[[attr_name]] <- categoryAttr[[attr_name]][categoryNames]
        }
    }

    if (!is.null(legacy_category_size) && is.symbol(categorySizeExpr)) {
        legacy_col <- as.character(categorySizeExpr)
        if (!legacy_col %in% names(categoryDf)) {
            stop(
                sprintf("`categorySize = \"%s\"` requires a numeric `%s` value per category.",
                        as.character(legacy_category_size), legacy_col),
                "\nA plain gene-set list carries no statistics. Either pass the",
                " enrichment result itself (e.g. `cnetplot(enrichResult, ...)`),",
                " or attach the values as a named numeric attribute, e.g. ",
                sprintf("`attr(x, \"%s\") <- c(A = 0.01, B = 0.02)`.", legacy_col),
                call. = FALSE
            )
        }
    }

    categorySizeValues <- rlang::eval_tidy(categorySizeExpr, data = categoryDf)
    if (!is.numeric(categorySizeValues)) {
        stop("`categorySizeBy` must evaluate to a numeric vector.")
    }
    if (length(categorySizeValues) == 1) {
        categorySizeValues <- rep(categorySizeValues, nrow(categoryDf))
    }
    if (length(categorySizeValues) != nrow(categoryDf)) {
        stop("`categorySizeBy` must return a scalar or one value per category.")
    }
    if (anyNA(categorySizeValues)) {
        stop("`categorySizeBy` returned NA values.")
    }

    g <- list2graph(x)

    is_cat <- V(g)$.isCategory
    vals <- setNames(categorySizeValues, categoryNames)
    V(g)$size[is_cat] <- vals[V(g)$name[is_cat]]

    V(g)$.hilight <- 1
    if (all(hilight != "none")) {
        # maybe color the edge ?

        y <- subset_cnet_list(x, hilight)
        V(g)$.hilight <- hilight_alpha
        V(g)$.hilight[V(g)$name %in% names(y)] <- 1
        V(g)$.hilight[V(g)$name %in% unlist(y)] <- 1
    }

    # V(g)$.color <- color_item
    # V(g)$.color[1:length(x)] <- color_category

    if (!is.null(foldChange)) {
        V(g)$foldChange <- foldChange[V(g)$name]
        fc_mapping <- aes(color = foldChange, alpha = I(.data$.hilight))
    } else {
        fc_mapping <- aes(color = I(color_item), alpha = I(.data$.hilight))
    }

    p <- do.call(ggplot, c(list(g, layout = layout), dots))

    ## restore original category size
    if (
        length(node_label) > 1 &&
            getOption("cnetplot_subset", default = FALSE)
    ) {
        if (exprIsItemNum) {
            keep_idx <- match(names(originalItemNum), p$data$label)
            keep_idx <- keep_idx[!is.na(keep_idx)]
            p$data$size[keep_idx] <- originalItemNum[p$data$label[keep_idx]]
        }
    }

    if (is.na(color_edge)) {
        p <- p
    } else if (isTRUE(color_edge == "category")) {
        ed <- get_edge_data(g, names = TRUE)
        names(ed)[1] <- "category"
        p <- p +
            do.call(
                geom_edge,
                c(
                    list(
                        aes(color = .data$category),
                        data = ed
                    ),
                    edge_layer_params
                )
            )
    } else {
        edge_layer_params$color <- color_edge
        p <- p + do.call(geom_edge, edge_layer_params)
    }

    if (size_category > 0) {
        p <- p +
            geom_point(
                aes(size = .data$size, alpha = I(.data$.hilight)),
                data = td_filter(.data$.isCategory),
                color = color_category
            )
    }
    if (size_item > 0) {
        vals <- categorySizeValues
        if (
            exprIsItemNum &&
                length(node_label) > 1 &&
                getOption("cnetplot_subset", default = FALSE)
        ) {
            originalItemNum <- originalItemNum[names(x)]
            vals <- originalItemNum
        }
        vals <- vals[is.finite(vals)]
        size_breaks <- if (length(vals) > 0) {
            pretty(vals, n = min(4, length(unique(vals))))
        } else {
            NULL
        }
        if (missing(categorySizeBy) && exprIsItemNum) {
            size_name <- "size"
        } else {
            size_name <- rlang::as_label(categorySizeExpr)
        }
        p <- p +
            ggnewscale::new_scale_color() +
            geom_point(
                fc_mapping,
                data = td_filter(!.data$.isCategory),
                size = 3 * size_item
            ) +
            scale_size(
                name = size_name,
                range = c(3, 8) * size_category,
                breaks = size_breaks
            )
    }

    if (length(node_label) > 1 || node_label != "none") {
        label_params <- if (is.null(node_label_size)) {
            list()
        } else {
            list(size = node_label_size)
        }

        if (
            length(node_label) > 1 ||
                node_label %in% c("exclusive", "share")
        ) {
            p <- p + do.call(geom_cnet_label, c(list(node_label = "category"), label_params))
        }

        p <- p + do.call(geom_cnet_label, c(list(node_label = node_label), label_params))
    }

    return(p)
}

# maybe a geom_cnet_point function

#' @method ggplot_add cnet_label
#' @export
ggplot_add.cnet_label <- function(object, plot, object_name, ...) {
    params <- object$params
    node_label <- object$node_label
    default_params <- list(bg.color = "white", bg.r = .1)
    params <- modifyList(default_params, params)

    default_mapping <- aes(label = .data$label)
    if (is.null(object$mapping)) {
        params$mapping <- default_mapping
    } else {
        params$mapping <- modifyList(default_mapping, object$mapping)
    }

    x <- graph2list(plot$plot_env$data)

    if (!is.null(object$data)) {
        d <- object$data
    } else if (length(node_label) > 1) {
        d <- td_filter(.data$label %in% node_label)
    } else if (node_label == "none") {
        d <- td_filter(FALSE)
    } else if (node_label == "all") {
        d <- NULL
    } else if (node_label == "category") {
        d <- td_filter(.data$.isCategory)
    } else if (node_label == "item") {
        d <- td_filter(!.data$.isCategory)
    } else if (node_label %in% c("exclusive", "share")) {
        node_to_label <- lapply(seq_along(x), function(i) {
            j <- x[[i]] %in% unlist(x[-i])
            if (node_label == "exclusive") {
                return(x[[i]][!j])
            } else {
                return(x[[i]][j])
            }
        })
        d <- td_filter(.data$label %in% unlist(node_to_label))
    } else {
        d <- td_filter(!!str2lang(paste("foldChange", node_label)))
    }

    params$data <- d

    layer <- do.call(geom_text_repel, params) # (aes(label=.data$label), data = d, bg.color="white", bg.r=.1)
    ggplot_add(layer, plot, object_name, ...)
}

#' add labels of cnetplot
#'
#' @title geom_cnet_label
#' @param mapping aes mapping, default is NULL
#' @param data plot data, default is NULL
#' @param node_label which type of node label to be displayed, see also [cnetplot]
#' @param ... parameters that passed to `geom_text_repel`
#' @export
#' @author Guangchuang Yu
geom_cnet_label <- function(
  mapping = NULL,
  data = NULL,
  node_label = "all",
  ...
) {
    structure(
        list(
            mapping = mapping,
            data = data,
            node_label = node_label,
            params = list(...)
        ),
        class = "cnet_label"
    )
}


translate_cnet_category_size <- function(name) {
    if (!rlang::is_string(name)) {
        stop("`categorySize` must be a single column name, e.g. \"pvalue\".",
             call. = FALSE)
    }
    key <- tolower(trimws(name))
    if (key %in% c("genenum", "count", "itemnum")) {
        return(rlang::sym("itemNum"))
    }
    rlang::sym(name)
}

subset_cnet_list <- function(x, showCategory) {
    if (!is.numeric(showCategory)) {
        return(x[names(x) %in% showCategory])
    }

    n <- length(x)
    if (length(showCategory) == 1) {
        showCategory <- seq(showCategory)
    }

    if (any(showCategory) > n) {
        msg <- sprintf(
            "any showCategory value that is large than %d will be removed.",
            n
        )
        message(msg)
    }

    showCategory <- showCategory[showCategory <= n]
    return(x[showCategory])
}

subset_cnet_list_item <- function(x, showItem = "all") {
    if (length(showItem) == 1 && showItem == "all") {
        return(x)
    }

    lapply(x, \(y) y[y %in% showItem])
}
