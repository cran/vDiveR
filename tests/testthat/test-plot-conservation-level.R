test_that("conservation-level labels render as supported rich text", {
    plot <- plot_conservation_level(
        proteins_1host,
        conservation_label = 1,
        alpha = 0.8,
        base_size = 15
    )

    expect_error(ggplot2::ggplotGrob(plot), NA)
})

test_that("full labels retain absent conservation levels", {
    plot <- plot_conservation_level(
        proteins_1host,
        conservation_label = 1
    )
    labels <- plot$layers[[3]]$data$label

    expect_true(all(grepl("HD: 0 \\(0\\.0 %\\)", labels)))
    expect_true(all(grepl("ED: 0 \\(0\\.0 %\\)", labels)))
})

test_that("multi-host plots retain and render the shared legend", {
    plot <- expect_no_warning(
        plot_conservation_level(
            protein_2hosts,
            conservation_label = 1,
            host = 2
        )
    )

    grob <- expect_no_warning(ggplot2::ggplotGrob(plot))
    grob <- grid::grid.force(grob)
    grob_names <- grid::grid.ls(grob, print = FALSE, recursive = TRUE)$name

    expect_true(any(grepl("guide-box", grob_names)))
})
