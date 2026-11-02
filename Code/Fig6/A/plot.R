plot_files <- c(plot_files, local({
  input <- prepare_plot_data("A")
  plots <- lapply(seq_along(input$cells), function(i) {
    make_tb_panel(input$cells[i], input$plot_data, input$statistics, show_y = i == 1L)
  })
  save_plot(combine_box_plots(plots, "A"), "A")
}))
