#' Modified theme_minimal for ggplots
#'
#' @description From Albert Rapp's blog - nice adjustments to default theme minimal
#'
#' https://github.com/AlbertRapp/3mw_newsletters/blob/master/37_Jul_19_2023/37_Jul_19_2023.qmd?utm_source=alberts-newsletter.beehiiv.com&utm_medium=newsletter&utm_campaign=3mw-create-your-own-ggplot-theme
#'
#' @return Returns a theme
#'
#' @author Albert Rapp
#'
#' @name theme_minimal_modified
#'
#' @examples
#' penguins = palmerpenguins::penguins |> 
#'   dplyr::filter(!is.na(sex))
#' 
#' basic_plot = penguins |> 
#'   ggplot2::ggplot(ggplot2::aes(x = bill_length_mm, y = body_mass_g, fill = species)) +
#'   ggplot2::geom_point(shape = 21, size = 5, alpha = 0.85, color = 'grey10') +
#'   ggplot2::labs(
#'     x = 'Bill length (in mm)',
#'     y = ggplot2::element_blank(),
#'     fill = 'Species',
#'     title = 'Penguins from the Palmer Archipelago',
#'     subtitle = 'Penguin weights (in g)',
#'     caption = 'Data: {palmerpenguins} R package'
#'   )
#'
#' basic_plot +
#'   ggplot2::theme_minimal(base_size = 16, base_family = 'Source Sans Pro') +
#'   ggplot2::scale_fill_manual(values = c('Adelie'= '#da2c38', 'Chinstrap'= '#FED18C', 'Gentoo'= '#30C5FF')) +
#'   theme_minimal_modified()
#'
#' @export
#'

theme_minimal_modified = function() {
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      color = 'grey90', linetype = 2
    ),
    plot.title.position = 'plot',
    plot.caption.position = 'plot',
    plot.title = element_text(
      family = 'Cambria', 
      size = rel(1.7),
      margin = margin(b = 7, unit = 'mm')
    ),
    plot.subtitle = element_text(size = rel(1.1)),
    text = element_text(color = 'grey20'),
    axis.text = element_text(color = 'grey30'),
    panel.background = element_rect(color = 'grey90')
  )
}

