#' Annotate a ggplot2 plot with text
#'
#' @description One can specify the relative position of the figure easily. Configurable margin, text and box justification. The added bonus in the following code is that you can specify which facet to annotate with something like facets=data.frame(cat1='blue', cat2='tall')
#'
#' From Stack Overflow comment https://stackoverflow.com/questions/22488563/ggplot2-annotate-layer-position-in-r
#'
#' @return Returns a geom layer 
#'
#' @author Rosen Matev
#'
#' @name annotate_textp
#'
#' @param label A string. The text to be plotted.
#' @param x Numeric. X axis position (min=0, max=1 - relative to the plot margins)
#' @param y Numeric. Y axis position (min=0, max=1 - relative to the plot margins)
#' @param hjust Numeric. The horizontal alignment  (0 = left, 0.5 = center, 1 = right). Default=0
#' @param size Numeric. Font size. Default={theme default}
#' @param alpha Numeric. Transparency (0:1). Default=NA
#'
#' @examples
#' p1 = ggplot2::ggplot(data.frame(x=1:10,y=1:10), ggplot2::aes(x,y)) + ggplot2::geom_point()
#' p1 + annotate_textp(label="x=1, y=0, hjust=1", 
#'                     x=1, y=0, hjust=1)
#' p1 + annotate_textp(label="x=0.1, y=0.9, hjust=0", 
#'                     x=0, y=0.9, hjust=0)
#' p1 + annotate_textp(label="x = 0.5, y=0.5, hjust=0.5, box_just=c(0.5,0.5)", 
#'                     x=0.5, y=0.5, hjust=0.5, box_just=c(0.5,0.5))
#' p1 + annotate_textp(label="x = 0.5, y=0.5, hjust=0.5, box_just=c(0.5,0.5), size=6, alpha=0.5", 
#'                     x=0.5, y=0.5, hjust=0.5, box_just=c(0.5,0.5), size=6, alpha=0.5)
#'
#' @export
#'

annotate_textp <- function(label, x, y, facets=NULL, hjust=0, vjust=0, color='black', alpha=NA,
                          family=thm$text$family, size=thm$text$size, fontface=1, lineheight=1.0,
                          box_just=ifelse(c(x,y)<0.5,0,1), margin=grid::unit(size/2, 'pt'), thm=ggplot2::theme_get()) {
  x = scales::squish_infinite(x)
  y = scales::squish_infinite(y)
  data = if (is.null(facets)) data.frame(x=NA) else data.frame(x=NA, facets)

  tg = grid::textGrob(
    label, x=0, y=0, hjust=hjust, vjust=vjust,
    gp=grid::gpar(col=scales::alpha(color, alpha), fontsize=size, fontfamily=family, fontface=fontface, lineheight=lineheight)
  )
  ts = grid::unit.c(grid::grobWidth(tg), grid::grobHeight(tg))
  vp = grid::viewport(x=x, y=y, width=ts[1], height=ts[2], just=box_just)
  tg = grid::editGrob(tg, x=ts[1]*hjust, y=ts[2]*vjust, vp=vp)
  inner = grid::grobTree(tg, vp=grid::viewport(width=grid::unit(1, 'npc')-margin*2, height=grid::unit(1, 'npc')-margin*2))

  ggplot2::layer(
    data = NULL,
    stat = ggplot2::StatIdentity,
    position = ggplot2::PositionIdentity,
    geom = ggplot2::GeomCustomAnn,
    inherit.aes = TRUE,
    params = list(
      grob=grid::grobTree(inner), 
      xmin=-Inf, 
      xmax=Inf, 
      ymin=-Inf, 
      ymax=Inf
    )
  )
}
