draw_key_roundedrect <- function(data, params, size) {
  lwd <- min(data$size, min(size)/4)
  grid::roundrectGrob(width = unit(1, 'npc') - unit(lwd, 'mm'),
                      height = unit(1, 'npc') - unit(lwd, 'mm'),
                      r = unit(data$roundness, 'mm'),
                      just = 'centre',
                      gp = gpar(col = data$colour,
                                fill = alpha(data$fill, data$alpha),
                                lty = data$linetype,
                                lwd = lwd * .pt))
}


#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
#' @include geom-rect.r
GeomRoundedTile <- ggproto("GeomRoundedTile", GeomRect,
  extra_params = c("na.rm"),

  setup_data = function(data, params) {
    data$width <- data$width %||% params$width %||% resolution(data$x, FALSE)
    data$height <- data$height %||% params$height %||% resolution(data$y, FALSE)

    epsilon = .7

    transform(data,
      xmin = x - ((width * epsilon) / 2),
      xmax = x + ((width * epsilon) / 2),  width = NULL,
      ymin = y - ((height * epsilon) / 2),
      ymax = y + ((height * epsilon) / 2), height = NULL
    )
  },

  default_aes = aes(fill = 'grey20', colour = NA, size = 0.0, linetype = 1,
                    alpha = NA, width = NA, height = NA, roundness = 0),

  required_aes = c('x', 'y'),

  draw_key = draw_key_roundedrect
)



geom_rounded_tile <- function(mapping = NULL, data = NULL,
															stat = 'identity', position = 'identity',
															...,
															na.rm = FALSE,
															show.legend = NA,
															inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRoundedTile,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

if (F) {
rounded_square <- function(x, y, size) {
  roundrectGrob(x = x, y = y, width = size, height = size,
                default.units = 'native',
                r = unit(0.0, 'snpc'),
                just = 'centre')
}

GeomRoundedSquare <- ggproto('GeomRoundedSquare', Geom,
  draw_panel = function(data, panel_scales, coord, img, na.rm = FALSE) {
    coords <- coord$transform(data, panel_scales)
    ggplot2:::ggname('geom_rounded_square',
                     rounded_square(coords$x, coords$y, coords$size))
  }, non_missing_aes = c('size'), required_aes = c('x', 'y'),
  default_aes = aes(size = 0.05),
  icon = function(.) {
  }, desc_params = list(), seealso = list(geom_point = GeomPoint$desc),
  examples = function(.) { })

geom_rounded_square <- function(mapping = NULL, data = NULL, stat = "identity",
                                position = "identity", na.rm = FALSE,
                                show.legend = NA, inherit.aes = TRUE, ...) {
  layer(data = data, mapping = mapping, stat = stat,
        geom = GeomRoundedSquare,
        position = position, show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, img = img, ...))
}
}
