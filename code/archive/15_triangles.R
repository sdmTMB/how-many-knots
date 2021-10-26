#' Convert inla.mesh to sp objects
#'
#' @param mesh An \code{\link{inla.mesh}} object
#' @return A list with \code{sp} objects for triangles and vertices:
#' \describe{
#' \item{triangles}{\code{SpatialPolygonsDataFrame} object with thetriangles in
#' the same order as in the original mesh, but each triangle looping through
#' the vertices in clockwise order (\code{sp} standard) instead of
#' counterclockwise order (\code{inla.mesh} standard). The \code{data.frame}
#' contains the vertex indices for each triangle, which is needed to link to
#' functions defined on the vertices of the triangulation.
#' \item{vertices}{\code{SpatialPoints} object with the vertex coordinates,
#' in the same order as in the original mesh.}
#' }
#' @export
inla.mesh2sp <- function(mesh) {
  crs <- inla.CRS(inla.CRSargs(mesh$crs))
  isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp' doesn't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before
calling inla.mesh2sp()."))
  }
  
  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(lapply(
      1:nrow(mesh$graph$tv),
      function(x) {
        tv <- mesh$graph$tv[x, , drop = TRUE]
        Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)],
                                       1:2,
                                       drop = FALSE])),
                 ID = x)
      }
    ),
    proj4string = crs
    ),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)
  
  list(triangles = triangles, vertices = vertices)
}


#' Overlay inla.mesh with points to identify points / triangle
#'
#' @param inla_mesh An \code{\link{inla.mesh}} object
#' @param coordinates A matrix or dataframe of coordinates
overlay = function(inla_mesh, coordinates) {
# convert mesh to list of triangles
# see https://groups.google.com/g/r-inla-discussion-group/c/z1n1exlZrKM
sp = inla.mesh2sp(inla_mesh)
# add ids
sp[[1]]@data$id = 1:nrow(sp[[1]]@data)
# turn points to SP data frame and add the ids
p <- SpatialPointsDataFrame(coords=coordinates, 
                            data.frame(ids=1:nrow(coordinates)))
#plot(sp[[1]])
#points(p,col="red")
res <- sp::over(p, sp[[1]])
return(res)
}
