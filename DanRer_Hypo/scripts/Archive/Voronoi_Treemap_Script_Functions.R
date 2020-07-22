#  Copyright (C) 2012 Paul Murrell
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/licenses/gpl.txt

# Extracting coordinates from "gpc.poly"

# It is possible for us to end up with a cell containing a hole
# (because when the cell boundary is either zero area or
#  contains a zero-area indent, gpclib::intersect() can
#  produce isolated islands)

# To handle these cases, we just ignore the holes (which should
# be zero-area anyway) and take the outer boundary.

# There is defensive code in there to warn if the unexpected
# happens and there is more than one outer boundary.
getpts <- function(x) {
  pts <- get.pts(x)
  nonholes <- sapply(pts, function(y) !y$hole)
  if (sum(nonholes) < 1 || sum(nonholes) > 1) {
    warning("Cell or region with more than one boundary")
  }
  border <- which(nonholes)[1]
  list(x=pts[[border]]$x,
       y=pts[[border]]$y)
}


#  Copyright (C) 2012 Paul Murrell
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/licenses/gpl.txt

# The sites and weights for the diagram must be in "sites.cin"

recordSites <- function(sites, weights) {
  write.table(cbind(sites$x, sites$y, weights), "sites.cin",
              row.names=FALSE, col.names=FALSE)
}

# Read the sites and weights from "sites.cin"
# Generate an additively weighted Voronoi diagram
# Write bounding cells to "diagram.txt"
awv <- function(s, w, region, debug=FALSE, debugCell=FALSE) {
  recordSites(s, w)
  result <- system("./voronoiDiagram > diagram.txt")
  if (result)
    stop(paste("Voronoi diagram failed with error", result))
  roughCells <- readCells(s)
  tolerance <- .0015
  # tolerance <- max(diff(range(s$x)), diff(range(s$y)))*.000001
  tidyCells <- tidyCells(roughCells, tolerance, debug, debugCell)
  trimCells(tidyCells, region)
}

# Read the bounding cells from "diagram.txt"
# Return a list of cell information

# OLD version which assumed that order of cells matched order of sites
OLDreadCells <- function() {
  diagInfo <- readLines("diagram.txt")
  starts <- grep("^Vertex", diagInfo)
  lengths <- diff(c(starts, length(diagInfo)+1))
  
  readCell <- function(start, length) {
    vline <- readLines("diagram.txt", n=start)[start]
    list(vertex=as.numeric(strsplit(vline, " ")[[1]][2:3]),
         border=read.table("diagram.txt", skip=start, nrows=length-1,
                           na.strings="nan"))
  }
  
  mapply(readCell, starts, lengths, SIMPLIFY=FALSE)
}

# NOTE that the result may NOT be in same order as original sites
# SO have to check order in here
readCells <- function(s) {
  diagInfo <- readLines("diagram.txt")
  starts <- grep("^Vertex", diagInfo)
  lengths <- diff(c(starts, length(diagInfo)+1))
  vertices <- read.table(textConnection(diagInfo[starts]))[-1]
  
  readCell <- function(start, length) {
    vline <- readLines("diagram.txt", n=start)[start]
    if (length > 1) {
      border <- read.table("diagram.txt", skip=start, nrows=length-1,
                           na.strings="nan")
    } else {
      # There has been a case where the vertex is in the file
      # but there are no edges for the site
      # Handle this as an empty cell to keep the algorithm going
      border <- NULL
    }
    list(vertex=as.numeric(strsplit(vline, " ")[[1]][2:3]),
         border=border)
  }
  
  result <- vector("list", length(s$x))
  for (i in 1:length(s$x)) {
    sx <- s$x[i]
    sy <- s$y[i]
    eps <- .01
    match <- which((abs(sx - vertices[, 1]) < eps) &
                     (abs(sy - vertices[, 2]) < eps))
    if (length(match)) {
      result[[i]] <- readCell(starts[match], lengths[match])
    } else {
      # Return empty cell to allow algorithm to continue
      result[[i]] <- list(vertex=c(sx, sy),
                          border=NULL)
      # stop("Site not present in voronoiDiagram result")
    }
  }
  result
}

# Tidy the cell information
# Return a list of lists with x/y
tidyCells <- function(cells, tolerance, debug=FALSE, debugCell=FALSE) {
  lapply(cells, tidyCell, tolerance, debug, debugCell)
}

debugCell <- FALSE

getCellBorder <- function(border, seenBefore=numeric(), tolerance,
                          debug, debugCell) {
  N <- nrow(border)
  direction <- 1
  vx <- numeric(N)
  vy <- numeric(N)
  vcount <- 1
  complete <- FALSE
  endCondition <- FALSE
  returning <- length(seenBefore)
  visited <- rep(FALSE, N)
  # Two start cases:
  # Segment that starts on boundary
  # No segments that start on boundary (just start with first segment)
  start <- which(border[, 1] == -2*scale |
                   border[, 1] == 2*scale |
                   border[, 2] == -2*scale |
                   border[, 2] == 2*scale)
  # If this is a return visit, it must be a start point
  # we have not seen before
  start <- start[!start %in% seenBefore]
  if (length(start)) {
    index <- start[1]
    startingX <- border[index, 1]
    startingY <- border[index, 2]
    seenBefore <- c(seenBefore, index)
  } else {
    # Look for segment that *ends* on boundary
    start <- which(border[, 3] == -2*scale |
                     border[, 3] == 2*scale |
                     border[, 4] == -2*scale |
                     border[, 4] == 2*scale)
    start <- start[!start %in% seenBefore]
    if (length(start)) {
      index <- start[1]
      direction <- -1
      startingX <- border[index, 3]
      startingY <- border[index, 4]
      seenBefore <- c(seenBefore, index)
    } else if (returning) {
      # If we're returning for another look
      # (means we ended on the boundary last time)
      # AND we have not found another boundary to start from
      # break straight back out
      complete <- TRUE
      vx <- numeric()
      vy <- numeric()
      endCondition <- "breakout"
    } else {
      # Start at first point
      index <- 1
      startingX <- border[index, 1]
      startingY <- border[index, 2]
    }
  }
  # Following segments:
  # First check whether next segment follows on from current one
  # If not, search all other segments for segment that follows on
  # If no match, find segment that *ends* at current end and
  #   follow that segment *backwards*
  while (!complete) {
    if (debug && debugCell)
      cat(paste(index, direction, "\n"))
    if (direction > 0) {
      x1 <- 1
      y1 <- 2
      x2 <- 3
      y2 <- 4
    } else {
      x1 <- 3
      y1 <- 4
      x2 <- 1
      y2 <- 2
    }
    vx[vcount] <- border[index, x1]
    vy[vcount] <- border[index, y1]
    vcount <- vcount + 1
    visited[index] <- TRUE
    if (all(visited)) {
      complete <- TRUE
    } else {
      # Try next segment
      nextIndex <- index + direction
      endX <- border[index, x2]
      endY <- border[index, y2]
      startX <- border[nextIndex, x1]
      startY <- border[nextIndex, y1]
      if (nextIndex > 0 && nextIndex <= N &&
          sim(endX, startX, tolerance) && sim(endY, startY, tolerance)) {
        index <- nextIndex
      } else {
        # Try segment that starts where we ended
        newIndex <- startsAt(border, endX, endY, direction, tolerance)
        # Do NOT use visited segment
        newIndex <- newIndex[!visited[newIndex]]
        if (length(newIndex)) {
          index <- newIndex[1]
        } else {
          # Try segment that ends where we ended
          direction <- -direction
          newIndex <- startsAt(border, endX, endY, direction,
                               tolerance)
          # Do NOT use visited segment
          newIndex <- newIndex[!visited[newIndex]]
          if (length(newIndex)) {
            index <- newIndex[1]
          } else {
            # Run out of segments
            # cat("Found no next edge\n")
            complete <- TRUE
            warning("Failed to trace cell")
            # recover()
          }
        }
      }
    }
    # Stopping conditions:
    # Hit boundary
    # Got back to start
    endCondition <- stopping(border[index, x2], border[index, y2],
                             startingX, startingY, tolerance)
    if (endCondition == "boundary") {
      # Add boundary nodes to border
      vx[vcount] <- border[index, x1]
      vx[vcount+1] <- border[index, x2]
      vy[vcount] <- border[index, y1]
      vy[vcount+1] <- border[index, y2]
      vcount <- vcount + 2
      seenBefore <- c(seenBefore, index)
      complete <- TRUE
    } else if (endCondition == "start") {
      complete <- TRUE
    }
  }
  list(x=vx[1:(vcount-1)], y=vy[1:(vcount-1)],
       end=endCondition, seen=seenBefore)
}

tidyCell <- function(cell, tolerance, debug=FALSE, debugCell=FALSE) {
  if (debug && debugCell)
    print(cell$vertex)
  border <- cell$border
  # Handle empty cells
  if (is.null(border))
    return(NULL)
  ok <- !apply(cell$border, 1, function(x) any(is.na(x)))
  border <- border[ok, ]
  
  result <- getCellBorder(border, tolerance=tolerance,
                          debug=debug, debugCell=debugCell)
  if (result$end == "boundary") {
    # Check for other boundaries
    results <- list(result)
    repeat {
      result <- getCellBorder(border, seenBefore=result$seen,
                              tolerance=tolerance,
                              debug=debug, debugCell=debugCell)
      if (result$end == "breakout")
        break
      results <- c(results, list(result))
    }
    if (length(results) == 1) {
      result <- results[[1]]
    } else {
      result <- results
      class(result) <- "multipleBorders"
    }
  }
  
  if (inherits(result, "multipleBorders") ||
      result$end == "boundary") {
    # Need to close the cell
    cellBorder <- closeCell(result, cell$vertex, tolerance)
  } else { # if (end == "start") {
    # This covers two cases:
    # end == "start" means that we have completed a loop
    # end == FALSE means that we ran out of points
    # The latter case may produce some unpredictable results
    cellBorder <- list(x=c(result$x, result$x[1]),
                       y=c(result$y, result$y[1]))
  }
  cellBorder
}

# SIDES

# 1 = left
# 2 = top
# 3 = right
# 4 = bottom

# CORNERS

# 1 = top-left
# 2 = top-right
# 3 = bottom-right
# 4 = bottom-left

side <- function(x, y) {
  if (x == -2*scale)
    return(1)
  if (y == 2*scale)
    return(2)
  if (x == 2*scale)
    return(3)
  if (y == -2*scale)
    return(4)
}

clockCorner <- function(side) {
  switch(side, 1, 2, 3, 4)
}

clockSide <- function(corner) {
  switch(corner, 2, 3, 4, 1)
}

antiCorner <- function(side) {
  switch(side, 4, 1, 2, 3)
}

antiSide <- function(corner) {
  switch(corner, 1, 2, 3, 4)
}

cornerX <- c(-2*scale, 2*scale, 2*scale, -2*scale)
cornerY <- c(2*scale, 2*scale, -2*scale, -2*scale)

closeClock <- function(x, y, start, end) {
  side <- end
  repeat {
    corner <- clockCorner(side)
    x <- c(x, cornerX[corner])
    y <- c(y, cornerY[corner])
    side <- clockSide(corner)
    if (side == start) {
      break
    }
  }
  x <- c(x, x[1])
  y <- c(y, y[1])
  list(x=x, y=y)
}

closeAnti <- function(x, y, start, end) {
  side <- end
  repeat {
    corner <- antiCorner(side)
    x <- c(x, cornerX[corner])
    y <- c(y, cornerY[corner])
    side <- antiSide(corner)
    if (side == start) {
      break
    }
  }
  x <- c(x, x[1])
  y <- c(y, y[1])
  list(x=x, y=y)
}

closeCell <- function(cell, vertex, tol) {
  UseMethod("closeCell")
}

stretchX <- function(x, y, N, side) {
  switch(side,
         c(x[1], max(x), x[N]),
         c(max(-2*scale, x[1] - scale*.05),
           x[which.min(y)],
           min(2*scale, x[N] + scale*.05)),
         c(x[1], min(x), x[N]),
         c(max(-2*scale, x[1] - scale*.05),
           x[which.max(y)],
           min(2*scale, x[N] + scale*.05)))           
}

stretchY <- function(x, y, N, side) {
  switch(side,
         c(max(-2*scale, y[1] - scale*.05),
           y[which.max(x)],
           min(2*scale, y[N] + scale*.05)),
         c(y[1], max(y), y[N]),
         c(max(-2*scale, y[1] - scale*.05),
           y[which.min(x)],
           min(2*scale, y[N] + scale*.05)),
         c(y[1], min(y), y[N]))
}

library(sp)
closeCell.default <- function(cell, vertex, tol) {
  # Two options:  go round clip region boundary clockwise or anit-clockwise
  # Try first one, check if vertex is "inside" the result
  # If not, do second one (and check that vertex is "inside" that result!)
  
  x <- cell$x
  y <- cell$y
  
  # It is possible to get here with a cell that STARTS on a
  # boundary but does not end on a boundary (because voronoiDiagram
  # is capable of producing this sort of thing sometimes;  I have seen it!)
  # This case is characterised by cell$end being FALSE
  # In that case, just hail mary and join start to end
  # (so that end is also on boundary)
  if (is.logical(cell$end)) {
    if (!cell$end) {
      x <- c(x, x[1])
      y <- c(y, y[1])
    }
  }
  
  # ASSUME that both first and last vertices are on boundary!
  N <- length(x)
  startSide <- side(x[1], y[1])
  endSide <- side(x[N], y[N])
  
  # Start and end on same side
  if (startSide == endSide) {
    # Special case:  startPOINT and endPOINT same
    if (sim(x[1], x[N], tol) &&
        sim(y[1], y[N], tol)) {
      newx <- stretchX(x, y, N, startSide)
      newy <- stretchY(x, y, N, startSide)
      x <- newx
      y <- newy
    }
    # Just join start to end
    cell <- list(x=c(x, x[1]), y=c(y, y[1]))
    if (point.in.polygon(vertex[1], vertex[2],
                         cell$x, cell$y) == 0) {
      boundRect <- suppressWarnings(as(list(x=c(-2*scale, -2*scale, 
                                                2*scale, 2*scale),
                                            y=c(-2*scale, 2*scale, 
                                                2*scale, -2*scale)),
                                       "gpc.poly"))
      # "Subtract" smallCell from bound rect to get largeCell
      cellPoly <- suppressWarnings(as(cell, "gpc.poly"))
      cellPoly <- intersect(append.poly(cellPoly, boundRect), boundRect)
      pts <- getpts(cellPoly)
      cell <- list(x=c(pts$x, pts$x[1]), y=c(pts$y, pts$y[1]))
      if (point.in.polygon(vertex[1], vertex[2],
                           cell$x, cell$y) == 0) {
        stop("Failed to close cell")
      }
    }
    
  } else {
    cell <- closeClock(x, y, startSide, endSide)
    if (point.in.polygon(vertex[1], vertex[2],
                         cell$x, cell$y) == 0) {
      cell <- closeAnti(x, y, startSide, endSide)
      if (point.in.polygon(vertex[1], vertex[2],
                           cell$x, cell$y) == 0) {
        stop("Failed to close cell")
      }
    }
  }
  cell
}

closeCell.multipleBorders <- function(cell, vertex, tol) {
  result <- lapply(cell, closeCell, vertex, tol)
  class(result) <- "multipleCells"
  result
}

# Does this constant need adjusting if 'scale' in init.R
# is altered??
# The .0015 is fairly sensitive
# Examples seen where .001 is (just) too small
# Examples seen where .005 is too big
sim <- function(a, b, tol=.0015) {
  abs(a - b) < tol
}

stopping <- function(endX, endY, startX, startY, tol) {
  if (endX == -2*scale || endX == 2*scale ||
      endY == -2*scale || endY == 2*scale) {
    # cat("Hit boundary\n")
    return("boundary")
  }
  if (sim(endX, startX, tol) && sim(endY, startY, tol)) {
    # cat("Back to start\n")
    return("start")
  }
  return(FALSE)
}

startsAt <- function(border, x, y, direction, tol) {
  if (direction > 0) {
    x1 <- 1
    y1 <- 2
  } else {
    x1 <- 3
    y1 <- 4
  }
  which(sim(border[, x1], x, tol) &
          sim(border[, y1], y, tol))
}

# Take polygons from Voronoi cells and intersect them with
# outer polygon (e.g., circle radius 1000)
# Return list of "gpc.poly"
convertCell <- function(cell) {
  UseMethod("convertCell")
}

convertCell.default <- function(cell) {
  # Handle empty cells
  if (is.null(cell)) {
    new("gpc.poly")
  } else {
    suppressWarnings(as(cell, "gpc.poly"))
  }
}

convertCell.multipleCells <- function(cell) {
  polys <- suppressWarnings(lapply(cell, as, "gpc.poly"))
  Reduce(intersect, polys)
}

trimCells <- function(cells, region) {
  polys <- lapply(cells, convertCell)
  lapply(polys, intersect, region)
}


#  Copyright (C) 2012 Paul Murrell
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/licenses/gpl.txt

library(gpclib) # Non-commercial only!

debugPlotGen <- function() {
  plotStarted <- FALSE
  x <- 1
  prevYoverall <- NA
  prevYmax <- NA
  
  function(yoverall, ymax) {
    if (!plotStarted) {
      plot.new()
      plot.window(c(0, 200), c(0, .2))
      abline(h=c(0, .001), col="grey", lty=c("solid", "dashed"))
      axis(1)
      axis(2)
      box()
      plotStarted <<- TRUE
    }
    segments(x, prevYoverall, x + 1, yoverall, col="grey")
    segments(x, prevYmax, x + 1, ymax, col="black")
    x <<- x + 1
    prevYoverall <<- yoverall
    prevYmax <<- ymax
  }
}

debugPlot <- debugPlotGen()

prevError <- 0

cellError <- function(a, target) {
  normA <- a/sum(a)
  diff <- abs(normA - target)
  max(diff)
}

breaking <- function(a, target, max=TRUE, debug=FALSE, dplot=FALSE) {
  if (max) {
    # Stop when largest individual cell error is less than .1%
    # (the default)
    err <- cellError(a, target)
    if (debug)
      cat(paste("Difference: ", err,
                " (", abs(err - prevError), ")",
                "\n", sep=""))
    stopping <- err < .001
    prevError <<- err
    if (!debug && dplot) {
      debugPlot(err, err)
    }
  } else {
    normA <- a/sum(a)
    diff <- abs(normA - target)
    # Stop when *change* in *total* cell error is tiny
    # (i.e., we are not improving the solution)
    if (debug)
      cat(paste("Difference: ", round(sum(diff), 2),
                " (", round(abs(sum(diff) - prevError), 3), ")",
                "\n", sep=""))
    stopping <- abs(sum(diff) - prevError) < .0001
    prevError <<- sum(diff)
    if (!debug && dplot) {
      debugPlot(sum(diff)/sum(target), max(diff/target))
    }
  }
  stopping
}

OLDadjustWeights <- function(w, a, target) {
  # The choice of 'eps' is important
  # Either too large (e.g., 10) or too small (e.g., .001)
  # can result in failure to converge
  eps <- scale/10000
  # Watch out for zero-area cells (set to target value)
  # a <- ifelse(a == 0, target/sum(target)*sum(a), a)
  # Watch out for zero-area cells (set to small value)
  a <- ifelse(a == 0, .01*sum(a), a)
  normA <- a/sum(a)
  w <- ifelse(abs(w) < eps,
              ifelse(w < 0, -eps, eps),
              w)
  w + abs(w)*((target - normA)/target)
}

# Instead of adjusting by a multiple of current weight
# adjust by multiple of average absolute weights
# This avoids problem of getting stuck at a tiny weight
# (and stabilizes the algorithm generally)
adjustWeights <- function(w, a, target) {
  # Watch out for zero-area cells (set to small value)
  a <- ifelse(a == 0, .01*sum(a), a)
  normA <- a/sum(a)
  w + mean(abs(w))*((target - normA)/target)
}

# Use centroid finding code in 'soiltexture' for now
# May be able to find a better source later (?)
library(soiltexture)
shiftSites <- function(s, k) {
  newSites <- mapply(function(poly, x, y) {
    # Handle empty polygons
    if (length(poly@pts)) {
      pts <- getpts(poly)
      TT.polygon.centroids(pts$x,
                           pts$y)
    } else {
      c(x=x, y=y)
    }
  },
  k, as.list(s$x), as.list(s$y),
  SIMPLIFY=FALSE)
  list(x=sapply(newSites, "[", "x"),
       y=sapply(newSites, "[", "y"))
}

OLDshiftWeights <- function(s, w) {
  n <- length(s$x)
  factor <- Inf
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        # Deviation from published algorithm here
        # to use abs(w) so that ensure non-overlapping
        # circles even when weights are negative
        if (FALSE) {
          f = sqrt((s$x[i] - s$x[j])^2 +
                     (s$y[i] - s$y[j])^2)/(abs(w[i]) + abs(w[j]))
        }
        # Original published algorithm
        if (FALSE) {
          f = sqrt((s$x[i] - s$x[j])^2 +
                     (s$y[i] - s$y[j])^2)/(w[i] + w[j])
        }
        # Another variation where max(abs(w)) used
        # to avoid circle overlapping other site
        if (FALSE) {
          f = sqrt((s$x[i] - s$x[j])^2 +
                     (s$y[i] - s$y[j])^2)/
            max(abs(w[i]), abs(w[j]))
        }
        
        if (f > 0 && f < factor)
          factor <- f
      }
    }
  }
  if (factor < 1) {
    w <- w*factor
  }
  w
}

# Variation from published algorithm
# Allow factor adjustment per pair of sites
shiftWeights <- function(s, w) {
  n <- length(s$x)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        # Deviation from published algorithm here
        # to use abs(w) so that ensure non-overlapping
        # circles even when weights are negative
        f = sqrt((s$x[i] - s$x[j])^2 +
                   (s$y[i] - s$y[j])^2)/(abs(w[i]) + abs(w[j]))
        
        if (f > 0 && f < 1) {
          w[i] <- w[i]*f
          w[j] <- w[j]*f
        }
      }
    }
  }
  w
}

# The algorithm can fail to converge sometimes so
# just give up after 'maxIteration's
allocate <- function(names, s, w, outer, target, maxIteration=200,
                     debug=FALSE, dplot=FALSE, debugCell=FALSE) {
  count <- 1
  debugPlot <<- debugPlotGen()
  repeat {
    k <- awv(s, w, outer, debug, debugCell)
    areas <- lapply(k, area.poly)
    if (debug) {
      drawRegions(list(names=names, k=k,
                       s=s, w=w, a=areas, t=target),
                  debug)
      info <- rbind(area=round(unlist(areas)/sum(unlist(areas)), 4),
                    target=round(target, 4),
                    weight=round(w, 1))
      colnames(info) <- names
      print(info)
    }
    if (count > maxIteration ||
        breaking(unlist(areas), target, debug=debug, dplot=dplot)) {
      return(list(names=names, k=k, s=s, w=w, a=areas, t=target))
    } else {
      w <- adjustWeights(w, unlist(areas), target)
      s <- shiftSites(s, k)
      w <- shiftWeights(s, w)
    }
    count <- count + 1
  }
}


#  Copyright (C) 2012 Paul Murrell
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/licenses/gpl.txt

drawPoly <- function(p, name, gp=gpar()) {
  if (length(p@pts)) {
    pts <- getpts(p)
    grid.polygon(pts$x, pts$y, default="native",
                 gp=gp, name=name)
  }
}

library(grid)

polyRangeX <- function(p) {
  if (length(p@pts)) {
    pts <- getpts(p)
    range(pts$x)
  } else {
    NA
  }
}

polyRangeY <- function(p) {
  if (length(p@pts)) {
    pts <- getpts(p)
    range(pts$y)
  } else {
    NA
  }
}

drawRegions <- function(result, label=FALSE, top=TRUE,
                        gp=gpar(), newpage=TRUE, debug=FALSE) {
  names <- result$names
  k <- result$k
  sites <- result$s
  weights <- result$w
  if (newpage) {
    grid.newpage()
    xrange <- range(unlist(lapply(k, polyRangeX)), na.rm=TRUE)
    yrange <- range(unlist(lapply(k, polyRangeY)), na.rm=TRUE)
    pushViewport(viewport(width=.8, height=.8,
                          xscale=xrange, yscale=yrange))
  }
  if (top) {
    invisible(mapply(drawPoly, k, names, MoreArgs=list(gp=gpar(lwd=8))))
    invisible(mapply(drawPoly, k, names,
                     MoreArgs=list(gp=gpar(lwd=6, col="grey90"))))
  } else {
    invisible(mapply(drawPoly, k, names, MoreArgs=list(gp=gp)))
  }
  if (label) {
    if (top) {
      cex=1
    } else {
      cex=.5
    }
    grid.text(result$names, sites$x, sites$y, default="native",
              gp=gpar(cex=cex))
    if (debug) {
      col <- ifelse(weights < 0, "red", "black")
      grid.circle(sites$x, sites$y, default="native",
                  r=unit(weights, "native"),
                  gp=gpar(col=col, fill=NA))
    }
  }
}


#  Copyright (C) 2012 Paul Murrell
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.gnu.org/licenses/gpl.txt

# Code to help debug steps in the Voronoi Treemap algorithm

# Plot a tidied cell
# Use semitrans blue so can see backtracks
# Draw nice big blob at start
# Draw slightly smaller blob at end (so can see if got back to start)

drawTidyCell <- function(...) {
  UseMethod("drawTidyCell")
}

drawTidyCell.default <- function(cell, 
                                 lineCol=rgb(0,0,1,.5), 
                                 lineWidth=1,
                                 startFill=rgb(0,0,1,.5),
                                 endFill=rgb(1,0,0,.5),
                                 newpage=TRUE) {
  if (is.null(cell))
    return()
  if (newpage) {
    grid.newpage()
    pushViewport(viewport(width=.9, height=.9,
                          xscale=c(-2*scale, 2*scale),
                          yscale=c(-2*scale, 2*scale)))
  }
  grid.circle(cell$x[1], cell$y[1], r=unit(2, "mm"),
              default.units="native",
              gp=gpar(col=NA, fill=startFill))
  grid.lines(cell$x, cell$y,  
             default.units="native",
             gp=gpar(col=lineCol, lwd=lineWidth))
  N <- length(cell$x)
  grid.circle(cell$x[N], cell$y[N], r=unit(1, "mm"),
              default.units="native",
              gp=gpar(col=NA, fill=endFill))
}

drawTidyCell.multipleCells <- function(cell, 
                                       lineCol=rgb(0,0,1,.5), 
                                       lineWidth=1,
                                       startFill=rgb(0,0,1,.5),
                                       endFill=rgb(1,0,0,.5),
                                       newpage=TRUE) {
  drawTidyCell(cell[[1]],
               lineCol, lineWidth, startFill, endFill, newpage)
  for (i in 2:length(cell)) {
    drawTidyCell(cell[[i]],
                 lineCol, lineWidth, startFill, endFill, newpage=FALSE)
  }
}

drawRoughCell <- function(cell,
                          lineCol=rgb(0,0,1,.5), 
                          lineWidth=1,
                          startFill=rgb(0,0,1,.5),
                          endFill=rgb(1,0,0,.5),
                          newpage=TRUE) {
  if (is.null(cell$border))
    return()
  if (newpage) {
    grid.newpage()
    pushViewport(viewport(width=.9, height=.9,
                          xscale=c(-2*scale, 2*scale),
                          yscale=c(-2*scale, 2*scale)))
  }
  grid.segments(cell$border[, 1], cell$border[, 2],
                cell$border[, 3], cell$border[, 4],
                default.units="native",
                gp=gpar(col=lineCol, lwd=lineWidth))
}

drawTidyCells <- function(cells) {
  for (i in 1:length(cells)) {
    grid.newpage()
    pushViewport(viewport(width=.9, height=.9))
    pushViewport(viewport(xscale=c(-2*scale, 2*scale),
                          yscale=c(-2*scale, 2*scale)))
    grid.rect(gp=gpar(col=NA, fill="grey80"))
    for (j in 1:length(cells)) {
      drawTidyCell(cells[[j]], 
                   lineCol="white",
                   lineWidth=3,
                   startFill=NA, endFill=NA, 
                   newpage=FALSE)
    }
    drawTidyCell(cells[[i]], newpage=FALSE)
  }
}

drawRoughCells <- function(cells) {
  for (i in 1:length(cells)) {
    grid.newpage()
    pushViewport(viewport(width=.9, height=.9))
    pushViewport(viewport(xscale=c(-2*scale, 2*scale),
                          yscale=c(-2*scale, 2*scale)))
    grid.rect(gp=gpar(col=NA, fill="grey80"))
    for (j in 1:length(cells)) {
      drawRoughCell(cells[[j]], 
                    lineCol="white",
                    lineWidth=3,
                    startFill=NA, endFill=NA, 
                    newpage=FALSE)
    }
    drawRoughCell(cells[[i]], newpage=FALSE)
  }
}

drawTidyCellsOnePage <- function(cells) {
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(3, 4)))
  for (i in 1:length(cells)) {
    pushViewport(viewport(layout.pos.col=(i-1)%%4 + 1,
                          layout.pos.row=(i-1)%/%4 + 1))
    grid.rect()
    drawTidyCell(cells[[i]], newpage=FALSE)
    popViewport()
  }
  popViewport()
}

drawSites <- function(s, w, scol="grey", r=unit(.5, "mm"),
                      weights=TRUE, newpage=TRUE) {
  if (newpage) {
    grid.newpage()
    pushViewport(viewport(xscale=c(-1*scale, 1*scale),
                          yscale=c(-1*scale, 1*scale)))
  }
  grid.circle(s$x, s$y, r=r, 
              gp=gpar(col=NA, fill=scol),
              default.units="native")
  if (weights) {
    col <- ifelse(w < 0, "red", "black")
    grid.circle(s$x, s$y, r=abs(w), 
                gp=gpar(col=col),
                default.units="native")
  }
}

drawRegion <- function(region, gp=gpar(), newpage=TRUE) {
  if (newpage) {
    grid.newpage()
    pushViewport(viewport(xscale=c(-2*scale, 2*scale),
                          yscale=c(-2*scale, 2*scale)))
  }
  pts <- getpts(region)
  grid.polygon(pts$x,
               pts$y,
               gp=gp,
               default.units="native")
}

