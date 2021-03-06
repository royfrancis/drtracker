% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drtracker.R
\name{ltrack}
\alias{ltrack}
\title{Track larval movement in 24 or 48 well assay plates.}
\usage{
ltrack(files = NULL, wells = 24, markededges = TRUE, fps = 25,
  mm = 5.4, framelinking = TRUE, follow = "none", trackfeatures = TRUE,
  activitydist = 5, spotcoverage = TRUE, spotalphahull = TRUE,
  alphavalue = 4, spotdensity = TRUE, spotmsd = FALSE, spotratio = TRUE,
  outercircleradius = 8.1, innercircleradius = 4.8,
  filenamediscard = ".txt", exportdata = TRUE, exportplot = TRUE,
  useexisting = FALSE)
}
\arguments{
\item{files}{A character or vector of paths or filenames. See details.}

\item{wells}{A numeric indicating plate format. Either 24 or 48.}

\item{markededges}{A logical indicating if corners of plates have been marked. See details.}

\item{fps}{A numeric indicating framerate of video in number of frames per second.}

\item{mm}{A numeric indicating number of pixels in 1 mm.}

\item{framelinking}{A logical indicating if spots should be linked frame to frame.}

\item{follow}{A character indicating if the algorithm must be followed at every step. Options are 'none', 'interactive' or 'track'. See details.}

\item{trackfeatures}{A logical indicating if track features should be computed.}

\item{activitydist}{A numeric indicating distance in mm after which the larvae is considered active.}

\item{spotcoverage}{A logical indicating if the coverge must be computed. See details.}

\item{spotalphahull}{A logical indicating if alphahull of point cloud is to be calculated. See details.}

\item{alphavalue}{A numeric indicating alpha value for alphahull. See details.}

\item{spotdensity}{A logical indicating if 2D kernal density is to be estimated. See details.}

\item{spotmsd}{A logical indicating if the minimum spanning distance should be computed. See details.}

\item{spotratio}{A logical indicating if spotratio should be computed. See details.}

\item{outercircleradius}{A numeric indicating radius of outer circle (radius of one well) in mm.}

\item{innercircleradius}{A numeric indicating radius of inner circle in mm.}

\item{filenamediscard}{A character for part of the filename to be removed.}

\item{exportdata}{A logical indicating if data tables must be exported to working directory.}

\item{exportplot}{A logical if results must be plotted and exported to working directory.}

\item{useexisting}{A logical indicating if existing files in working directory
with same filename should be recomputed and overwritten. If set to FALSE,
all modules are computed and existing files are quietly overwritten. If set to TRUE, existing files are
not modified, but non-existing files are computed and exported. Default set to FALSE.}
}
\value{
Returns the features dataframe.\cr
'Features' contains row number, plate name, well, id, total distance in pixels and mm, mean speed in pixels per sec (pps) and
in mm per second (mmps), max speed in pixels per sec (pps) and mm per sec (mmps), duration of
the whole sequence in frames (fr) and seconds (sec), activity (active seconds), framerate (fps) and calibration,
number of pixels in one mm (mm). If \code{coverage=TRUE}, then coverage is added in pxsq and mmsq. If \code{alphahull=TRUE}, then
alphahull area is added in pxsq and mmsq. If \code{spotratio=TRUE}, then outer spots, inner spots and spots ratio is added.
If \code{msd=TRUE}, then msd_px and msd_mm are added.
If more than one file was selected, a Combined-Features file is also exported.
If \code{exportplot = T}, then 8 figures are exported: EdgeMarks,
AssignedWells, Tracks, spotCoverage, spotMsd, spotAlphahull, spotDensity and spotRatio.
}
\description{
Tracks single larvae in 24 or 48 well assay plates from xy spot data
and computes distance and speed. Exports data as text files and generates plots.
}
\details{
The quality of tracks almost completely depends on the image thresholding and subsequent xy data.\cr
\strong{files}\cr
The input files must be tab-delimited decimal as dot (.) text files. The file
must contain a minimum of 3 columns named x, y and slice. x is a numeric
indicating x coordinate and y is a numeric indicating y coordinate of each spot.
slice indicates the frame number for each spot. Extra columns are not used.\cr
\strong{markededges}\cr
\code{markededges = T} indicates that the four corners of the plate have been marked in frame 1. They
will be used for plate alignment and plotting and will be removed from analyses.\cr
\strong{coverage}\cr
Area covered by each larvae in their respective well across the whole duration.
Computed from convex hull of points. The polygon area is computed based on the
function \code{polyarea} from package \code{pracma}.\cr
\strong{msd}\cr
The minimum spanning distance based on the minimum spanning network of point cloud.
Computed using \code{spantree} from package \code{vegan}. Very slow. Set to FALSE by default.\cr
\strong{alphahull}\cr
The alphahull of point cloud based on the alphavalue. The function \code{ahull} from
package \code{alphahull}. Smaller alphavalue produces more gaps in spot cloud.\cr
\strong{spotdensity}\cr
Computes 2D kernal density of spots and generates an image file.\cr
\strong{spotratio}\cr
Spots in the periphery of the well and spots in the centre of the well are computed. A ratio is computed.
\strong{filenamediscard}\cr
The file name of the input file is used on plots and text files for identification.
Part of the filename to be removed can be indicated here. '.txt' is removed by default.\cr
\strong{follow}\cr
Set as 'none', 'interactive' or 'track'.
In 'interactive' mode, a plot is shown at every frame and waits for user input.
In 'track' mode, the track path creation is shown in real-time.
}

