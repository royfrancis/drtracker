# LARVAE TRACKING AND FEATURE CALCULATION
# Roy Mathew Francis
# v1.1.0
# 21-Dec-2015

# Assumptions
# 24 or 48 well plate
# plate is constant throughout the film. absolutely no movement.
# edges and larvae positions are marked in frame 1
# single larvae per well
# if larval position is missing in a frame, old position is duplicated
# if two positions are available, nearest is chosen.
# kernal density colours are by plate

# To Do
# In package detection, better detection
# Faster. Cpp funtions, multicore design

# load libraries
# library(alphahull) #ahull()
# library(Cairo) #graphics
# library(fossil) #dino.mst
# library(ggplot2) #plots
# library(plyr) #rbind.fill()
# library(RANN) #nn2()
# library(reshape2) #dcast()
# library(SDMTools) #pnt.in.poly()

#check packages
pkgs <- c("alphahull","fossil","ggplot2","plyr","RANN","reshape2","SDMTools")
if(any(!pkgs %in% installed.packages()))
{
  stop(paste0("Packages '",paste0(pkgs[which(!pkgs %in% installed.packages())],collapse=", "),"' not installed.\n"))
}
rm(pkgs)

# FUNCTION LTRACK
#' Track larval movement in 24 or 48 well assay plates.
#' @description Tracks single larvae in 24 or 48 well assay plates from xy spot data
#' and computes distance and speed. Exports data as text files and generates plots.
#' @param files A character or vector of paths or filenames. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param markededges A logical indicating if corners of plates have been marked. See details.
#' @param fps A numeric indicating framerate of video in number of frames per second.
#' @param mm A numeric indicating number of pixels in 1 mm.
#' @param framelinking A logical indicating if spots should be linked frame to frame.
#' @param follow A character indicating if the algorithm must be followed at every step. Options are 'none', 'interactive' or 'track'. See details.
#' @param trackfeatures A logical indicating if track features should be computed.
#' @param activitydist A numeric indicating distance in mm after which the larvae is considered active.
#' @param spotcoverage A logical indicating if the coverge must be computed. See details.
#' @param spotalphahull A logical indicating if alphahull of point cloud is to be calculated. See details.
#' @param alphavalue A numeric indicating alpha value for alphahull. See details.
#' @param spotdensity A logical indicating if 2D kernal density is to be estimated. See details.
#' @param spotmsd A logical indicating if the minimum spanning distance should be computed. See details.
#' @param spotratio A logical indicating if spotratio should be computed. See details.
#' @param filenamediscard A character for part of the filename to be removed.
#' @param exportdata A logical indicating if data tables must be exported to working directory.
#' @param exportplot A logical if results must be plotted and exported to working directory.
#' @param useexisting A logical indicating if existing files in working directory
#' with same filename should be recomputed and overwritten. If set to FALSE,
#' all modules are computed and existing files are quietly overwritten. If set to TRUE, existing files are
#' not modified, but non-existing files are computed and exported. Default set to FALSE.

#' @return Returns the features dataframe.\cr
#' 'Features' contains row number, plate name, well, id, total distance in pixels and mm, mean speed in pixels per sec (pps) and
#' in mm per second (mmps), max speed in pixels per sec (pps) and mm per sec (mmps), duration of
#' the whole sequence in frames (fr) and seconds (sec), activity (active seconds), framerate (fps) and calibration,
#' number of pixels in one mm (mm). If \code{coverage=TRUE}, then coverage is added in pxsq and mmsq. If \code{alphahull=TRUE}, then
#' alphahull area is added in pxsq and mmsq. If \code{spotratio=TRUE}, then outer spots, inner spots and spots ratio is added.
#' If \code{msd=TRUE}, then msd_px and msd_mm are added.
#' If more than one file was selected, a Combined-Features file is also exported.
#' If \code{exportplot = T}, then 3-7 figures are exported: EdgeMarks (if edge marks are used),
#' wells & spots, wells & tracks, coverage, msd, alphahull and distance.
#' @details
#' The quality of tracks almost completely depends on the image thresholding and subsequent xy data.\cr
#' \strong{files}\cr
#' The input files must be tab-delimited decimal as dot (.) text files. The file
#' must contain a minimum of 3 columns named x, y and slice. x is a numeric
#' indicating x coordinate and y is a numeric indicating y coordinate of each spot.
#' slice indicates the frame number for each spot. Extra columns are not used.\cr
#' \strong{markededges}\cr
#' \code{markededges = T} indicates that the four corners of the plate have been marked in frame 1. They
#' will be used for plate alignment and plotting and will be removed from analyses.\cr
#' \strong{coverage}\cr
#' Area covered by each larvae in their respective well across the whole duration.
#' Computed from convex hull of points. The polygon area is computed based on the
#' function \code{polyarea} from package \code{pracma}.\cr
#' \strong{msd}\cr
#' The minimum spanning distance based on the minimum spanning network of point cloud.
#' Computed using \code{spantree} from package \code{vegan}. Very slow. Set to FALSE by default.
#' \strong{alphahull}\cr
#' The alphahull of point cloud based on the alphavalue. The function \code{ahull} from
#' package \code{alphahull}. Smaller alphavalue produces more gaps in spot cloud.
#' \strong{spotdensity}\cr
#' Computes 2D kernal density of spots and generates an image file.
#' \strong{spotratio}\cr
#' Spots in the periphery of the well and spots in the centre of the well are computed. A ratio is computed.
#' \strong{filenamediscard}\cr
#' The file name of the input file is used on plots and text files for identification.
#' Part of the filename to be removed can be indicated here. '.txt' is removed by default.\cr
#' \strong{follow}\cr
#' Set as 'none', 'interactive' or 'track'.
#' In 'interactive' mode, a plot is shown at every frame and waits for user input.
#' In 'track' mode, the track path creation is shown in real-time.
#' @export
#' @import alphahull
#' @import fossil
#' @import RANN
#' @import reshape2
#' @import SDMTools
#'
ltrack <- function(files = NULL, wells = 24, markededges = TRUE, fps = 25, mm = 5.4,
                   framelinking=TRUE, follow = "none",
                   trackfeatures = TRUE, activitydist = 5,
                   spotcoverage = TRUE,
                   spotalphahull = TRUE, alphavalue = 4,
                   spotdensity = TRUE,
                   spotmsd = FALSE,
                   spotratio = TRUE,
                   filenamediscard = ".txt", exportdata = TRUE, exportplot = TRUE,
                   useexisting = FALSE)
{
  currtime <- Sys.time()

  cat(paste0(files,"\n"))
  # Argument checks
  if(any(is.null(files))) stop("No input files.\n")
  if(any(is.na(files)) | length(files)==0) stop("No input files.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.logical(markededges)) stop("Argument 'markededges' is not set correctly. Set as TRUE or FALSE.\n")
  if(!is.numeric(fps)) stop("Argument 'fps' not set correctly. Assign a number.\n")
  if(!is.numeric(mm)) stop("Argument 'mm' not set correctly. Assign a number.\n")
  if(!is.character(filenamediscard)) stop("Argument 'filenamediscard' not set correctly. Assign a character.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")

  if(!is.logical(spotcoverage)) stop("Argument 'coverage' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(spotalphahull)) stop("Argument 'coverage' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(spotdensity)) stop("Argument 'coverage' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(spotmsd)) stop("Argument 'coverage' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(spotratio)) stop("Argument 'coverage' not set correctly. Set as TRUE or FALSE.\n")

  if((follow != "none") && (follow != "interactive") && (follow != "track")) stop("Set argument 'follow' as 'none', 'interactive' or 'track'.\n")

  cat(paste0("Note values in use; Wells: ",wells,", Framerate: ",fps," fps, 1 mm = ",mm," pixels.\n"))

  #looping over selected files
  fnames <- sub(filenamediscard,"",basename(files))
  flen <- length(files)
  cdatalist <- vector("list", length = flen)
  fileloop=1
  for(fileloop in 1:flen)
  {
    cat(paste0("\nFile ",fileloop," of ",flen,".\n"))
    cat(paste0("Current file: ",files[fileloop],".\n"))

    dfa <- data.frame(row = 1:wells, plate = rep(fnames[fileloop],wells), well = 1:wells, stringsAsFactors=FALSE)

    #COMPUTATION
    #assign wells
    if(useexisting)
    {
      if(paste0(fnames[fileloop],"-AssignedWells.txt") %in% list.files())
      {
        cat(paste0("Found ",fnames[fileloop],"-AssignedWells.txt.\n"))
        aw <- read.delim(paste0(fnames[fileloop],"-AssignedWells.txt"), header = TRUE, stringsAsFactors = FALSE)
      }else{
        cat(paste0("Assigning wells...\n"))
        aw <- assignWells(files = files[fileloop], wells = wells, markededges = markededges,
                          filenamediscard = filenamediscard, exportdata = TRUE,
                          exportplot = exportplot, quiet = TRUE)
      }
    }
    if(!useexisting)
    {
      cat(paste0("Assigning wells...\n"))
      aw <- assignWells(files = files[fileloop], wells = wells, markededges = markededges,
                        filenamediscard = filenamediscard, exportdata = TRUE,
                        exportplot = exportplot, quiet = TRUE)
    }

    #framelinking
    if(framelinking)
    {
      if(useexisting)
      {
        if(paste0(fnames[fileloop],"-Tracks.txt") %in% list.files())
        {
          cat(paste0("Found ",fnames[fileloop],"-Tracks.txt.\n"))
          ff <- read.delim(paste0(fnames[fileloop],"-Tracks.txt"), header = TRUE, stringsAsFactors = FALSE)
        }else{
          cat(paste0("Linking frames...\n"))
          ff <- linkFrames(dframe = aw, wells = wells, exportdata = TRUE, exportplot = exportplot,
                           follow = follow, quiet = TRUE)
        }
      }
      if(!useexisting)
      {
        cat(paste0("Linking frames...\n"))
        ff <- linkFrames(dframe = aw, wells = wells, exportdata = TRUE, exportplot = exportplot,
                         follow = follow, quiet = TRUE)
      }
    }

    #stop if tracks file not found
    if(!paste0(fnames[fileloop],"-EdgeMarks.txt") %in% list.files()) stop(paste0("File ",fnames[fileloop],"-EdgeMarks.txt not found.\n"))

    if(trackfeatures)
    {
      cat(paste0("Computing track features...\n"))
      #track features
      ff <- trackFeatures(files = paste0(fnames[fileloop],"-Tracks.txt"), wells = wells, fps = fps, mm = mm,
                           activitydist = activitydist, exportdata = exportdata,
                           exportplot = exportplot, quiet = TRUE)
      dfa <- cbind(dfa,ff[,-c(1,3)])
    }

    #spot coverage
    if(spotcoverage)
    {
      cat(paste0("Computing spot coverage...\n"))
      ff <- spotCoverage(files = paste0(fnames[fileloop],"-Tracks.txt"), wells = wells, mm = mm,
                          exportdata = exportdata, exportplot = exportplot, quiet = TRUE)
      dfa <- cbind(dfa,ff[,-c(1:3)])
    }

    #spot alphahull
    if(spotalphahull)
    {
      cat(paste0("Computing spot alphahulls...\n"))
      ff <- spotAlphahull(files = paste0(fnames[fileloop],"-Tracks.txt"), wells = wells, mm = mm, alphavalue = alphavalue,
                           exportdata = exportdata, exportplot = exportplot, quiet = TRUE)
      dfa <- cbind(dfa,ff[,-c(1:3)])
    }

    #spot msd
    if(spotmsd)
    {
      cat(paste0("Computing spot minimum spanning tree...\n"))
      ff <- spotMsd(files = paste0(fnames[fileloop],"-Tracks.txt"), wells = wells, mm = mm,
                     exportdata = exportdata, exportplot = exportplot, quiet = TRUE)
      dfa <- cbind(dfa,ff[,-c(1:3)])
    }

    #spot density
    if(spotdensity)
    {
      cat(paste0("Computing spot density...\n"))
      spotDensity(files = paste0(fnames[fileloop],"-Tracks.txt"), wells = wells, mm = mm,
                         exportdata = exportdata, exportplot = exportplot, quiet = TRUE)
    }

    #spot ratio
    if(spotratio)
    {
      cat(paste0("Computing spot ratio...\n"))
      ff <- spotRatio(files = paste0(fnames[fileloop],"-Tracks.txt"), wells = wells, mm = mm,
                       exportdata = exportdata, exportplot = exportplot, quiet = TRUE)
      dfa <- cbind(dfa,ff[,-c(1:3)])
    }

    #---------------------------------------------------------------------------
    if(exportdata)
    {
      write.table(x = dfa, file = paste0(fnames[fileloop],"-Features.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fnames[fileloop],"-Features.txt exported.\n"))
    }

    cdatalist[[fileloop]] <- dfa

  }

  cdatadf <- plyr::rbind.fill(cdatalist)

  if(exportdata)
  {
    if(flen > 1)
    {
      write.table(x = cdatadf, file = paste0("Combined-Features.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0("Combined-Features.txt exported.\n"))
    }
  }

  cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
}

#-------------------------------------------------------------------------------

# FUNCTION ASSIGNWELLS
#' Takes xy data and assigns to wells for 24 or 48 well plate.
#' @description Takes xy data and assigns to wells for 24 or 48 well plate.
#' @param files A character or vector of paths or filenames. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param markededges A logical indicating if corners of plates have been marked. See details.
#' @param filenamediscard A character for part of the filename to be removed. '.txt' is removed by default.
#' @param exportdata A logical indicating if data table must be exported as a text file to the working directory.
#' @param exportplot A logical if results must be plotted as a figure and exported to working directory.
#' @param quiet A logical indicating if messages should be printed to console during the run. If \code{FALSE}, all output to console is killed except progress bar.
#' @return Returns a dataframe with original data (x,y,slice) along with a new columns row, plate and well.
#' plate column is the filename without filenamediscard. The well column contains the well assignment.
#' If \code{exportdata=TRUE}, then a text file with the dataframe is exported in the working directory for each input file.
#' @details
#' \strong{files}\cr
#' A character or vector of paths to files. On windows, use \code{choose.files()} for interactive selection.
#' The input files must be tab-delimited decimal as dot (.) text files. The file
#' must contain a minimum of 3 columns named x, y and slice. x is a numeric
#' indicating x coordinate and y is a numeric indicating y coordinate of each spot.
#' slice indicates the frame number for each spot. Any extra columns are not used.\cr
#' \strong{markededges}\cr
#' \code{markededges=T} indicates that the four corners of the plate have been marked in frame 1. They
#' will be used for plate alignment and plotting and will be removed from other analyses.\cr
#' @export
#' @import RANN

assignWells <- function(files = NULL, wells = 24, markededges = TRUE, filenamediscard = ".txt",
                        exportdata = TRUE, exportplot = TRUE, quiet = FALSE)
{
  if(!quiet) currtime <- Sys.time()

  # Argument checks
  if(any(is.null(files))) stop("No input files.\n")
  if(any(any(is.na(files)) | any(length(files)==0))) stop("No input files.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.logical(markededges)) stop("Argument 'markededges' is not set correctly. Set as TRUE or FALSE.\n")
  if(!is.character(filenamediscard)) stop("Argument 'filenamediscard' not set correctly. Assign a character.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")

  #looping over selected files
  flen <- length(files)
  dframelist <- vector("list",length = flen)
  fileloop=1
  for(fileloop in 1:flen)
  {
    if(!quiet) cat(paste0("Assigning wells. File ",fileloop," of ",flen,".\n"))
    if(!quiet) cat(paste0("Current file: ",files[fileloop],".\n"))

    fname <- sub(filenamediscard,"",basename(files[fileloop]))
    dframe <- read.delim(files[fileloop],header = T,stringsAsFactors = F)
    colnames(dframe) <- tolower(colnames(dframe))
    if(!all(c("x","y","slice") %in% colnames(dframe))) stop("Columns 'x', 'y' or 'slice' not found in input file.\n")
    if(any(c("id","linktype","frame") %in% colnames(dframe))) warning("Input file may be incorrect.\n")
    #extradf <- subset(dframe, select = -c(x,y,slice))
    dframe <- dframe[,c("x","y","slice")]
    #round all coordinates
    #dframe$x <- round(dframe$x,2)
    #dframe$y <- round(dframe$y,2)

    #-----------------------------------------------------------------------------

    if(wells == 24)
    {
      prows <- 4
      pcols <- 6
    }

    if(wells == 48)
    {
      prows <- 6
      pcols <- 8
    }

    if(wells == 96)
    {
      prows <- 8
      pcols <- 12
    }

    # creating well layout
    verticalrange <- range(dframe$y)
    horizontalrange <- range(dframe$x)
    vertpos <- cumsum(c(verticalrange[1],rep((diff(verticalrange)/prows),prows)))
    horizpos <- cumsum(c(horizontalrange[1],rep((diff(horizontalrange)/pcols),pcols)))

    h1 <- vector()
    h2 <- vector()
    v1 <- vector()
    v2 <- vector()
    for(i in 1:prows)
    {
      for(j in 1:pcols)
      {
        h1 <- c(h1,horizpos[j])
        h2 <- c(h2,horizpos[j+1])
        v1 <- c(v1,vertpos[i])
        v2 <- c(v2,vertpos[i+1])
      }
    }

    wellsdf <- data.frame(well = factor(1:wells),h1 = h1,h2 = h2,v1 = v1,v2 = v2, stringsAsFactors = FALSE)
    wellsdf$x <- wellsdf$h1+((wellsdf$h2-wellsdf$h1)/2)
    wellsdf$y <- wellsdf$v1+((wellsdf$v2-wellsdf$v1)/2)
    wellsdf$labx <- wellsdf$h1+(0.08*(diff(horizontalrange)/prows))
    wellsdf$laby <- wellsdf$v1+(0.15*(diff(verticalrange)/pcols))

    # handling marked edges
    if(markededges)
    {
      dframe1 <- subset(dframe,dframe$slice == 1)
      edges <- chull(x = dframe1$x,y = dframe1$y)
      edgedf <- dframe[edges,]
      rm(dframe1)
      # remove edge spots
      dframe <- dframe[-edges,]
      if(!quiet) cat(paste0(length(edges)," edge spots removed.\n"))

      if(exportdata)
      {
        write.table(x = edgedf, file = paste0(fname,"-EdgeMarks.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
        cat(paste0(fname,"-EdgeMarks.txt exported.\n"))
      }

      # plot edge spots
      if(exportplot)
      {
        p <- ggplot(edgedf,aes(x,y))+
          geom_point(shape = 1,size = 3,na.rm = TRUE)+
          ggtitle(paste0(wells," Well Plate: ",fname,"  |  ",length(edges)," EdgeMarks"))+
          scale_y_reverse()+
          theme_bw(base_size = 5)+
          theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 4.5,face="bold"),
                axis.title = element_blank(),panel.border = element_blank(),
                axis.text = element_text(colour = "grey40"),
                axis.ticks = element_line(colour = "grey40",size = 0.3))
        ggsave(filename = paste0(fname,"-EdgeMarks.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fname,"-EdgeMarks.png exported.\n"))
      }
    }

    # compute well positions using nn2
    #cat(paste0("Computing well positions.\n"))
    nrowdf <- nrow(dframe)
    if(!quiet) cat(paste0("Assigning ",nrowdf," spots to ",wells," wells...\n"))
    pb <- txtProgressBar(min = 0, max = nrowdf, style = 3)
    wellvec <- vector(length = nrowdf)

    for(i in 1:nrowdf)
    {
      wellvec[i] <- as.numeric(RANN::nn2(wellsdf[,c("x","y")],dframe[i,c("x","y")],k = 1)$nn.idx)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    #add well assignment to df
    dframe$well <- as.numeric(as.character((wellvec)))
    #bind rows, plate name, df and extradf
    dframe <- cbind(data.frame(row=1:nrowdf,plate=rep(fname,nrowdf)),dframe)
    #if(ncol(extradf) > 0) dframe <- cbind(dframe,extradf)

    if(exportdata)
    {
      write.table(x = dframe, file = paste0(fname,"-AssignedWells.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fname,"-AssignedWells.txt exported.\n"))
    }

    # plot spot layout after well assignment
    if(exportplot)
    {
      if(wells == 24) txtsz <- 1.5
      if(wells == 48) txtsz <- 1.3
      if(wells == 96) txtsz <- 1.1
      if(wells == 24) cent <- 7
      if(wells == 48) cent <- 4
      if(wells == 96) cent <- 2
      p <- ggplot()+
        geom_blank(data = edgedf,aes(x,y))+
        geom_rect(data = wellsdf,aes(xmin = h1,xmax = h2,ymin = v1,ymax = v2),colour = "grey90",fill = "white",size = 0.3,na.rm = TRUE)+
        geom_point(data = dframe,aes(x,y,colour = factor(well)),size = 0.3,na.rm = TRUE)+
        geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)
      p <- p + scale_x_continuous(expand = c(0,0))+
        scale_y_reverse(expand = c(0,0))+
        theme_bw(base_size = 5)+
        labs(title = paste0(wells, " Well Plate: ",fname,"  |  Wells & Spots"))+
        coord_fixed()+
        theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 4.5,face="bold"),legend.position = "none",
              axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
              axis.line = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
              panel.border = element_blank())
      ggsave(filename = paste0(fname,"-AssignedWells.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
      cat(paste0(fname,"-AssignedWells.png exported.\n"))
    }

    dframelist[[fileloop]] <- dframe
  }

  if(!quiet) cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
  if(flen > 1) return(plyr::rbind.fill(dframelist))
  if(flen == 1) return(dframe)
}

#-------------------------------------------------------------------------------

# FUNCTION LINKFRAMES
#' Takes xy data with wells and links spots frame to frame.
#' @description Takes xy data with wells and links spots frame to frame.
#' @param files A character or vector of paths or filenames. A text file exported from \code{assignWells()} is suitable. See details.
#' @param dframe A dataframe output from \code{assignWells()}. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param filenamediscard A character for part of the filename to be removed. '-AssignedWells.txt' is removed by default.
#' @param exportdata A logical indicating if data table must be exported as a text file to the working directory.
#' @param exportplot A logical if results must be plotted as a figure and exported to working directory.
#' @param follow A character indicating if the algorithm must be followed at every step. See details.
#' @param quiet A logical indicating if messages should be printed to console during the run. If \code{FALSE}, all output to console is killed except progress bar.
#' @return Returns a dataframe with original data (plate,x,y,well) along with a new columns id, linktype and frame.
#' plate column is the filename without filenamediscard. The well column contains the well assignment. If number
#' of input files > 1, then the dataframe RESULTS FOR THE LAST INPUT FILE is returned. If \code{exportdata=TRUE}, then a text file ending in '-Tracks' with
#' dataframe is exported in the working directory for each input file. This '-Tracks' text file is used as input for several other functions.
#' @details
#' Note that the 'xx-EdgeMarks.txt' file must be present in the working directory for the plot exports to work.
#' \strong{Input}\cr
#' The input for this function can be a vector of filenames passed to argument \code{files} or
#' a dataframe passed to argument \code{dframe}. Use any one and not both.
#' \strong{files}\cr
#' A character or vector of paths to files. On windows, use \code{choose.files()} for interactive selection.
#' The input files must be tab-delimited decimal as dot (.) text files. The file
#' must contain a minimum of 4 columns named plate, x, y and well. plate is a character
#' denoting the filename. x is a numeric indicating x coordinate and y is a numeric
#' indicating y coordinate of each spot. well indicates the well number to which each spot
#' is assigned. Any extra columns are not used.\cr
#' \strong{dframe}\cr
#' A dataframe containing a minimum of 4 columns named plate, x, y and well. plate is a character
#' denoting the filename. x is a numeric indicating x coordinate and y is a numeric
#' indicating y coordinate of each spot. well indicates the well number to which each spot
#' is assigned. Any extra columns are not used.\cr
#' \strong{follow}\cr
#' Set as 'none', 'interactive' or 'track'.
#' In 'interactive' mode, a plot is shown at every frame and waits for user input.
#' In 'track' mode, the track path creation is shown in real-time.
#' @export
#' @import RANN

linkFrames <- function(files = NULL, dframe = NULL, wells = 24, filenamediscard = "-AssignedWells.txt",
                       exportdata = TRUE, exportplot = TRUE, follow = "none", quiet = FALSE)
{
  if(!quiet) currtime <- Sys.time()
  if(is.null(files) && is.null(dframe)) stop("Arguments 'files' and 'dframe' are both empty. Use any one.\n")
  if(!is.null(files) && !is.null(dframe)) stop("Arguments 'files' and 'dframe' are both in use. Use any one.\n")
  if(is.null(dframe) && !is.null(files)) {typefiles <- TRUE}else{typefiles <- FALSE}

  # Argument checks
  if(typefiles) {if(any(is.null(files)) | any(is.na(files)) | length(files)==0) stop("No input files.\n")}
  if(!typefiles)
  {
    colnames(dframe) <- tolower(colnames(dframe))
    if(!all(c("plate","x","y","well") %in% colnames(dframe))) stop("Columns 'plate', 'x', 'y' or 'well' not found in input file.\n")
  }

  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")

  if((follow != "none") && (follow != "interactive") && (follow != "track")) stop("Set argument 'follow' as 'none', 'interactive' or 'track'.\n")

  #if files input
  if(typefiles)
  {
    flen <- length(files)
    fnames <- sub(filenamediscard,"",basename(files))
  }

  #if dframe input
  if(!typefiles)
  {
    fnames <- sub(filenamediscard,"",unique(as.character(dframe$plate)))
    flen <- length(fnames)
  }

  fileloop=1
  for(fileloop in 1:flen)
  {
    if(!quiet) cat(paste0("Linking frames. File ",fileloop," of ",flen,".\n"))
    if(!quiet) cat(paste0("Current file: ",fnames[fileloop],".\n"))

    if(typefiles)
    {
      dframefile <- read.delim(files[fileloop],header = TRUE,stringsAsFactors = FALSE)
      colnames(dframefile) <- tolower(colnames(dframefile))
      if(!all(c("plate","x","y","well") %in% colnames(dframefile))) stop("Columns 'plate', 'x', 'y' or 'well' not found in input file.\n")
    }

    if(!typefiles)
    {
      dframefile <- subset(dframe,dframe$plate == fnames[fileloop])
    }

    #loop over frames
    framelist <- split(dframefile,dframefile$slice)
    dflist <- vector("list",length = length(framelist))
    numberofpoints <- length(framelist)
    pb <- txtProgressBar(min = 0, max = numberofpoints, style = 3)
    frameloop=1
    for(frameloop in 1:numberofpoints)
    {
      setTxtProgressBar(pb, frameloop)
      #if(frameloop==51) stop("Error.")
      if(frameloop == numberofpoints) break
      #get exact points from frame 1
      if (frameloop == 1)
      {
        currentframe <- framelist[[frameloop]]
        currentframe$id <- 1:nrow(currentframe)
        #ids <- currentframe$id
      }
      nextframe <- framelist[[frameloop+1]]

      #loop over points
      xylist <- vector("list",length = nrow(currentframe))
      pointloop = 1
      for(pointloop in 1:nrow(currentframe))
      {
        #cat(paste0("Current point: ",pointloop,"\n"))
        # select a point
        searchpointpos <- currentframe[pointloop,]

        # look for a point in same well in next frame
        foundpointtemp <- subset(nextframe,well == searchpointpos$well)

        #check if any match was found
        if(nrow(foundpointtemp) > 0)
        {
          #if one point was found
          if(nrow(foundpointtemp) == 1)
          {
            foundpointpos <- foundpointtemp
            foundpointpos$id <- searchpointpos$id
            foundpointpos$linktype <- "Single"
          }

          #if more than one point was found
          if(nrow(foundpointtemp) > 1)
          {
            if(nrow(foundpointtemp) > 1)
            {
              foundpointpos <- foundpointtemp[as.numeric(RANN::nn2(foundpointtemp[,c("x","y")],searchpointpos[,c("x","y")],k = 1)$nn.idx),]
              foundpointpos$id <- searchpointpos$id
              foundpointpos$linktype <- "Nearest"
            }
          }

        }else{
          # if no match was found, duplicate point in next frame
          foundpointpos <- searchpointpos
          foundpointpos$slice <- frameloop
          foundpointpos$linktype <- "Duplicated"
        }

        # add new pos to xylist
        xylist[[pointloop]] <- foundpointpos
        #if((follow == "interactive") | (follow == "track")) segments(x0 = searchpointpos$x, y0 = searchpointpos$y,x1 = foundpointpos$x,y1 = foundpointpos$y,lty = 3)
      }
      tempdf <- plyr::rbind.fill(xylist)

      #plot for follow
      if(follow == "interactive")
      {
        pframe <- rbind(currentframe[,c("x","y","slice","id")],tempdf[,c("x","y","slice","id")])
        pframe$frame <- c(rep(currentframe$slice[1],nrow(currentframe)),rep(currentframe$slice[1]+1,nrow(tempdf)))
        #pframe$linktype <- factor(c(rep("None",nrow(currentframe)),as.character(tempdf$linktype)),levels=c("None","Single","Nearest","Duplicated"))
        p <- ggplot()+
          geom_point(data=pframe,aes(x=x,y=y,shape=factor(frame),col=factor(frame)),size=6,alpha=0.6,na.rm = TRUE)+
          geom_line(data=pframe,aes(x=x,y=y,group=factor(id)),linetype=2,na.rm = TRUE)+
          geom_text(data=currentframe,aes(x=x+20,y=y,label=id),size=3,na.rm = TRUE)+
          scale_shape_manual(values=c(7,15,8,9))+
          scale_fill_manual(values=c("coral","steelblue","darkolivegreen3"))+
          ggtitle(paste0("Current Frame: ",frameloop))+
          labs(shape="Frame",colour="Frame")+
          theme_bw()+
          theme(plot.title=element_text(hjust=0,line=0.5),axis.title=element_blank(),
                panel.border=element_blank())

        print(p)
        readline("Press enter to continue:")
      }

      #add first frame
      if(frameloop == 1)
      {

        dflist[[frameloop]] <- currentframe
        dflist[[frameloop]]$id <- 1:nrow(currentframe)
        dflist[[frameloop]]$linktype <- rep(NA,nrow(currentframe))
        dflist[[frameloop]]$frame <- rep(1,nrow(currentframe))
      }

      currentframe <- tempdf
      #add frame number
      tempdf$frame <- rep(frameloop+1,nrow(currentframe))
      #save df to list
      dflist[[frameloop+1]] <- tempdf

      if(follow == "track" && ((frameloop %% 10) == 0))
      {
        dftrack <- plyr::arrange(plyr::rbind.fill(dflist[1:frameloop]),frame)
        p <- ggplot()+
          geom_path(data=dftrack,aes(x=x,y=y,group=factor(id),colour=factor(well)),linetype=1,size=0.8)+
          #geom_text(data=dftrack,aes(x=x+20,y=y,label=id),size=3)+
          ggtitle(paste0("Current Frame: ",frameloop))+
          theme_bw()+
          theme(plot.title=element_text(hjust=0,line=0.5),axis.title=element_blank(),
                panel.border=element_blank(),legend.position="none")

        print(p)
        #Sys.sleep(0.2)
        rm(dftrack)
      }

    }
    close(pb)
    df1 <- plyr::rbind.fill(dflist)

    #correct for missing wells
    wellunique <- sort(as.numeric(unique(as.character(df1$well))))
    if(length(wellunique) != wells)
    {
      #compare wells in df1 to full seq of wells
      wdiff <- setdiff(1:wells, wellunique)
      wlen <- length(wdiff)
      #if missing wells, add NA rows
      if(wlen > 0)
      {
        tempdf <- data.frame(row=rep(NA,wlen),plate=rep(fnames[fileloop],wlen),x=rep(NA,wlen),
                   y=rep(NA,wlen),slice=rep(NA,wlen),well=wdiff,id=rep(NA,wlen),
                   linktype=rep(NA,wlen),frame=rep(NA,wlen))
        df1 <- rbind(df1,tempdf)
        rm(tempdf)
      }
    }
    #replace rows
    df1$row <- 1:nrow(df1)

    if(exportdata)
    {
      write.table(x = df1, file = paste0(fnames[fileloop],"-Tracks.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fnames[fileloop],"-Tracks.txt exported.\n"))
    }

    if(exportplot)
    {
      if(paste0(fnames[fileloop],"-EdgeMarks.txt") %in% list.files())
      {
        edgedf <- read.delim(paste0(fnames[fileloop],"-EdgeMarks.txt"),header=T,stringsAsFactors = FALSE)
        wellsdf <- gridWells(dframe = edgedf, wells = wells)
        df1$linktype <- factor(as.character(df1$linktype),levels=c("Single","Duplicated","Nearest"))

        if(wells == 24)
        {
          txtsz <- 1.5
          txtszbp <- 1.8
        }

        if(wells == 48)
        {
          txtsz <- 1.3
          txtszbp <- 1.5
        }

        if(wells == 96)
        {
          txtsz <- 1.1
          txtszbp <- 1.2
        }

        if(wells == 24) cent <- 7
        if(wells == 48) cent <- 4
        if(wells == 96) cent <- 2

        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          geom_path(data = df1,aes(x = x,y = y,colour = linktype,group = factor(id)),size = 0.25,alpha=0.8,na.rm = TRUE)+
          geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)+
          theme_bw()+
          scale_x_continuous(expand = c(0,0)) +
          scale_y_reverse(expand = c(0,0))+
          scale_colour_manual(values=c("#31a354","#3182bd","#de2d26"))+
          labs(title = paste0(wells, " Well Plate: ",fnames[fileloop],"  |  Wells & Tracks"))+
          coord_fixed()+
          theme_bw(base_size = 5)+
          theme(legend.position = "none",axis.ticks = element_blank(),legend.text = element_blank(),
                axis.text = element_blank(),axis.title = element_blank(),legend.title = element_blank(),
                plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 4.5,face="bold"),panel.border = element_blank(),
                panel.background = element_blank(),panel.grid = element_blank())
        ggsave(filename = paste0(fnames[fileloop],"-Tracks.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fnames[fileloop],"-Tracks.png exported.\n"))
      }else{
        warning(paste0("Tracks plot not exported since file ",paste0(fnames[fileloop],"EdgeMarks.txt") ," not found in working directory.\n"))
      }
    }
  }

  if(!quiet) cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
  return(df1)
}

#-------------------------------------------------------------------------------

# FUNCTION TRACKFEATURES
#' Computes features for tracks such as distance, speed etc.
#' @description Computes features for tracks such as distance, speed etc.
#' @param files A character or vector of paths or filenames. A text file exported from \code{linkFrames()} is suitable. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param fps A numeric indicating framerate of video in number of frames per second.
#' @param mm A numeric indicating number of pixels in 1 mm.
#' @param activitydist A numeric indicating distance in mm after which the larvae is considered active.
#' @param filenamediscard A character for part of the filename to be removed. '-Tracks.txt' is removed by default.
#' @param exportdata A logical indicating if data table must be exported as a text file to the working directory.
#' @param exportplot A logical if results must be plotted as a figure and exported to working directory.
#' @param quiet A logical indicating if messages should be printed to console during the run. If \code{FALSE}, all output to console is killed except progress bar.
#' @return Returns a dataframe with columns: plate name, id, well, total distance in pixels and mm, mean speed in pixels per sec (pps) and
#' in mm per second (mmps), max speed in pixels per sec (pps) and mm per sec (mmps), activity (active seconds), duration of
#' the whole sequence in frames (fr) and seconds (sec), framerate (fps) and calibration,
#' number of pixels in one mm (mm).
#' @details
#' \strong{files}\cr
#' A character or vector of paths to files. On windows, use \code{choose.files()} for interactive selection.
#' The input files must be tab-delimited decimal as dot (.) text files. The file
#' must contain a minimum of 5 columns named plate, x, y, well and frame. plate is a character
#' denoting the filename. x is a numeric indicating x coordinate and y is a numeric
#' indicating y coordinate of each spot. well indicates the well number to which each spot
#' is assigned. frame indicates frame number for each spot. Any extra columns are not used.\cr
#' The number of wells are determined from input data.
#' @export

trackFeatures <- function(files = NULL, wells = 24, fps = 25, mm = 5.4,activitydist = 5,
                          filenamediscard = "-Tracks.txt", exportdata = TRUE,
                          exportplot = TRUE, quiet = FALSE)
{
  if(!quiet) currtime <- Sys.time()

  # Argument checks
  if(any(is.null(files))) stop("No input files.\n")
  if(any(is.na(files)) | length(files)==0) stop("No input files.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.numeric(fps)) stop("Argument 'fps' not set correctly. Assign a number.\n")
  if(!is.numeric(mm)) stop("Argument 'mm' not set correctly. Assign a number.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")

  if(!quiet) cat(paste0("Note values in use; Framerate: ",fps," fps, 1 mm = ",mm," pixels.\n"))

  flen <- length(files)
  fnames <- sub(filenamediscard,"",basename(files))
  fileloop=1
  for(fileloop in 1:flen)
  {
    if(!quiet) cat(paste0("Computing distance, speed and activity.. File ",fileloop," of ",flen,".\n"))
    if(!quiet) cat(paste0("Current file: ",fnames[fileloop],".\n"))

    df1 <- read.delim(files[fileloop],header = TRUE,stringsAsFactors = FALSE)
    colnames(df1) <- tolower(colnames(df1))
    if(!all(c("plate","x","y","well","frame") %in% colnames(df1))) stop("Columns 'plate', 'x', 'y', 'well' or 'frame' not found in input file.\n")
    pb <- txtProgressBar(min = 0, max = wells, style = 3)
    welllist <- vector("list",length = wells)
    secspeedppslist <- vector("list",length=wells)
    #loop over id
    loop1=1
    for(loop1 in 1:wells)
    {
      setTxtProgressBar(pb, loop1)
      currentwell <- subset(df1, df1$well == loop1)

      if(nrow(currentwell) > 2)
      {
        framelen <- nrow(currentwell)
        distancevector <- vector(length = framelen)
        loop2 = 1
        secspeed <- vector(length = fps)
        secspeedpps <- vector()
        secspeedmps <- vector()
        k = 1
        #loop over frames
        while(loop2 < framelen)
        {
          if(k == (fps+1))
          {
            # save speed pps
            secspeedpps <- c(secspeedpps,sum(secspeed))
            # save speed mps
            secspeedmps <- c(secspeedmps,(sum(secspeed)/mm))
            secspeed <- vector(length = fps)
            k <- 1
          }
          #calculate euclidian distance
          d <- sqrt(diff(c(currentwell$x[loop2],currentwell$x[loop2+1]))^2 + diff(c(currentwell$y[loop2],currentwell$y[loop2+1]))^2)
          secspeed[k] <- d
          distancevector[loop2] <- d
          loop2 <- loop2+1
          k = k+1
        }
        # sum total distance
        totdist <- sum(distancevector)

        # create track stats df
        welllist[[loop1]] <- data.frame(plate = fnames[fileloop],id = as.numeric(currentwell$id[1]),well = loop1,
                                        dist_px = round(totdist,0),dist_mm = round(totdist/mm,2),speed_mean_pps = round(mean(secspeedpps),2),
                                        speed_mean_mmps = round(mean(secspeedmps),2),speed_max_pps = round(max(secspeedpps),2),
                                        speed_max_mmps = round(max(secspeedmps),2),duration_fr = nrow(currentwell),
                                        duration_sec = round(nrow(currentwell)/fps,0), activity = length(which(secspeedmps > activitydist)),
                                        fps = fps,mm = mm,stringsAsFactors = FALSE)
      }else{
        # create track stats df
        welllist[[loop1]] <- data.frame(plate = fnames[fileloop],id = as.numeric(currentwell$id[1]),well = loop1,
                                        dist_px = NA,dist_mm = NA,speed_mean_pps = NA,
                                        speed_mean_mmps = NA,speed_max_pps = NA,
                                        speed_max_mmps = NA,duration_fr = NA,
                                        duration_sec = NA, activity = NA,
                                        fps = fps,mm = mm,stringsAsFactors = FALSE)
      }


    }
    close(pb)
    df2 <- plyr::arrange(plyr::rbind.fill(welllist),well)

    if(exportdata)
    {
      write.table(x = df2, file = paste0(fnames[fileloop],"-TrackFeatures.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fnames[fileloop],"-TrackFeatures.txt exported.\n"))
    }
  }

  if(!quiet) cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
  return(df2)
}

#-------------------------------------------------------------------------------

# FUNCTION GRIDWELLS
#' Internal: Computes plate layout and wells from marked data.
#' @description Internal: Computes plate layout and wells from marked data.
#' @param dframe A dataframe with x and y coordinates.
#' @param wells An integer indicating number of wells on the plate.
#' @return Returns a dataframe with columns:
#' well,horizontal position 1 (h1), horizontal position 2 (h2), vertical position 1 (v1),
#' vertical position 2 (v2), well midpoint x (x), well midpoint y (y), top left label position x (labx),
#' top left label position y (laby) and bottom label y (laby1).
#' @details
#' The 4 marked edges of the plate is used to compute well layout positions. The resulting output is used in plotting and other functions.
#' @export

gridWells <- function(dframe = NULL, wells = NULL)
{
  # Argument checks
  if(is.null(dframe)) stop("No input dataframe.\n")
  if(!is.data.frame(dframe)) stop("Argument 'dframe' is not a dataframe.\n")
  colnames(dframe) <- tolower(colnames(dframe))
  if(!all(c("x","y") %in% colnames(dframe))) stop("Columns 'x' or 'y' not found in input dataframe.\n")
  wells <- as.integer(wells)[1]
  if(is.null(wells) | is.na(wells) | length(wells) == 0) stop("Argument 'wells' is empty.\n")
  if(!is.integer(wells)) stop("Argument 'wells' is not an integer. Set as 24 or 48.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")

  #determine rows and columns
  if(wells == 24)
  {
    prows <- 4
    pcols <- 6
  }

  if(wells == 48)
  {
    prows <- 6
    pcols <- 8
  }

  if(wells == 96)
  {
    prows <- 8
    pcols <- 12
  }

  # creating well layout
  verticalrange <- range(dframe$y)
  horizontalrange <- range(dframe$x)
  vertpos <- cumsum(c(verticalrange[1],rep((diff(verticalrange)/prows),prows)))
  horizpos <- cumsum(c(horizontalrange[1],rep((diff(horizontalrange)/pcols),pcols)))

  h1 <- vector()
  h2 <- vector()
  v1 <- vector()
  v2 <- vector()
  for(i in 1:prows)
  {
    for(j in 1:pcols)
    {
      h1 <- c(h1,horizpos[j])
      h2 <- c(h2,horizpos[j+1])
      v1 <- c(v1,vertpos[i])
      v2 <- c(v2,vertpos[i+1])
    }
  }

  wellsdf <- data.frame(well = factor(1:wells),h1 = h1,h2 = h2,v1 = v1,v2 = v2, stringsAsFactors = FALSE)

  #mid point of blocks
  wellsdf$x <- wellsdf$h1+((wellsdf$h2-wellsdf$h1)/2)
  wellsdf$y <- wellsdf$v1+((wellsdf$v2-wellsdf$v1)/2)

  #top left label positions
  wellsdf$labx <- wellsdf$h1+(0.08*(diff(horizontalrange)/prows))
  wellsdf$laby <- wellsdf$v1+(0.15*(diff(verticalrange)/pcols))

  #bottom label positions
  if(wells == 24) wellsdf$laby1 <- wellsdf$y + 47
  if(wells == 48) wellsdf$laby1 <- wellsdf$y + 30
  if(wells == 96) wellsdf$laby1 <- wellsdf$y + 20

  return(wellsdf)
}

#-------------------------------------------------------------------------------

# FUNCTION SPOTCOVERAGE
#' Computes coverage (convex hull) of xy spot data.
#' @description Computes coverage (convex hull) of xy spot data.
#' @param files A character or vector of paths or filenames. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param filenamediscard A character for part of the filename to be removed. '-Tracks.txt' is removed by default.
#' @param mm A numeric indicating number of pixels in 1 mm.
#' @param exportdata A logical indicating if data table must be exported as a text file to the working directory.
#' @param exportplot A logical if results must be plotted as a figure and exported to working directory.
#' @param quiet A logical indicating if messages should be printed to console during the run. If \code{FALSE}, all output to console is killed except progress bar.
#' @return Returns a dataframe with columns: row, plate names, wells, coverage (convex hull area) in  pixel square and millimetre square.
#' If \code{exportdata=TRUE}, the dataframe is export as a text file in the working directory for each input file.
#' @details
#' Note that the 'xx-EdgeMarks.txt' file must be present in the working directory for the plot exports to work.
#' \strong{files}\cr
#' A character or vector of paths to files. On windows, use \code{choose.files()} for interactive selection.
#' The input files must be tab-delimited decimal as dot (.) text files. The file
#' must contain a minimum of 3 columns named x, y and well. x is a numeric
#' indicating x coordinate and y is a numeric indicating y coordinate of each spot. well indicates well for each spot.
#' Any extra columns are not used. Use output from \code{linkFrames()}.\cr
#' A text file with edge marks (filename-EdgeMarks.txt) must be present in the working directory for plotting. The edge marks file is exported by function \code{assignWells()}.
#' @export

spotCoverage <- function(files = NA, wells = 24, filenamediscard = "-Tracks.txt", mm = 5.4,
                         exportdata = TRUE, exportplot = TRUE, quiet = FALSE)
{
  if(!quiet) currtime <- Sys.time()

  # internal function to compute polygon area
  polyarea <- function (x, y)
  {
    if (length(x) == 0 && length(y) == 0)
      return(0)
    if (!(is.numeric(x) || is.complex(x)) || !(is.numeric(y) || is.complex(y)))
      stop("Arguments 'x' and 'y' must be real or complex.")
    if (is.null(dim(x)))
      x <- matrix(x, length(x), 1)
    if (is.null(dim(y)))
      y <- matrix(y, length(y), 1)
    if (any(dim(x) != dim(y)))
      stop("Matrices 'x' and 'y' must be of same size.")
    n <- nrow(x)
    m <- ncol(x)
    z <- numeric(m)
    for (i in 1:m) {
      xi <- x[, i]
      yi <- y[, i]
      p1 <- sum(xi[1:(n - 1)] * yi[2:n]) + xi[n] * yi[1]
      p2 <- sum(xi[2:n] * yi[1:(n - 1)]) + xi[1] * yi[n]
      z[i] <- 0.5 * (p1 - p2)
    }
    return(z)
  }

  # Argument checks
  if(any(is.null(files)) | any(is.na(files)) | length(files)==0) stop("No input files.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.numeric(mm)) stop("Argument 'mm' not set correctly. Use a numeric value.\n")
  if(!is.character(filenamediscard)) stop("Argument 'filenamediscard' not set correctly. Assign a character.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(quiet)) stop("Argument 'quiet' not set correctly. Set as TRUE or FALSE.\n")

  #looping over selected files
  fnames <- sub(filenamediscard,"",basename(files))
  flen <- length(files)
  dframelist <- vector("list",length = flen)
  fileloop=1
  for(fileloop in 1:flen)
  {
    if(!quiet) cat(paste0("Computing convex hull. File ",fileloop," of ",flen,".\n"))
    if(!quiet) cat(paste0("Current file: ",files[fileloop],".\n"))


    dframe <- read.delim(files[fileloop],header = T,stringsAsFactors = F)
    colnames(dframe) <- tolower(colnames(dframe))
    if(!all(c("x","y","well") %in% colnames(dframe))) stop("Columns 'x', 'y' or 'well' not found in input file.\n")
    #extradf <- subset(dframe, select = -c(x,y,slice))
    #dframe <- dframe[,c("x","y","slice")]

    if(exportplot) covlisthulls <- vector("list",length=wells)
    covareapxsq <- vector(length=wells)
    covareammsq <- vector(length=wells)
    pb <- txtProgressBar(min = 0, max = wells, style = 3)
    coverageloop = 1
    for(coverageloop in 1:wells)
    {
      setTxtProgressBar(pb, coverageloop)
      covwell <- subset(dframe,dframe$well == coverageloop)[,c("x","y","well")]

      if(nrow(covwell) > 2)
      {
        if(exportplot) covlisthulls[[coverageloop]] <- covwell[chull(x = covwell$x,y = covwell$y),]
        areatemp <- round(abs(polyarea(covlisthulls[[coverageloop]]$x,covlisthulls[[coverageloop]]$y)),0)
        covareapxsq[coverageloop] <- areatemp
        covareammsq[coverageloop] <- areatemp/(mm^2)
      }else{
        if(exportplot) covlisthulls[[coverageloop]] <- data.frame(x=NA,y=NA,well=coverageloop,stringsAsFactors=FALSE)
        covareapxsq[coverageloop] <- NA
        covareammsq[coverageloop] <- NA
      }

    }
    close(pb)
    if(exportplot) covdfhulls <- plyr::rbind.fill(covlisthulls)

    covdfarea <- data.frame(row = 1:wells,plate = rep(fnames[fileloop],wells),well = 1:wells,
                            coverage_pxsq = covareapxsq, coverage_mmsq = covareammsq,stringsAsFactors = FALSE)

    if(exportdata)
    {
      write.table(x = covdfarea, file = paste0(fnames[fileloop],"-SpotCoverage.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fnames[fileloop],"-SpotCoverage.txt exported.\n"))
    }

    if(exportplot)
    {
      if(paste0(fnames[fileloop],"-EdgeMarks.txt") %in% list.files())
      {
        edgedf <- read.delim(paste0(fnames[fileloop],"-EdgeMarks.txt"),header=T,stringsAsFactors = FALSE)
        wellsdf <- gridWells(dframe = edgedf, wells = wells)
        covdata <- base::merge(covdfarea,covdfhulls,by="well")
        covwells <- base::merge(wellsdf,covdfarea,by = "well")

        if(wells == 24)
        {
          txtsz <- 1.5
          txtszbp <- 1.8
        }

        if(wells == 48)
        {
          txtsz <- 1.3
          txtszbp <- 1.5
        }

        if(wells == 96)
        {
          txtsz <- 1.1
          txtszbp <- 1.2
        }

        if(wells == 24) cent <- 7
        if(wells == 48) cent <- 4
        if(wells == 96) cent <- 2

        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          geom_rect(data = wellsdf,aes(xmin = h1,xmax = h2,ymin = v1,ymax = v2),colour = "grey90",fill = "white",size = 0.3,na.rm = TRUE)+
          geom_point(data = dframe,aes(x,y),colour="lightgrey",size = 0.25,alpha=0.85,na.rm = TRUE)+
          geom_polygon(data = covdfhulls,aes(x = x,y = y,colour = factor(well),group = factor(well)),fill=NA,size = 0.25,na.rm = TRUE)+
          geom_point(data = covdfhulls,aes(x = x,y = y,colour = factor(well)),size = 0.8,alpha=0.9,na.rm = TRUE)+
          geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)+
          theme_bw()+
          scale_x_continuous(expand = c(0,0)) +
          scale_y_reverse(expand = c(0,0))+
          labs(title = paste0(wells, " Well Plate: ",fnames[fileloop],"  |  Coverage"))+
          coord_fixed()+
          theme_bw(base_size = 5)+
          theme(legend.position = "none",axis.ticks = element_blank(),legend.text = element_blank(),
                axis.text = element_blank(),axis.title = element_blank(),legend.title = element_blank(),
                plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 4.5,face="bold"),panel.border = element_blank(),
                panel.background = element_blank(),panel.grid = element_blank())
        ggsave(filename = paste0(fnames[fileloop],"-SpotCoverage.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fnames[fileloop],"-SpotCoverage.png exported.\n"))
      }else{
        warning(paste0("Coverage plot not exported since file ",paste0(fnames[fileloop],"EdgeMarks.txt") ," not found in working directory.\n"))
      }
    }
  }
  if(!quiet) cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
  return(covdfarea)
}

#-------------------------------------------------------------------------------

# FUNCTION SPOTALPHAHULL
#' Computes alphahull (concave hull) of xy spot data.
#' @description Computes alphahull (concave hull) of xy spot data.
#' @param files A character or vector of paths or filenames. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param filenamediscard A character for part of the filename to be removed. '-Tracks.txt' is removed by default.
#' @param mm A numeric indicating number of pixels in 1 mm.
#' @param alphavalue A numeric indicating alpha value for alphahull. See details.
#' @param exportdata A logical indicating if data table must be exported as a text file to the working directory.
#' @param exportplot A logical if results must be plotted as a figure and exported to working directory.
#' @param quiet A logical indicating if messages should be printed to console during the run. If \code{FALSE}, all output to console is killed except progress bar.
#' @return Returns a dataframe with columns: row, plate names, wells, alphahull (concave hull area) in  pixel square and millimetre square.
#' If \code{exportdata=TRUE}, the dataframe is export as a text file in the working directory for each input file.
#' @details
#' Note that the 'xx-EdgeMarks.txt' file must be present in the working directory for the plot exports to work.
#' \strong{files}\cr
#' A character or vector of paths to files. On windows, use \code{choose.files()} for interactive selection.
#' The input files must be tab-delimited decimal as dot (.) text files. The file
#' must contain a minimum of 3 columns named x, y and well. x is a numeric
#' indicating x coordinate and y is a numeric indicating y coordinate of each spot. well indicates well for each spot.
#' Any extra columns are not used. Use output from \code{linkFrames()}.\cr
#' A text file with edge marks (filename-EdgeMarks.txt) must be present in the working directory for plotting. The edge marks file is exported by function \code{assignWells()}.
#' \strong{alphavalue}\cr
#' A default value of 4 is used. Smaller values produce more gaps in the spot cloud.
#' @export
#' @import alphahull

spotAlphahull <- function(files = NA, wells = 24, filenamediscard = "-Tracks.txt", mm = 5.4, alphavalue = 4,
                          exportdata = TRUE, exportplot = TRUE, quiet = FALSE)
{
  if(!quiet) currtime <- Sys.time()

  # Argument checks
  if(any(is.null(files)) | any(is.na(files)) | length(files)==0) stop("No input files.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.numeric(mm)) stop("Argument 'mm' not set correctly. Use a numeric value.\n")
  if(!is.numeric(alphavalue)) stop("Argument 'alphavalue' not set correctly. Use a numeric value.\n")
  if(!is.character(filenamediscard)) stop("Argument 'filenamediscard' not set correctly. Assign a character.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(quiet)) stop("Argument 'quiet' not set correctly. Set as TRUE or FALSE.\n")

  #looping over selected files
  fnames <- sub(filenamediscard,"",basename(files))
  flen <- length(files)
  dframelist <- vector("list",length = flen)
  fileloop=1
  for(fileloop in 1:flen)
  {
    if(!quiet) cat(paste0("Computing alpha hull. File ",fileloop," of ",flen,".\n"))
    if(!quiet) cat(paste0("Current file: ",files[fileloop],".\n"))

    dframe <- read.delim(files[fileloop],header = T,stringsAsFactors = F)
    colnames(dframe) <- tolower(colnames(dframe))
    if(!all(c("x","y","well") %in% colnames(dframe))) stop("Columns 'x', 'y' or 'well' not found in input file.\n")

    ahvecpxsq <- vector(length=wells)
    if(exportplot) ahverlist <- vector("list",length=wells)
    pb <- txtProgressBar(min = 0, max = wells, style = 3)
    for(ahloop in 1:wells)
    {
      setTxtProgressBar(pb, ahloop)
      #ahwell <- ahlist[[ahloop]][,c("x","y")]
      ahwell <- subset(dframe,dframe$well == ahloop)[,c("x","y")]

      if(nrow(ahwell) > 2)
      {
        #remove duplicates
        ahwell <- ahwell[-which(duplicated(round(ahwell,1))),]
        #jitter
        ahwell$x <- jitter(ahwell$x,0.1,0.1)
        ahwell$y <- jitter(ahwell$y,0.1,0.1)

        #collinearity check
        s1 <- (ahwell$y[2]-ahwell$y[1])/(ahwell$x[2]-ahwell$x[1])
        s2 <- (ahwell$y[3]-ahwell$y[2])/(ahwell$x[3]-ahwell$x[2])
        if(!is.na(s1) && !is.na(s2))
        {
          if(s1 == s2) ahwell <- ahwell[sample(1:nrow(ahwell)),]
        }
      }

      if(nrow(ahwell) > 2)
      {
        ah1 <- alphahull::ahull(x=ahwell$x,ahwell$y,alpha=alphavalue)
        ahvecpxsq[ahloop] <- alphahull::areaahull(ah1)
        if(exportplot) ahverlist[[ahloop]] <- as.data.frame(ah1$ashape.obj$edges)[,c("x1","y1","x2","y2")]
        if(exportplot) ahverlist[[ahloop]]$well <- rep(ahloop,nrow(ahverlist[[ahloop]]))
        #segments(x0 = t1$x1,y0 = t1$y1,x1 = t1$x2,y1=t1$y2,col="red")
      }

      if(nrow(ahwell) < 3)
      {
        ahvecpxsq[ahloop] <- NA
        if(exportplot) ahverlist[[ahloop]] <- data.frame(x1=NA, y1=NA, x2=NA, y2=NA, well=ahloop, stringsAsFactors = FALSE)
      }
    }
    close(pb)
    if(exportplot) ahdf <- plyr::rbind.fill(ahverlist)

    ahdfarea <- data.frame(row = 1:wells,plate = rep(fnames[fileloop],wells),well = 1:wells,
                           alphahull_pxsq = round(ahvecpxsq,2), alphahull_mmsq = round(ahvecpxsq/(mm^2),2),stringsAsFactors = FALSE)

    if(exportdata)
    {
      write.table(x = ahdfarea, file = paste0(fnames[fileloop],"-SpotAlphahull.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fnames[fileloop],"-SpotAlphahull.txt exported.\n"))
    }

    if(exportplot)
    {
      if(paste0(fnames[fileloop],"-EdgeMarks.txt") %in% list.files())
      {
        edgedf <- read.delim(paste0(fnames[fileloop],"-EdgeMarks.txt"),header=T,stringsAsFactors = FALSE)
        wellsdf <- gridWells(dframe = edgedf, wells = wells)

        if(wells == 24)
        {
          txtsz <- 1.5
          txtszbp <- 1.8
        }

        if(wells == 48)
        {
          txtsz <- 1.3
          txtszbp <- 1.5
        }

        if(wells == 96)
        {
          txtsz <- 1.1
          txtszbp <- 1.2
        }

        if(wells == 24) cent <- 7
        if(wells == 48) cent <- 4
        if(wells == 96) cent <- 2

        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          geom_point(data = dframe,aes(x,y),colour="lightgrey",size = 0.25,alpha=0.85,na.rm = TRUE)+
          geom_segment(data=ahdf,aes(x=x1,xend=x2,y=y1,yend=y2,colour=factor(well)),size=0.2,na.rm = TRUE)+
          geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)+
          scale_x_continuous(expand = c(0,0))+
          scale_y_reverse(expand = c(0,0))+
          theme_bw(base_size = 5)+
          labs(title = paste0(wells, " Well Plate: ",fnames[fileloop],"  |  Alpha Hulls at alpha ",alphavalue))+
          coord_fixed()+
          theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 4.5, face="bold"),legend.position = "none",
                axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
                axis.line = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
                panel.border = element_blank())
        ggsave(filename = paste0(fnames[fileloop],"-SpotAlphahull.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fnames[fileloop],"-SpotAlphahull.png exported.\n"))
      }else{
        warning(paste0("Alphahull plot not exported since file ",paste0(fnames[fileloop],"EdgeMarks.txt") ," not found in working directory.\n"))
      }
    }

  }
  if(!quiet) cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
  return(ahdfarea)
}

#-------------------------------------------------------------------------------

# FUNCTION SPOTMSD
#' Computes minimum spanning distance of xy spot data.
#' @description Computes minimum spanning distance of xy spot data.
#' @param files A character or vector of paths or filenames. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param filenamediscard A character for part of the filename to be removed. '-Tracks.txt' is removed by default.
#' @param mm A numeric indicating number of pixels in 1 mm.
#' @param exportdata A logical indicating if data table must be exported as a text file to the working directory.
#' @param exportplot A logical if results must be plotted as a figure and exported to working directory.
#' @param quiet A logical indicating if messages should be printed to console during the run. If \code{FALSE}, all output to console is killed except progress bar.
#' @return Returns a dataframe with columns: row, plate names, wells, msd in  pixels and millimetre.
#' If \code{exportdata=TRUE}, the dataframe is export as a text file in the working directory for each input file.
#' @details
#' Note that the 'xx-EdgeMarks.txt' file must be present in the working directory for the plot exports to work.
#' \strong{files}\cr
#' A character or vector of paths to files. On windows, use \code{choose.files()} for interactive selection.
#' The input files must be tab-delimited decimal as dot (.) text files. The file
#' must contain a minimum of 3 columns named x, y and well. x is a numeric
#' indicating x coordinate and y is a numeric indicating y coordinate of each spot. well indicates well assignment for each spot.
#' Any extra columns are not used. Use output from \code{linkFrames()}.\cr
#' A text file with edge marks (filename-EdgeMarks.txt) must be present in the working directory for plotting. The edge marks file is exported by function \code{assignWells()}.
#' @export
#' @import fossil
#'
spotMsd <- function(files = NA, wells = 24, filenamediscard = "-Tracks.txt", mm = 5.4,
                    exportdata = TRUE, exportplot = TRUE, quiet = FALSE)
{
  if(!quiet) currtime <- Sys.time()

  # Argument checks
  if(any(is.null(files)) | any(is.na(files)) | length(files)==0) stop("No input files.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.numeric(mm)) stop("Argument 'mm' not set correctly. Use a numeric value.\n")
  if(!is.character(filenamediscard)) stop("Argument 'filenamediscard' not set correctly. Assign a character.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(quiet)) stop("Argument 'quiet' not set correctly. Set as TRUE or FALSE.\n")

  #looping over selected files
  fnames <- sub(filenamediscard,"",basename(files))
  flen <- length(files)
  dframelist <- vector("list",length = flen)
  fileloop=1
  for(fileloop in 1:flen)
  {
    if(!quiet) cat(paste0("Computing minimum spanning tree. File ",fileloop," of ",flen,".\n"))
    if(!quiet) cat(paste0("Current file: ",files[fileloop],".\n"))

    dframe <- read.delim(files[fileloop],header = T,stringsAsFactors = F)
    colnames(dframe) <- tolower(colnames(dframe))
    if(!all(c("x","y","well") %in% colnames(dframe))) stop("Columns 'x', 'y' or 'well' not found in input file.\n")

    msdvec <- vector(length=wells)
    msdlinelist <- vector("list",length=wells)
    pb <- txtProgressBar(min = 0, max = wells, style = 3)
    msdloop=1
    for(msdloop in 1:wells)
    {
      setTxtProgressBar(pb, msdloop)
      msdwell <- subset(dframe,dframe$well == msdloop)[,c("x","y")]
      #msdwell <- msdlist[[msdloop]][,c("x","y")]
      #remove duplicates
      if(nrow(msdwell)>1) msdwell <- msdwell[-which(duplicated(round(msdwell,1))),]

      if(nrow(msdwell) > 2)
      {
        #jitter
        msdwell$x <- jitter(msdwell$x,0.1,0.1)
        msdwell$y <- jitter(msdwell$y,0.1,0.1)

        #compute msn from euclidean distance
        msdobj <- fossil::dino.mst(stats::dist(msdwell))
        msdobj1 <- as.data.frame(as.table(t(msdobj)),stringsAsFactors=FALSE)
        msdobj2 <- subset(msdobj1,msdobj1$Freq!=0)
        msdobj3 <- cbind(msdwell[msdobj2$Var1,],msdwell[msdobj2$Var2,])
        colnames(msdobj3) <- c("x0","y0","x1","y1")
        msdobj3$well <- rep(msdloop,nrow(msdobj3))
        #calculate distance between point pairs
        msdobj3$msd_px <- sqrt(((msdobj3$x1-msdobj3$x0)^2)+((msdobj3$y1-msdobj3$y0)^2))
        msdlinelist[[msdloop]] <- msdobj3
        #tot distance of all point pairs
        msdvec[msdloop] <- sum(msdobj3$msd_px)
      }else{
        msdvec[msdloop] <- NA
        msdlinelist[[msdloop]] <- data.frame(x0=NA,y0=NA,x1=NA,y1=NA,well=msdloop,msd_px=NA,stringsAsFactors=FALSE)
      }

    }
    close(pb)
    msddf <- plyr::rbind.fill(msdlinelist)

    msddf1 <- data.frame(row = 1:wells,plate = rep(fnames[fileloop],wells),well = 1:wells,
                         msd_px = round(msdvec,2), msd_mm = round(msdvec/mm,2),stringsAsFactors = FALSE)

    if(exportdata)
    {
      write.table(x = msddf1, file = paste0(fnames[fileloop],"-SpotMsd.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fnames[fileloop],"-SpotMsd.txt exported.\n"))
    }

    if(exportplot)
    {
      if(paste0(fnames[fileloop],"-EdgeMarks.txt") %in% list.files())
      {
        edgedf <- read.delim(paste0(fnames[fileloop],"-EdgeMarks.txt"),header=T,stringsAsFactors = FALSE)
        wellsdf <- gridWells(dframe = edgedf, wells = wells)

        if(wells == 24)
        {
          txtsz <- 1.5
          txtszbp <- 1.8
        }

        if(wells == 48)
        {
          txtsz <- 1.3
          txtszbp <- 1.5
        }

        if(wells == 96)
        {
          txtsz <- 1.1
          txtszbp <- 1.2
        }

        if(wells == 24) cent <- 7
        if(wells == 48) cent <- 4
        if(wells == 96) cent <- 2

        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          geom_segment(data = msddf,aes(x = x0,xend = x1,y = y0,yend = y1,group = well,colour=factor(well)),size=0.2,na.rm = TRUE)+
          geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)+
          scale_x_continuous(expand = c(0,0))+
          scale_y_reverse(expand = c(0,0))+
          theme_bw(base_size = 5)+
          labs(title = paste0(wells, " Well Plate: ",fnames[fileloop],"  |  Minimum Spanning Distance"))+
          coord_fixed()+
          theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 4.5, face="bold"),legend.position = "none",
                axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
                axis.line = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
                panel.border = element_blank())
        ggsave(filename = paste0(fnames[fileloop],"-SpotMsd.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fnames[fileloop],"-SpotMsd.png exported.\n"))
      }else{
        warning(paste0("Msd plot not exported since file ",paste0(fnames[fileloop],"EdgeMarks.txt") ," not found in working directory.\n"))
      }
    }
  }
  if(!quiet) cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
  return(msddf1)
}

#-------------------------------------------------------------------------------

# FUNCTION SPOTDENSITY
#' Computes 2D kernal density of xy spot data.
#' @description Computes 2D kernal density of xy spot data.
#' @param files A character or vector of paths or filenames. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param filenamediscard A character for part of the filename to be removed. '-Tracks.txt' is removed by default.
#' @param mm A numeric indicating number of pixels in 1 mm.
#' @param exportdata A logical indicating if data table must be exported as a text file to the working directory.
#' @param exportplot A logical if results must be plotted as a figure and exported to working directory.
#' @param quiet A logical indicating if messages should be printed to console during the run. If \code{FALSE}, all output to console is killed except progress bar.
#' @return Returns Nothing for now. If \code{exportplot=T}, a figure with plotted with spot densities.
#' @details
#' Note that the 'xx-EdgeMarks.txt' file must be present in the working directory for the plot exports to work.
#' \strong{files}\cr
#' A character or vector of paths to files. On windows, use \code{choose.files()} for interactive selection.
#' The input files must be tab-delimited decimal as dot (.) text files. The file
#' must contain a minimum of 3 columns named x, y and well. x is a numeric
#' indicating x coordinate and y is a numeric indicating y coordinate of each spot. well indicates well assignment for each spot.
#' Any extra columns are not used. Use output from \code{linkFrames()}.\cr
#' A text file with edge marks (filename-EdgeMarks.txt) must be present in the working directory for plotting. The edge marks file is exported by function \code{assignWells()}.
#' @export
#' @import fossil
#'
spotDensity <- function(files = NA, wells = 24, filenamediscard = "-Tracks.txt", mm = 5.4,
                        exportdata = TRUE, exportplot = TRUE, quiet = FALSE)
{
  if(!quiet) currtime <- Sys.time()

  # Argument checks
  if(any(is.null(files)) | any(is.na(files)) | length(files)==0) stop("No input files.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.numeric(mm)) stop("Argument 'mm' not set correctly. Use a numeric value.\n")
  if(!is.character(filenamediscard)) stop("Argument 'filenamediscard' not set correctly. Assign a character.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(quiet)) stop("Argument 'quiet' not set correctly. Set as TRUE or FALSE.\n")

  #looping over selected files
  fnames <- sub(filenamediscard,"",basename(files))
  flen <- length(files)
  dframelist <- vector("list",length = flen)
  fileloop=1
  for(fileloop in 1:flen)
  {
    if(!quiet) cat(paste0("Computing spot density. File ",fileloop," of ",flen,".\n"))
    if(!quiet) cat(paste0("Current file: ",files[fileloop],".\n"))

    dframe <- read.delim(files[fileloop],header = T,stringsAsFactors = F)
    colnames(dframe) <- tolower(colnames(dframe))
    if(!all(c("x","y","well") %in% colnames(dframe))) stop("Columns 'x', 'y' or 'well' not found in input file.\n")

    kdlist <- vector("list",length=wells)
    pb <- txtProgressBar(min = 0, max = wells, style = 3)
    kdloop = 1
    for(kdloop in 1:wells)
    {
      setTxtProgressBar(pb, kdloop)
      kdwell <- subset(dframe,dframe$well == kdloop)[,c("x","y")]
      #kdwell <- covlist[[kdloop]][,c("x","y","well")]

      if(nrow(kdwell) > 2)
      {
        #declare constants
        nbin <- 128
        cols <- c("midnightblue","#00FEFF","#45FE4F","#FCFF00","#FF9400", "#FF3100")
        #dim colours
        #cols<-c("#324AA2","#77C9D4","#71C637","#F7DC13","#F97419","#E83C14")
        #use densCols() output to get density at each point
        kdwell$col <- as.character(grDevices::densCols(kdwell$x,kdwell$y, nbin=nbin,colramp=colorRampPalette(cols)))
        kdwell$well <- rep(kdloop,nrow(kdwell))
        kdlist[[kdloop]] <- kdwell
      }else{
        kdlist[[kdloop]] <- data.frame(x=NA,y=NA,kd=NA,col=NA,well=kdloop,stringsAsFactors = FALSE)
      }

    }
    close(pb)
    kddf <- plyr::rbind.fill(kdlist)

#     if(exportdata)
#     {
#       write.table(x = , file = paste0(fnames[fileloop],"-SpotKernalDensity.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
#       if(!quiet) cat(paste0(fnames[fileloop],"-SpotKernalDensity.txt exported.\n"))
#     }

    if(exportplot)
    {
      if(paste0(fnames[fileloop],"-EdgeMarks.txt") %in% list.files())
      {
        edgedf <- read.delim(paste0(fnames[fileloop],"-EdgeMarks.txt"),header=T,stringsAsFactors = FALSE)
        wellsdf <- gridWells(dframe = edgedf, wells = wells)

        if(wells == 24)
        {
          txtsz <- 1.5
          txtszbp <- 1.8
        }

        if(wells == 48)
        {
          txtsz <- 1.3
          txtszbp <- 1.5
        }

        if(wells == 96)
        {
          txtsz <- 1.1
          txtszbp <- 1.2
        }

        if(wells == 24) cent <- 7
        if(wells == 48) cent <- 4
        if(wells == 96) cent <- 2

        colvals <- as.character(levels(factor(kddf$col,ordered=F)))

        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          geom_point(data = kddf,aes(x,y,colour=factor(col)),size = 0.3,alpha=0.9,na.rm = TRUE)+
          scale_colour_manual(values=colvals)+
          geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)+
          scale_x_continuous(expand = c(0,0))+
          scale_y_reverse(expand = c(0,0))+
          theme_bw(base_size = 5)+
          labs(title = paste0(wells, " Well Plate: ",fnames[fileloop],"  |  2D Kernal Density"))+
          coord_fixed()+
          theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 4.5, face = "bold"),legend.position = "none",
                axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
                axis.line = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
                panel.border = element_blank())
        ggsave(filename = paste0(fnames[fileloop],"-SpotDensity.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fnames[fileloop],"-SpotDensity.png exported.\n"))
      }else{
        warning(paste0("Density plot not exported since file ",paste0(fnames[fileloop],"EdgeMarks.txt") ," not found in working directory.\n"))
      }
    }
  }
  if(!quiet) cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
  #return(ahdfarea)
}

#-------------------------------------------------------------------------------

# FUNCTION SPOTRATIO
#' Computes larval presence in center vs edges of the wells.
#' @description Counts spots in the center vs edges of the wells and computes a ratio.
#' @param files A character or vector of paths or filenames. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param filenamediscard A character for part of the filename to be removed. '-Tracks.txt' is removed by default.
#' @param mm A numeric indicating number of pixels in 1 mm.
#' @param exportdata A logical indicating if data table must be exported as a text file to the working directory.
#' @param exportplot A logical if results must be plotted as a figure and exported to working directory.
#' @param quiet A logical indicating if messages should be printed to console during the run. If \code{FALSE}, all output to console is killed except progress bar.
#' @return Returns a dataframe with columns: row, plate names, wells, number of outer spots, inner spots and spots ratio.
#' If \code{exportdata=TRUE}, the dataframe is export as a text file in the working directory for each input file.
#' @details
#' Note that the 'xx-EdgeMarks.txt' file must be present in the working directory for the plot exports to work.
#' \strong{files}\cr
#' A character or vector of paths to files. On windows, use \code{choose.files()} for interactive selection.
#' The input files must be tab-delimited decimal as dot (.) text files. The file
#' must contain a minimum of 3 columns named x, y and well. x is a numeric
#' indicating x coordinate and y is a numeric indicating y coordinate of each spot. well indicates well assignment for each spot.
#' Any extra columns are not used. Use output from \code{linkFrames()}.\cr
#' A text file with edge marks (filename-EdgeMarks.txt) must be present in the working directory for plotting. The edge marks file is exported by function \code{assignWells()}.
#' @export
#' @import SDMTools
#'
spotRatio<- function(files = NA, wells = 24, filenamediscard = "-Tracks.txt", mm = 5.4,
                     exportdata = TRUE, exportplot = TRUE, quiet = FALSE)
{
  if(!quiet) currtime <- Sys.time()

  #function circlepos
  circlepos <- function (x, y, radius, nv = 100)
  {
    ymult <- 1
    angle.inc <- 2 * pi/nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (length(radius) < length(x))
      radius <- rep(radius, length.out = length(x))
    if (length(col) < length(radius))
      col <- rep(col, length.out = length(radius))
    for (circle in 1:length(radius)) {
      xv <- cos(angles) * radius[circle] + x[circle]
      yv <- sin(angles) * radius[circle] * ymult + y[circle]
    }
    return(data.frame(x = xv, y = yv, stringsAsFactors = FALSE))
  }

  # Argument checks
  if(any(is.null(files)) | any(is.na(files)) | length(files)==0) stop("No input files.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.numeric(mm)) stop("Argument 'mm' not set correctly. Use a numeric value.\n")
  if(!is.character(filenamediscard)) stop("Argument 'filenamediscard' not set correctly. Assign a character.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(quiet)) stop("Argument 'quiet' not set correctly. Set as TRUE or FALSE.\n")

  #looping over selected files
  fnames <- sub(filenamediscard,"",basename(files))
  flen <- length(files)
  dframelist <- vector("list",length = flen)
  fileloop=1
  for(fileloop in 1:flen)
  {
    if(!quiet) cat(paste0("Computing spot ratio. File ",fileloop," of ",flen,".\n"))
    if(!quiet) cat(paste0("Current file: ",files[fileloop],".\n"))

    if(paste0(fnames[fileloop],"-EdgeMarks.txt") %in% list.files())
    {
      edgedf <- read.delim(paste0(fnames[fileloop],"-EdgeMarks.txt"),header=T,stringsAsFactors = FALSE)
      dframe <- read.delim(files[fileloop],header = T,stringsAsFactors = F)
      colnames(dframe) <- tolower(colnames(dframe))
      if(!all(c("x","y","well") %in% colnames(dframe))) stop("Columns 'x', 'y' or 'well' not found in input file.\n")
      wellsdf <- gridWells(dframe = edgedf, wells = wells)

      outerpointslist <- vector("list",length=wells)
      innerpointslist <- vector("list",length=wells)
      outercirclelist <- vector("list",length=wells)
      innercirclelist <- vector("list",length=wells)
      opvec <- vector(length=wells)
      ipvec <- vector(length=wells)
      pb <- txtProgressBar(min = 0, max = wells, style = 3)
      bsloop = 1
      for(bsloop in 1:wells)
      {
        setTxtProgressBar(pb, bsloop)
        bswell <- subset(dframe,dframe$well == bsloop)[,c("x","y")]

        if(nrow(bswell) > 2)
        {
          #calculate circle polygon for well
          outercircle <- circlepos(wellsdf$x[bsloop],wellsdf$y[bsloop],radius = 44)
          innercircle <- circlepos(wellsdf$x[bsloop],wellsdf$y[bsloop],radius = 26)
          #num of points in large circle
          outerpoints <- SDMTools::pnt.in.poly(bswell,outercircle)
          outerpoints <- subset(outerpoints, outerpoints$pip == 1)
          #num of points in small circle
          innerpoints <- SDMTools::pnt.in.poly(bswell,innercircle)
          innerpoints <- subset(innerpoints, innerpoints$pip == 1)
          #add number of op and ip
          opt <- nrow(outerpoints) - nrow(innerpoints)
          #convert negative values to zero
          if(opt < 0){opvec[bsloop] <- 0}else{opvec[bsloop] <- opt}
          ipvec[bsloop] <- nrow(innerpoints)

          #add wells to points
          outerpoints$well <- rep(bsloop,nrow(outerpoints))
          innerpoints$well <- rep(bsloop,nrow(innerpoints))
          #add wells to circles
          outercircle$well <- rep(bsloop,nrow(outercircle))
          innercircle$well <- rep(bsloop,nrow(innercircle))

          #add points to list
          outerpointslist[[bsloop]] <- outerpoints
          innerpointslist[[bsloop]] <- innerpoints
          #add circle to list
          outercirclelist[[bsloop]] <- outercircle
          innercirclelist[[bsloop]] <- innercircle
        }else{
          opvec[bsloop] <- NA
          ipvec[bsloop] <- NA
          outerpointslist[[bsloop]] <- data.frame(x=NA,y=NA,well=bsloop,stringsAsFactors=FALSE)
          innerpointslist[[bsloop]] <- data.frame(x=NA,y=NA,well=bsloop,stringsAsFactors=FALSE)
          outercirclelist[[bsloop]] <- data.frame(x=NA,y=NA,well=bsloop,stringsAsFactors=FALSE)
          innercirclelist[[bsloop]] <- data.frame(x=NA,y=NA,well=bsloop,stringsAsFactors=FALSE)
        }

      }
      close(pb)
      opdf <- plyr::rbind.fill(outerpointslist)
      ipdf <- plyr::rbind.fill(innerpointslist)
      ocdf <- plyr::rbind.fill(outercirclelist)
      icdf <- plyr::rbind.fill(innercirclelist)
      spotsdf <- data.frame(row = 1:wells,plate = rep(fnames[fileloop],wells),well = 1:wells,
                            spots_outer = opvec, spots_inner = ipvec, spots_ratio = round((opvec/(opvec+ipvec)) - (ipvec/(opvec+ipvec)),3),
                            stringsAsFactors = FALSE)

      if(exportdata)
      {
        write.table(x = spotsdf, file = paste0(fnames[fileloop],"-SpotRatio.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
        cat(paste0(fnames[fileloop],"-SpotRatio.txt exported.\n"))
      }

      if(exportplot)
      {
        #wellsdf$spots_ratio <- spotsdf$spots_ratio
        t1 <- merge(spotsdf,opdf,by="well")

        if(wells == 24)
        {
          txtsz <- 1.5
          txtszbp <- 1.8
        }

        if(wells == 48)
        {
          txtsz <- 1.3
          txtszbp <- 1.5
        }

        if(wells == 96)
        {
          txtsz <- 1.1
          txtszbp <- 1.2
        }

        if(wells == 24) cent <- 7
        if(wells == 48) cent <- 4
        if(wells == 96) cent <- 2

        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          geom_polygon(data = ocdf,aes(x = x,y = y,group = factor(well)),colour="#a2a18f",fill=NA,size = 0.25,alpha = 0.6,na.rm = TRUE)+
          geom_polygon(data = icdf,aes(x = x,y = y,group = factor(well)),colour="#a2a18f",fill=NA,size = 0.25,alpha = 0.6,na.rm = TRUE)+
          geom_point(data = t1,aes(x,y,colour = spots_ratio),size = 0.3,alpha=0.9,na.rm = TRUE)+
          geom_point(data = ipdf,aes(x,y),colour="grey80",size = 0.3,alpha=0.9,na.rm = TRUE)+
          geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)+
          #geom_text(data = wellsdf,aes(x = x,y = laby1,label = round(spots_ratio,3)),size = txtsz,colour = "grey20",alpha = 0.6,fontface = "bold",na.rm = TRUE)+
          scale_x_continuous(expand = c(0,0))+
          scale_y_reverse(expand = c(0,0))+
          theme_bw(base_size = 5)+
          labs(title = paste0(wells, " Well Plate: ",fnames[fileloop],"  |  Spot Ratio"))+
          coord_fixed()+
          theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 4.5, face = "bold"),legend.position = "none",
                axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
                axis.line = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
                panel.border = element_blank())
        ggsave(filename = paste0(fnames[fileloop],"-SpotRatio.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fnames[fileloop],"-SpotRatio.png exported.\n"))
      }

    }else{
      warning(paste0("Spots Ratio not computed since file ",paste0(fnames[fileloop],"EdgeMarks.txt") ," not found in working directory.\n"))
    }
  }
  if(!quiet) cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
  return(spotsdf)
}
