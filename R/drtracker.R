# LARVAE TRACKING AND FEATURE CALCULATION
# Roy Mathew Francis
# v1.0.4
# 27-Oct-2015

# Assumptions
# 24 or 48 well plate
# plate is constant throughout the film. absolutely no movement.
# edges and larvae positions are marked in frame 1
# single larvae per well
# if larval position is missing in a frame, old position is duplicated
# if two positions are available, nearest is chosen.

#load libraries
#library(plyr) #rbind.fill()
#library(RANN) #nn2()
#library(ggplot2)
#library(fossil) #dino.mst
#library(alphahull) #ahull()
#library(reshape2) #dcast()

#temp <- ltrack("test1.txt",msd=F)

# FUNCTION LTRACK
#' Track larval movement in 24 or 48 well assay plates.
#' @description Tracks single larvae in 24 or 48 well assay plates from xy spot data
#' and computes distance and speed. Exports data as text files and generates plots.
#' @param files A character or vector of paths or filenames. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param markededges A logical indicating if corners of plates have been marked. See details.
#' @param fps A numeric indicating framerate of video in number of frames per second.
#' @param mm A numeric indicating number of pixels in 1 mm.
#' @param activitydist A numeric indicating distance in mm after which the larvae is considered active.
#' @param coverage A logical indicating if the coverge must be computed. See details.
#' @param msd A logical indicating if the minimum spanning distance should be computed. See details.
#' @param alphahull A logical indicating if alphahull of point cloud is to be calculated. See details.
#' @param alphavalue A numeric indicating alpha value for alphahull. See details.
#' @param filenamediscard A character for part of the filename to be removed.
#' @param exportdata A logical indicating if data tables must be exported to working directory.
#' @param exportplot A logical if results must be plotted and exported to working directory.
#' @param follow A character indicating if the algorithm must be followed at every step.
#' Options are 'interactive', 'track' or 'none'. See details.
#' @param centerwell A logical indicating if the well number must be printed at center of the well on plots.
#' @return Returns a list containing two components; 'RawData' and 'Features'.\cr
#' The 'RawData' is a dataframe for one input file or a list for 2 or more input files.
#' The 'RawData' contains xy coordinates, slice, well, id and linktype per spot. The
#' id links the spots as tracks. The linktype shows the method used to link a spot.
#' 'Single' means that there was only one spot in the next frame. This may be the correct
#' spot or a stray spot. 'Nearest' means there was more than one spot in the next frame
#' and the nearest spot was chosen. 'Duplicate' means that no spot was present in the next frame
#' and the previous spot was duplicated in next frame.\cr
#' 'Features' contains plate name, id, well, total distance in pixels and mm, mean speed in pixels per sec (pps) and
#' in mm per second (mmps), max speed in pixels per sec (pps) and mm per sec (mmps), activity (proportion of time), duration of
#' the whole sequence in frames (fr) and seconds (sec), framerate (fps) and calibration,
#' number of pixels in one mm (mm). If \code{coverage=TRUE}, then coverage_pxsq is added.
#' If \code{msd=TRUE}, then msd_px and msd_mm are added. If \code{alphahull=TRUE}, then alphahull_pxsq is added.\cr
#' If \code{exportdata = T}, then two text files are exported: Tracks and Features.
#' If more than one file was selected, a Combined-Features file is also exported.
#' If \code{exportplot = T}, then 3-7 figures are exported: edgespots (if edge spots are used),
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
#' COmputed using \code{spantree} from package \code{vegan}.
#' \strong{alphahull}\cr
#' The alphahull of point cloud based on the alphavalue. The function \code{ahull} from
#' package \code{alphahull}.
#' \strong{filenamediscard}\cr
#' The file name of the input file is used on plots and text files for identification.
#' Part of the filename to be removed can be indicated here. '.txt' is removed by default.\cr
#' \strong{follow}\cr
#' In 'interactive' mode, a plot is shown at every frame and waits for user input.
#' In 'track' mode, the track path creation is shown in real-time.
#' @import plyr
#' @import RANN
#' @import alphahull
#' @import fossil
#' @export
#'
ltrack <- function(files = NULL, wells = 24, markededges = TRUE, fps = 25, mm = 5.4, activitydist = 5,
                   coverage = TRUE,msd = TRUE, alphahull = TRUE, alphavalue = 4,
                   filenamediscard = ".txt", exportdata = TRUE, exportplot = TRUE,
                   follow = "none", centerwell = FALSE)
{
  currtime <- Sys.time()
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

  print(files)
  # Argument checks
  if(any(is.null(files)) | any(is.na(files)) | length(files)==0) stop("No input files.\n")
  if((wells != 24) && (wells != 48)) stop("Set 'wells' as 24 or 48.\n")
  if(!is.logical(markededges)) stop("Argument 'markededges' is not set correctly. Set as TRUE or FALSE.\n")
  if(!is.numeric(fps)) stop("Argument 'fps' not set correctly. Assign a number.\n")
  if(!is.numeric(mm)) stop("Argument 'mm' not set correctly. Assign a number.\n")
  if(!is.character(filenamediscard)) stop("Argument 'filenamediscard' not set correctly. Assign a character.\n")
  if(!is.logical(exportdata)) stop("Argument 'exportdata' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(exportplot)) stop("Argument 'exportplot' not set correctly. Set as TRUE or FALSE.\n")
  if(!is.logical(coverage)) stop("Argument 'coverage' not set correctly. Set as TRUE or FALSE.\n")

  if((follow != "none") && (follow != "interactive") && (follow != "track")) stop("Set argument 'follow' as 'none', 'interactive' or 'track'.\n")

  cat(paste0("Note values in use; Wells: ",wells,", Framerate: ",fps," fps, 1 mm = ",mm," pixels.\n"))

  #looping over selected files
  flen <- length(files)
  distspeedlist <- vector("list",length = flen)
  fileloop=1
  for(fileloop in 1:flen)
  {
    cat(paste0("Processing file ",fileloop," of ",flen," files.\n"))
    cat(paste0("Current file: ",files[fileloop],".\n"))

    fname <- sub(filenamediscard,"",basename(files[fileloop]))
    dframe <- read.delim(files[fileloop],header = T,stringsAsFactors = F)
    colnames(dframe) <- tolower(colnames(dframe))
    if(!all(c("x","y","slice") %in% colnames(dframe))) stop("Columns 'x', 'y' or 'slice' not found in input file.\n")
    if(all(c("id","linktype","frame") %in% colnames(dframe))) warning("Input file may be incorrect.\n")
    dframe <- dframe[,c("x","y","slice")]
    #round all coordinates
    dframe$x <- round(dframe$x,2)
    dframe$y <- round(dframe$y,2)

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

    rm(vertpos,horizpos,h1,h2,v1,v2,i,j)

    #-----------------------------------------------------------------------------

    # handling marked edges
    if(markededges)
    {
      edges <- chull(x = dframe$x,y = dframe$y)
      edgedf <- dframe[edges,]
      # remove edge spots
      dframe <- dframe[-edges,]
      cat(paste0(length(edges)," edge spots removed.\n"))

      # plot edge spots
      if(exportplot)
      {
        p <- ggplot(edgedf,aes(x,y))+
          geom_point(shape = 1,size = 3,na.rm = TRUE)+
          ggtitle(paste0(wells," Well Plate: ",fname,"  |  ",length(edges)," EdgeSpots"))+
          theme_bw(base_size = 5)+
          theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 5),
                axis.title = element_blank(),panel.border = element_blank(),
                axis.text = element_text(colour = "grey40"),
                axis.ticks = element_line(colour = "grey40",size = 0.3))
        ggsave(filename = paste0(fname,"-0-EdgeSpots.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fname,"-0-EdgeSpots.png exported.\n"))
      }

      rm(p,edges)
    }

    #-----------------------------------------------------------------------------

    # compute well positions using nn2
    cat(paste0("Computing well positions.\n"))
    nrowdf <- nrow(dframe)
    cat(paste0("Assigning ",nrowdf," spots to ",wells," wells...\n"))
    pb <- txtProgressBar(min = 0, max = nrowdf, style = 3)
    wellvec <- vector(length = nrowdf)

    for(i in 1:nrowdf)
    {
      wellvec[i] <- as.numeric(RANN::nn2(wellsdf[,c("x","y")],dframe[i,c("x","y")],k = 1)$nn.idx)
      setTxtProgressBar(pb, i)
    }
    dframe$well <- as.numeric(as.character((wellvec)))
    close(pb)
    rm(pb,wellvec,nrowdf,i)

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
        geom_point(data = dframe,aes(x,y,colour = factor(well)),size = 0.3,na.rm = TRUE)
        if(centerwell) p <- p + geom_text(data = wellsdf,aes(x = x,y = y,label = well),size = cent,colour = "steelblue",alpha = 0.4,fontface = "bold")
        if(!centerwell) p <- p + geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)
        p <- p + scale_x_continuous(expand = c(0,0))+
        scale_y_reverse(expand = c(0,0))+
        theme_bw(base_size = 5)+
        labs(title = paste0(wells, " Well Plate: ",fname,"  |  Wells & Spots"))+
        coord_fixed()+
        theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 5),legend.position = "none",
              axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
              axis.line = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
              panel.border = element_blank())
      ggsave(filename = paste0(fname,"-1-Wells.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
      cat(paste0(fname,"-1-Wells.png exported.\n"))
      rm(p,prows,pcols)
    }

    #---------------------------------------------------------------------------
    #Coverage

    if(coverage)
    {
      cat(paste0("Computing coverage...\n"))
      #covlist <- split(dframe,dframe$well)
      #covlen <- length(covlist)
      covlisthulls <- vector("list",length=wells)
      covlistarea <- vector("list",length=wells)
      pb <- txtProgressBar(min = 0, max = wells, style = 3)
      coverageloop = 1
      for(coverageloop in 1:wells)
      {
        setTxtProgressBar(pb, coverageloop)
        covwell <- subset(dframe,dframe$well == coverageloop)[,c("x","y","well")]
        #covwell <- covlist[[coverageloop]][,c("x","y","well")]

        if(nrow(covwell) > 3)
        {
          covlisthulls[[coverageloop]] <- covwell[chull(x = covwell$x,y = covwell$y),]
          areatemp <- round(abs(polyarea(covlisthulls[[coverageloop]]$x,covlisthulls[[coverageloop]]$y)),0)
          covlistarea[[coverageloop]] <- data.frame(well=coverageloop,coverage_pxsq=areatemp,stringsAsFactors = FALSE)
        }else{
          covlisthulls[[coverageloop]] <- data.frame(x=NA,y=NA,well=coverageloop,stringsAsFactors=FALSE)
          covlistarea[[coverageloop]] <- data.frame(well=coverageloop,coverage_pxsq=NA,stringsAsFactors = FALSE)
        }

      }
      close(pb)
      covdfhulls <- plyr::rbind.fill(covlisthulls)
      covdfarea <- plyr::rbind.fill(covlistarea)
      rm(covlisthulls,covlistarea,coverageloop,areatemp,pb)
    }

    #---------------------------------------------------------------------------
    #Minimum spanning network

    if(msd)
    {
      cat(paste0("Computing minimum spanning tree and distance...\n"))
      #msdlist <- split(dframe,dframe$well)
      #msdlen <- length(msdlist)
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
        if(nrow(msdwell)>2) msdwell <- msdwell[-which(duplicated(round(msdwell,1))),]

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
          msdvec[msdloop] <- 0
          msdlinelist[[msdloop]] <- data.frame(x0=NA,y0=NA,x1=NA,y1=NA,well=msdloop,msd_px=NA,stringsAsFactors=FALSE)
        }

      }
      close(pb)
      msddf <- plyr::rbind.fill(msdlinelist)
      rm(msdlinelist,msdwell,msdloop,pb)

      if(exportplot)
      {
        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          #geom_rect(data = wellsdf,aes(xmin = h1,xmax = h2,ymin = v1,ymax = v2),colour = "grey90",fill = "white",size = 0.3)+
          #geom_point(data = dframe,aes(x,y,colour=factor(well)),size = 0.3)+
          #geom_path(data = msddf, aes(x,y,group = well,colour=factor(well)), size=0.2)+
          geom_segment(data=msddf,aes(x=x0,xend=x1,y=y0,yend=y1,group=well,colour=factor(well)),size=0.2,na.rm = TRUE)
        if(centerwell) p <- p + geom_text(data = wellsdf,aes(x = x,y = y,label = well),size = cent,colour = "steelblue",alpha = 0.4,fontface = "bold")
        if(!centerwell) p <- p + geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)
        p <- p + scale_x_continuous(expand = c(0,0))+
          scale_y_reverse(expand = c(0,0))+
          theme_bw(base_size = 5)+
          labs(title = paste0(wells, " Well Plate: ",fname,"  |  Minimum Spanning Distance"))+
          coord_fixed()+
          theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 5),legend.position = "none",
                axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
                axis.line = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
                panel.border = element_blank())
        ggsave(filename = paste0(fname,"-2-Msd.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fname,"-2-Msd.png exported.\n"))
      }
      rm(msddf)
    }

    #---------------------------------------------------------------------------
    # alphahull

    if(alphahull)
    {
      cat(paste0("Computing alphahull...\n"))
      #ahlist <- split(dframe,dframe$well)
      #ahlen <- length(ahlist)
      ahvec <- vector(length=wells)
      ahverlist <- vector("list",length=wells)
      pb <- txtProgressBar(min = 0, max = wells, style = 3)
      for(ahloop in 1:wells)
      {
        setTxtProgressBar(pb, ahloop)
        #ahwell <- ahlist[[ahloop]][,c("x","y")]
        ahwell <- subset(dframe,dframe$well == ahloop)[,c("x","y")]

        if(nrow(ahwell) > 3)
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

          ah1 <- alphahull::ahull(x=ahwell$x,ahwell$y,alpha=alphavalue)
          ahvec[ahloop] <- alphahull::areaahull(ah1)
          ahverlist[[ahloop]] <- as.data.frame(ah1$ashape.obj$edges)[,c("x1","y1","x2","y2")]
          ahverlist[[ahloop]]$well <- rep(ahloop,nrow(ahverlist[[ahloop]]))
          #segments(x0 = t1$x1,y0 = t1$y1,x1 = t1$x2,y1=t1$y2,col="red")
        }else{
          ahvec[ahloop] <- NA
          ahverlist[[ahloop]] <- data.frame(x1=NA, y1=NA, x2=NA, y2=NA, well=ahloop, stringsAsFactors = FALSE)
        }
      }
      close(pb)
      ahdf <- plyr::rbind.fill(ahverlist)
      rm(ahverlist,ahwell,ahloop,pb)

      if(exportplot)
      {
        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          geom_point(data = dframe,aes(x,y),colour="lightgrey",size = 0.25,alpha=0.7,na.rm = TRUE)+
          geom_segment(data=ahdf,aes(x=x1,xend=x2,y=y1,yend=y2,colour=factor(well)),size=0.2,na.rm = TRUE)
        if(centerwell) p <- p + geom_text(data = wellsdf,aes(x = x,y = y,label = well),size = cent,colour = "steelblue",alpha = 0.4,fontface = "bold")
        if(!centerwell) p <- p + geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)
        p <- p + scale_x_continuous(expand = c(0,0))+
          scale_y_reverse(expand = c(0,0))+
          theme_bw(base_size = 5)+
          labs(title = paste0(wells, " Well Plate: ",fname,"  |  Alpha Hulls at alpha ",alphavalue))+
          coord_fixed()+
          theme(plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 5),legend.position = "none",
                axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
                axis.line = element_blank(),panel.background = element_blank(),panel.grid = element_blank(),
                panel.border = element_blank())
        ggsave(filename = paste0(fname,"-3-AlphaHulls.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fname,"-3-AlphaHulls.png exported.\n"))
      }
      rm(ahdf)
    }

    #---------------------------------------------------------------------------
    #Frame linking

    cat(paste0("Linking frames...\n"))
    framelist <- split(dframe,dframe$slice)
    rm(dframe)

    #loop over frames
    #removepointpos <- data.frame()
    dflist <- vector("list",length = length(framelist))
    numberofpoints <- length(framelist)
    pb <- txtProgressBar(min = 0, max = numberofpoints, style = 3)
    frameloop=1
    for(frameloop in 1:numberofpoints)
    {
      setTxtProgressBar(pb, frameloop)
      #if(frameloop==51) stop("Error.")
      if(frameloop == numberofpoints) break
      #if((frameloop %% 100) == 0) cat(paste0("Frame: ",frameloop,"\n"))
      #get exact points from frame 1
      if (frameloop == 1)
      {
        currentframe <- framelist[[frameloop]]
        currentframe$id <- 1:nrow(currentframe)
        #ids <- currentframe$id
      }
      nextframe <- framelist[[frameloop+1]]

      #loop over points
      #positionlist <- vector("list",length = nrow(currentframe))
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
        #         plot(x = currentframe$x,y = currentframe$y,pch = "+",xlab = "",ylab = "",cex = 1.5)
        #         points(x = nextframe$x,y = nextframe$y,pch = 16,col = "#3182bdB3",cex = 1.5)
        #         title(main = paste0("Frame: ",frameloop),line = 0.5,adj = 0)
        #         legend(x = 470,y=480,legend = c("Current frame","Next frame"),col = c("black","#3182bdB3"),pch = c(3,16),ncol = 2,bty = "n",xpd = T)
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

      currentframe <- tempdf
      #add frame number
      tempdf$frame <- rep(frameloop,nrow(currentframe))
      #save df to list
      dflist[[frameloop]] <- tempdf

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
    rm(dflist,xylist,tempdf,pb,currentframe,nextframe,foundpointpos,foundpointtemp,
       searchpointpos,pointloop,numberofpoints,framelist,frameloop)

    if(exportdata)
    {
      write.table(x = df1, file = paste0(fname,"-Tracks.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fname,"-Tracks.txt exported.\n"))
    }

    #---------------------------------------------------------------------------

    #Distance and speed calculations
    cat(paste0("Computing distance, speed and activity...\n"))
    ##ids <- as.numeric(levels(factor(df1$id)))
    #idlen <- length(ids)
    pb <- txtProgressBar(min = 0, max = wells, style = 3)
    welllist <- vector("list",length = wells)
    secspeedppslist <- vector("list",length=wells)
    #loop over id
    loop1=1
    for(loop1 in 1:wells)
    {
      setTxtProgressBar(pb, loop1)
      currentwell <- subset(df1, df1$well == loop1)

      if(nrow(currentwell) > 1)
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

        #speed pps by well
        #well1 <- as.numeric(as.character(currentwell$well))[loop1]
        #secspeedppslist[[loop1]] <- data.frame(well=rep(well1,length(secspeedpps)),id=rep(loop1,
        #                            length(secspeedpps)),speedpps=secspeedpps,speedmmps=secspeedmps,stringsAsFactors = FALSE)

        # create track stats df
        welllist[[loop1]] <- data.frame(plate = fname,id = as.numeric(currentwell$id[1]),well = loop1,
                                       dist_px = round(totdist,0),dist_mm = round(totdist/mm,2),speed_mean_pps = round(mean(secspeedpps),2),
                                       speed_mean_mmps = round(mean(secspeedmps),2),speed_max_pps = round(max(secspeedpps),2),
                                       speed_max_mmps = round(max(secspeedmps),2),duration_fr = nrow(currentwell),
                                       duration_sec = round(nrow(currentwell)/fps,0), activity = round(length(which(secspeedmps > activitydist))/length(secspeedmps),3),
                                       fps = fps,mm = mm,stringsAsFactors = FALSE)
      }else{
        # create track stats df
        welllist[[loop1]] <- data.frame(plate = fname,id = as.numeric(currentwell$id[1]),well = loop1,
                                       dist_px = NA,dist_mm = NA,speed_mean_pps = NA,
                                       speed_mean_mmps = NA,speed_max_pps = NA,
                                       speed_max_mmps = NA,duration_fr = NA,
                                       duration_sec = NA, activity = NA,
                                       fps = fps,mm = mm,stringsAsFactors = FALSE)
      }


    }
    close(pb)
    distspeedtemp <- plyr::arrange(plyr::rbind.fill(welllist),well)
    rm(totdist,currentwell,framelen,distancevector,loop2,d,loop1,pb,secspeed,secspeedpps,secspeedmps,secspeedppslist,welllist)

    #---------------------------------------------------------------------------
    #combining data
    #merge coverage
    if(coverage) distspeedtemp <- base::merge(x = distspeedtemp,y = covdfarea, by = "well")
    #merge msd
    if(msd)
    {
      distspeedtemp$msd_px <- round(msdvec,0)
      distspeedtemp$msd_mm <- round(msdvec/mm,2)
      rm(msdvec)
    }
    #merge alphahull
    if(alphahull)
      {
      distspeedtemp$alphahull_pxsq <- round(ahvec,2)
        rm(ahvec)
      }

    if(exportdata)
    {
      write.table(x = distspeedtemp, file = paste0(fname,"-Features.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fname,"-Features.txt exported.\n"))
    }

    #---------------------------------------------------------------------------
    #plots

    #create track statistics plot
    if(exportplot)
    {
      cat(paste0("Generating plots.\n"))

      txtdf <- distspeedtemp[,c("id","well","dist_px")]
      #move text for labelling on plot
      txtdf$dist_px[is.na(txtdf$dist_px)] <- 0
      txtdf$dist_px[txtdf$dist_px > 0] <- round((0.90*txtdf$dist_px[txtdf$dist_px > 0]),0)
      txtdf$dist_px[txtdf$dist_px == 0] <- NA

      if(wells == 24) txtszbp <- 1.8
      if(wells == 48) txtszbp <- 1.5
      if(wells == 96) txtszbp <- 1.2

      # barplot distances
      p <- ggplot()+
        geom_bar(data = distspeedtemp,aes(x = well,y = dist_px),stat = "identity",fill = "steelblue",na.rm = TRUE)+
        scale_x_continuous(expand = c(0,0),breaks = distspeedlist[[fileloop]]$well)+
        scale_y_continuous(expand = c(0,0))+
        geom_text(data = txtdf,aes(x = well,y = dist_px,label = dist_px),colour = "grey80",size = txtszbp)+
        theme_bw(base_size = 5)+
        ggtitle(paste0(wells, " Well Plate: ",fname,"  |  Distance"))+
        ylab("Distance (Pixels)")+
        xlab("Well")+
        coord_flip()+
        theme(plot.title = element_text(hjust = 0,lineheight = 1.2,colour = "grey40",size = 5),
              axis.title = element_text(colour = "grey40"),axis.text = element_text(colour = "grey40"),
              axis.ticks = element_line(colour = "grey40",size = 0.3),panel.border = element_blank())
      ggsave(filename = paste0(fname,"-6-Distances.png"),plot = p,height = 12,width = 9,units = "cm",dpi = 300,type = "cairo")
      cat(paste0(fname,"-6-Distances.png exported.\n"))
      rm(p)

      wellsdf1 <- base::merge(txtdf,wellsdf,by = "well")
      if(wells == 24) wellsdf1$laby1 <- wellsdf1$y + 47
      if(wells == 48) wellsdf1$laby1 <- wellsdf1$y + 30
      if(wells == 96) wellsdf1$laby1 <- wellsdf1$y + 20
      #wellsdf1$laby1 <- wellsdf1$y+(0.38*(diff(verticalrange)/prows))
      df1.1 <- df1
      df1.1$order <- 1:nrow(df1.1)
      df2 <- base::merge(df1.1,txtdf,by="id",all.x=T,all.y=F)
      rm(df1.1)
      df2 <- plyr::arrange(df2,order)
      df2$linktype <- factor(df2$linktype,levels=c("Single","Duplicated","Nearest"))

      # tracks plot
      p <- ggplot()+
        geom_blank(data = edgedf,aes(x,y))+
        #geom_rect(data = wellsdf,aes(xmin = h1,xmax = h2,ymin = v1,ymax = v2),colour = "grey90",fill = "white",size = 0.3)+
        #geom_point(data = df1,aes(x = x,y = y,colour = factor(id),group = factor(id)),size = 0.2)+
        geom_path(data = df2,aes(x = x,y = y,colour = linktype,group = factor(id)),size = 0.25,alpha=0.8,na.rm = TRUE)

        if(centerwell) p <- p + geom_text(data = wellsdf,aes(x = x,y = y,label = well),size = cent,colour = "steelblue",alpha = 0.4,fontface = "bold")
        if(!centerwell) p <- p + geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "steelblue",alpha = 0.4,fontface = "bold",hjust=0.5,vjust=0.5)

        p <- p + geom_text(data = wellsdf1,aes(x = x,y = laby1,label = dist_px),size = txtsz,colour = "grey20",alpha = 0.6,fontface = "bold",na.rm = TRUE)+
        theme_bw()+
        scale_x_continuous(expand = c(0,0)) +
        scale_y_reverse(expand = c(0,0))+
        scale_colour_manual(values=c("#31a354","#3182bd","#de2d26"))+
        labs(title = paste0(wells, " Well Plate: ",fname,"  |  Wells & Tracks"))+
        coord_fixed()+
        theme_bw(base_size = 5)+
        theme(legend.position = "none",axis.ticks = element_blank(),legend.text = element_blank(),
              axis.text = element_blank(),axis.title = element_blank(),legend.title = element_blank(),
              plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 5),panel.border = element_blank(),
              panel.background = element_blank(),panel.grid = element_blank())
      ggsave(filename = paste0(fname,"-4-Tracks.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
      cat(paste0(fname,"-4-Tracks.png exported.\n"))

      if(coverage)
      {
        covdata <- base::merge(covdfarea,covdfhulls,by="well")
        covwells <- base::merge(wellsdf1,covdfarea,by = "well")
        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          geom_rect(data = wellsdf,aes(xmin = h1,xmax = h2,ymin = v1,ymax = v2),colour = "grey90",fill = "white",size = 0.3,na.rm = TRUE)+
          geom_point(data = df1,aes(x,y),colour="lightgrey",size = 0.25,alpha=0.6,na.rm = TRUE)+
          geom_polygon(data = covdata,aes(x = x,y = y,colour = coverage_pxsq,group = factor(well)),fill=NA,size = 0.3,na.rm = TRUE)+
          geom_point(data = covdata,aes(x = x,y = y,colour = coverage_pxsq,fill = coverage_pxsq,group = factor(well)),size = 0.8,alpha=0.9,na.rm = TRUE)
        #well number large center
        if(centerwell) p <- p + geom_text(data = wellsdf,aes(x = x,y = y,label = well),size = cent,colour = "grey20",alpha = 0.4,fontface = "bold",na.rm = TRUE)
        #well number small corner
        #geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "grey20",alpha = 0.3,fontface = "bold")+
        #distance in pixels
        p <- p +
          geom_text(data = covwells,aes(x = x,y = laby1,label = coverage_pxsq),size = txtsz,colour = "grey20",alpha = 0.6,fontface = "bold",na.rm = TRUE)+
          theme_bw()+
          scale_x_continuous(expand = c(0,0)) +
          scale_y_reverse(expand = c(0,0))+
          labs(title = paste0(wells, " Well Plate: ",fname,"  |  Coverage"))+
          coord_fixed()+
          theme_bw(base_size = 5)+
          theme(legend.position = "none",axis.ticks = element_blank(),legend.text = element_blank(),
                axis.text = element_blank(),axis.title = element_blank(),legend.title = element_blank(),
                plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 5),panel.border = element_blank(),
                panel.background = element_blank(),panel.grid = element_blank())
        ggsave(filename = paste0(fname,"-5-Coverage.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fname,"-5-Coverage.png exported.\n"))
        rm(covdata,covwells)
      }

      distspeedlist[[fileloop]] <- distspeedtemp
      rm(p,wellsdf1,txtdf,df2)
    }

    cat(paste0("\n---------------------------------------------------------\n\n"))
  }

  distspeeddata <- plyr::rbind.fill(distspeedlist)
  rm(distspeedlist)

  if(exportdata)
  {
    if(flen > 1)
    {
      write.table(x = distspeeddata, file = paste0("Combined-Features.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0("Combined-Features.txt exported.\n"))
    }
  }
  cat(paste0("Completed in ",format((Sys.time() - currtime),format="%S",digits=3),".\n"))
  return(list(RawData = df1,Features = distspeeddata))
}
