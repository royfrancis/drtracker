# LARVAE TRACKING AND DISTANCE CALCULATION
# v1.0.0
# 12-Oct-2015

# Assumptions
# 24 or 48 well plate
# plate is constant throughout the film. absolutely no movement.
# edges and larvae positions are marked in frame 1
# single larvae per well
# if larval position is missing in a frame, old position is duplicated
# if two positions are available, nearest is chosen.

#Issues
# fix label positions for different well plates

#load libraries
#library(plyr) #rbind.fill
#library(RANN) #nn2
#library(ggplot2)

# FUNCTION LTRACK
#' Track larval movement in 24 or 48 well assay plates.
#' @description Tracks single larvae in 24 or 48 well assay plates from xy spot data
#' and computes distance and speed. Exports data as text files and generates plots.
#' @param files A character or vector of paths or filenames. See details.
#' @param wells A numeric indicating plate format. Either 24 or 48.
#' @param markededges A logical indicating if corners of plates have been marked. See details.
#' @param fps A numeric indicating framerate of video in number of frames per second.
#' @param mm A numeric indicating number of pixels in 1 mm.
#' @param coverage A logical indicating if the coverge must be computed. See details.
#' @param filenamediscard A character for part of the filename to be removed.
#' @param exportdata A logical indicating if data tables must be exported to working directory.
#' @param exportplot A logical if results must be plotted and exported to working directory.
#' @param follow A character indicating if the algorithm must be followed at every step.
#' Options are 'interactive', 'track' or 'none'. See details.
#' @param centerwell A logical indicating if the well number must be printed at center of the well on plots.
#' @return Returns a list containing two components; 'RawData' and 'Track statistics'.\cr
#' The 'RawData' is a dataframe for one input file or a list for 2 or more input files.
#' The 'RawData' contains xy coordinates, slice, well, id and comment per spot. The
#' id links the spots as tracks. The comment shows the method used to link a spot.
#' 'Single' means that there was only one spot in the next frame. This may be the correct
#' spot or a stray spot. 'Nearest' means there was more than one spot in the next frame
#' and the nearest spot was chosen. 'Duplicate' means that no spot was present in the next frame
#' and the previous spot was duplicated in next frame.\cr
#' 'Track Statistics' contains plate name, id, well, total distance in pixels and mm, mean speed in pixels per frame(ppf) and
#' in mm per second (mps), max speed in pixels per sec (pps) and mm per sec (mps), duration of
#' the whole sequence in frames (fr) and seconds (sec), framerate (fps) and calibration,
#' number of pixels in one mm (mm). If \code{coverage=TRUE}, then coverage_pxsq and coverage_mmsq are added.\cr
#' If \code{exportdata = T}, then three text files are exported: Wells, Tracks and TrackStatistics.
#' If more than one file was selected, a Combined-TrackStatistics file is also exported.
#' If \code{exportplot = T}, then three-five figures are exported: EdgeSpots (if edge spots are used),
#' Wells & Spots, Wells & Tracks, Coverage and Track Statistics.
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
#' \strong{filenamediscard}\cr
#' The file name of the input file is used on plots and text files for identification.
#' Part of the filename to be removed can be indicated here. '.txt' is removed by default.\cr
#' \strong{follow}\cr
#' In 'interactive' mode, a plot is shown at every frame and waits for user input.
#' In 'track' mode, the track path creation is shown in real-time.
#' @import plyr
#' @import RANN
#' @export
#'
ltrack <- function(files = NULL, wells=24, markededges = TRUE, fps = 25, mm = 5.4, coverage = TRUE,
                   filenamediscard = ".txt", exportdata = TRUE, exportplot = TRUE,
                   follow = "none", centerwell = FALSE)
{
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
  for(fileloop in 1:flen)
  {
    cat(paste0("Processing file ",fileloop," of ",flen," files.\n"))
    cat(paste0("Current file: ",files[fileloop],".\n"))

    fname <- sub(filenamediscard,"",basename(files[fileloop]))
    dframe <- read.delim(files[fileloop],header = T,stringsAsFactors = F)
    colnames(dframe) <- tolower(colnames(dframe))
    if(all(c("id","comment","frame") %in% colnames(dframe))) warning("Input file may be incorrect.\n")
    dframe <- dframe[,c("x","y","slice")]
    dframe$x <- as.integer(round(dframe$x,0))
    dframe$y <- as.integer(round(dframe$y,0))

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
          geom_point(shape = 1,size = 3)+
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
        geom_rect(data = wellsdf,aes(xmin = h1,xmax = h2,ymin = v1,ymax = v2),colour = "grey90",fill = "white",size = 0.3)+
        geom_point(data = dframe,aes(x,y,colour = well),size = 0.3)
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
      covlist <- split(dframe,dframe$well)
      covlen <- length(covlist)
      covlisthulls <- vector("list",length=covlen)
      covlistarea <- vector("list",length=covlen)
      pb <- txtProgressBar(min = 0, max = covlen, style = 3)
      for(coverageloop in 1:covlen)
      {
        setTxtProgressBar(pb, coverageloop)
        covwell <- covlist[[coverageloop]]
        covlisthulls[[coverageloop]] <- covwell[chull(x = covwell$x,y = covwell$y),c("x","y","well")]
        covlistarea[[coverageloop]] <- data.frame(well=unique(covlisthulls[[coverageloop]]$well),area=round(abs(polyarea(covlisthulls[[coverageloop]]$x,covlisthulls[[coverageloop]]$y)),0),stringsAsFactors = FALSE)
      }
      close(pb)
      covdfhulls <- plyr::rbind.fill(covlisthulls)
      covdfarea <- plyr::rbind.fill(covlistarea)
      rm(covlist,covlen,covlisthulls,covlistarea,coverageloop,pb)
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
      #xylist1 <- vector("list",length = nrow(currentframe))
      #lenvec <- vector("numeric",length = nrow(currentframe))
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
            foundpointpos$comment <- "Single"
          }

          #if more than one point was found
          if(nrow(foundpointtemp) > 1)
          {
#             if(follow == "interactive")
#             {
#               cols <- alphacol(rainbow(nrow(foundpointtemp)))
#               points(foundpointtemp$x,foundpointtemp$y,pch = 15,cex = 1.2,col = cols)
#             }

            # check if blacklist > 0
            #         if(nrow(removepointpos) > 0)
            #         {
            #           for(k in 1:nrow(removepointpos))
            #           {
            #             if(any(duplicated(rbind(removepointpos[,c("x","y")],foundpointtemp[,c("x","y")]))))
            #             {
            #               foundpointtemp <- foundpointtemp[-k,]
            #             }
            #           }
            #         }

            if(nrow(foundpointtemp) > 1)
            {
              #logicalcheck <- (((foundpointtemp$x - searchpointpos$x)^2 + (foundpointtemp$y - searchpointpos$y)^2) <= radius^2)
              #           if(sum(logicalcheck) == 1)
              #           {
              #             foundpointpos <- foundpointtemp[which(logicalcheck),]
              #             foundpointpos$id <- searchpointpos$id
              #             foundpointpos$comment <- "Nearest"
              #           }

              foundpointpos <- foundpointtemp[as.numeric(RANN::nn2(foundpointtemp[,c("x","y")],searchpointpos[,c("x","y")],k = 1)$nn.idx),]
              foundpointpos$id <- searchpointpos$id
              foundpointpos$comment <- "Nearest"

              #           else
              #             {
              #             cols <- alphacol(rainbow(nrow(foundpointtemp)))
              #             if(interactive)
              #             {
              #               points(foundpointtemp$x,foundpointtemp$y,pch=15,cex=1.2,col=cols)
              #               text(searchpointpos$x,searchpointpos$y,labels=1:nrow(foundpointtemp),offset=0.6,pos=4,cex=1,col=cols)
              #             }else
              #             {
              #               plot(x=currentframe$x,y=currentframe$y,pch="+",xlab="",ylab="",cex=1.5)
              #               title(main=paste0("Frame: ",frameloop),line=0.5,adj=0)
              #               points(foundpointtemp$x,foundpointtemp$y,pch=15,cex=1.2,col=cols)
              #               text(foundpointtemp$x,foundpointtemp$y,labels=1:nrow(foundpointtemp),offset=0.6,pos=4,cex=1,col=cols)
              #             }
              #
              #             cat("Blacklist:\n")
              #             print(removepointpos)
              #             cat("Current Conflicts:\n")
              #             print(foundpointtemp)
              #
              #             #cat("Remove conflicts for this loop or all loops?\n1.This Loop 2.All Loops Enter a number:\n")
              #             #selecttype <- as.numeric(scan(n=1))
              #             selecttype <- 1
              #             if(selecttype == 1)
              #             {
              #               cat("Select point to keep. Enter a number:\n")
              #               keeppoint <- as.numeric(scan(n=1))
              #               foundpointpos <- foundpointtemp[keeppoint,]
              #               foundpointpos$id <- searchpointpos$id
              #               foundpointpos$comment <- "Selected"
              #             }
              #             if(selecttype == 2)
              #             {
              #               cat("Select point to keep. Enter a number:\n")
              #               keeppoint <- as.numeric(scan(n=1))
              #               foundpointpos <- foundpointtemp[keeppoint,]
              #               foundpointpos$id <- searchpointpos$id
              #               foundpointpos$comment <- "Selected"
              #               removepointpos <- rbind(removepointpos,foundpointtemp[-keeppoint,])
              #             }
            }
          }

        }else{
          # if no match was found, duplicate point in next frame
          foundpointpos <- searchpointpos
          foundpointpos$slice <- frameloop
          foundpointpos$comment <- "Duplicated"
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
        #pframe$comment <- factor(c(rep("None",nrow(currentframe)),as.character(tempdf$comment)),levels=c("None","Single","Nearest","Duplicated"))
        p <- ggplot()+
          geom_point(data=pframe,aes(x=x,y=y,shape=factor(frame),col=factor(frame)),size=6,alpha=0.6)+
          geom_line(data=pframe,aes(x=x,y=y,group=factor(id)),linetype=2)+
          geom_text(data=currentframe,aes(x=x+20,y=y,label=id),size=3)+
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
    cat(paste0("Computing distance and speed...\n"))
    ids <- as.numeric(levels(factor(df1$id)))
    idlen <- length(ids)
    pb <- txtProgressBar(min = 0, max = idlen, style = 3)
    idlist <- vector("list",length = idlen)
    #loop over id
    for(idloop in 1:idlen)
    {
      setTxtProgressBar(pb, idloop)
      currentid <- subset(df1, id == ids[idloop])
      framelen <- nrow(currentid)
      distancevector <- vector(length = framelen)
      idframeloop = 1
      secspeed <- vector(length = fps)
      secspeedpps <- vector()
      secspeedmps <- vector()
      k = 1
      #loop over frames
      while(idframeloop < framelen)
      {
        if(k == (fps+1))
        {
          # save speed pps
          secspeedpps <- c(secspeedpps,secspeed)
          # save speed mps
          secspeedmps <- c(secspeedmps,(secspeed/mm))
          secspeed <- vector(length = fps)
          k <- 1
        }
        #calculate euclidian distance
        d <- sqrt(diff(c(currentid$x[idframeloop],currentid$x[idframeloop+1]))^2 + diff(c(currentid$y[idframeloop],currentid$y[idframeloop+1]))^2)
        secspeed[k] <- d
        distancevector[idframeloop] <- d
        idframeloop <- idframeloop+1
        k = k+1
      }
      # sum total distance
      totdist <- sum(distancevector)
      # create track stats df
      idlist[[idloop]] <- data.frame(plate = fname,id = as.numeric(ids[idloop]),well = as.numeric(as.character(currentid$well))[idloop],
                                     dist_px = round(totdist,0),dist_mm = round(totdist/mm,2),speed_mean_ppf = round(totdist/framelen,2),
                                     speed_mean_mps = round(((totdist/mm)/(framelen/fps)),2),speed_max_pps = round(max(secspeedpps),2),
                                     speed_max_mps = round(max(secspeedmps),2),duration_fr = nrow(currentid),
                                     duration_sec = round(nrow(currentid)/fps,0),fps = fps,mm = mm,stringsAsFactors = FALSE)
    }
    close(pb)
    distspeedtemp <- plyr::arrange(plyr::rbind.fill(idlist),well)
    if(!coverage) distspeedlist[[fileloop]] <- distspeedtemp
    if(coverage) distspeedlist[[fileloop]] <- base::merge(x = distspeedtemp,y = covdfarea, by = "well")

    if(exportdata)
    {
      write.table(x = distspeedlist[[fileloop]], file = paste0(fname,"-TrackStatistics.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0(fname,"-TrackStatistics.txt exported.\n"))
    }

    rm(ids,idlen,idlist,totdist,currentid,framelen,distancevector,idframeloop,d,idloop,pb,secspeed,secspeedpps,secspeedmps)

    #---------------------------------------------------------------------------
    #plots

    #create track statistics plot
    if(exportplot)
    {
      cat(paste0("Generating plots.\n"))

      txtdf <- distspeedlist[[fileloop]][,c("id","well","dist_px")]
      txtdf$dist_px[txtdf$dist_px > 0] <- round((0.90*txtdf$dist_px[txtdf$dist_px > 0]),0)
      txtdf$dist_px[txtdf$dist_px == 0] <- NA

      if(wells == 24) txtszbp <- 1.8
      if(wells == 48) txtszbp <- 1.5
      if(wells == 96) txtszbp <- 1.2

      # barplot distances
      p <- ggplot()+
        geom_bar(data = distspeedlist[[fileloop]],aes(x = well,y = dist_px),stat = "identity",fill = "steelblue")+
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
      ggsave(filename = paste0(fname,"-4-Distances.png"),plot = p,height = 12,width = 9,units = "cm",dpi = 300,type = "cairo")
      cat(paste0(fname,"-4-Distances.png exported.\n"))
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

      # tracks plot
      p <- ggplot()+
        geom_blank(data = edgedf,aes(x,y))+
        #geom_rect(data = wellsdf,aes(xmin = h1,xmax = h2,ymin = v1,ymax = v2),colour = "grey90",fill = "white",size = 0.3)+
        #geom_point(data = df1,aes(x = x,y = y,colour = factor(id),group = factor(id)),size = 0.2)+
        geom_path(data = df2,aes(x = x,y = y,colour = dist_px,group = factor(id)),size = 0.25)
        #well number large center
        if(centerwell) p <- p + geom_text(data = wellsdf,aes(x = x,y = y,label = well),size = cent,colour = "grey20",alpha = 0.4,fontface = "bold")
        #well number small corner
        #geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "grey20",alpha = 0.3,fontface = "bold")+
        #distance in pixels
        p <- p + geom_text(data = wellsdf1,aes(x = x,y = laby1,label = dist_px),size = txtsz,colour = "grey20",alpha = 0.6,fontface = "bold")+
        theme_bw()+
        scale_x_continuous(expand = c(0,0)) +
        scale_y_reverse(expand = c(0,0))+
        labs(title = paste0(wells, " Well Plate: ",fname,"  |  Wells & Tracks"))+
        coord_fixed()+
        theme_bw(base_size = 5)+
        theme(legend.position = "none",axis.ticks = element_blank(),legend.text = element_blank(),
              axis.text = element_blank(),axis.title = element_blank(),legend.title = element_blank(),
              plot.title = element_text(lineheight = 1.2,hjust = 0,colour = "grey40",size = 5),panel.border = element_blank(),
              panel.background = element_blank(),panel.grid = element_blank())
      ggsave(filename = paste0(fname,"-2-Tracks.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
      cat(paste0(fname,"-2-Tracks.png exported.\n"))

      if(coverage)
      {
        covdata <- base::merge(covdfarea,covdfhulls,by="well")
        covwells <- base::merge(wellsdf1,covdfarea,by = "well")
        p <- ggplot()+
          geom_blank(data = edgedf,aes(x,y))+
          geom_rect(data = wellsdf,aes(xmin = h1,xmax = h2,ymin = v1,ymax = v2),colour = "grey90",fill = "white",size = 0.3)+
          geom_polygon(data = covdata,aes(x = x,y = y,colour = area,fill = area,group = factor(well)),size = 0.3,alpha=0.7)+
          geom_point(data = covdata,aes(x = x,y = y,colour = area,fill = area,group = factor(well)),size = 0.8,alpha=0.9)
        #well number large center
        if(centerwell) p <- p + geom_text(data = wellsdf,aes(x = x,y = y,label = well),size = cent,colour = "grey20",alpha = 0.4,fontface = "bold")
        #well number small corner
        #geom_text(data = wellsdf,aes(x = labx,y = laby,label = well),size = txtsz,colour = "grey20",alpha = 0.3,fontface = "bold")+
        #distance in pixels
        p <- p +
          geom_text(data = covwells,aes(x = x,y = laby1,label = area),size = txtsz,colour = "grey20",alpha = 0.6,fontface = "bold")+
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
        ggsave(filename = paste0(fname,"-3-Coverage.png"),plot = p,height = 8,width = 12,units = "cm",dpi = 300,type = "cairo")
        cat(paste0(fname,"-3-Coverage.png exported.\n"))
        rm(covdata,covwells)
      }

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
      write.table(x = distspeeddata, file = paste0("Combined-TrackStatistics.txt"),quote = FALSE,row.names = FALSE,sep = "\t",dec = ".")
      cat(paste0("Combined-TrackStatistics.txt exported.\n"))
    }
  }
  cat(paste0("Done.\n"))
  return(list(RawData = df1,TrackStatistics = distspeeddata))
}
