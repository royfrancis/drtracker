library(testthat)
library(drtracker)

#devtools::test()

#Preparation
deleteoutput = TRUE
#create a new folder and set as wd
currwd <- getwd()
dir.create(paste(currwd,"/drtrackerDemo",sep=""))
setwd(paste(currwd,"/drtrackerDemo",sep=""))
#read sample file from drtracker package
ll <- system.file("files/testfile.txt",package="drtracker")

#-------------------------------------------------------------------------------

#ASSIGNWELLS
context("assignWells")
aw <- assignWells(files = ll,wells = 24,markededges = T,exportdata = T,exportplot = T)

test_that("Is output dataframe?",{
  expect_equal(class(aw),"data.frame")
})

test_that("dataframe columns",{
  expect_equal(all(c("row","plate","x","y","slice","well") %in% colnames(aw)),TRUE)
})
rm(aw)

test_that("exportdata",{
  expect_equal("testfile-AssignedWells.txt" %in% list.files(),TRUE)
  expect_equal("testfile-EdgeMarks.txt" %in% list.files(),TRUE)
})

test_that("exportplot",{
  expect_equal("testfile-EdgeMarks.png" %in% list.files(),TRUE)
  expect_equal("testfile-AssignedWells.png" %in% list.files(),TRUE)
})

test_that("Error: no input",{
  expect_error(assignWells())
})

#-------------------------------------------------------------------------------

#LINKFRAMES
context("linkFrames")
lf <- linkFrames(files = "testfile-AssignedWells.txt",wells = 24,exportdata = T, exportplot = T)

test_that("Is output dataframe?",{
  expect_equal(class(lf),"data.frame")
})

test_that("dataframe columns",{
  expect_equal(all(c("row","plate","x","y","slice","well","id","linktype","frame") %in% colnames(lf)),TRUE)
})
rm(lf)

test_that("exportdata",{
  expect_equal("testfile-Tracks.txt" %in% list.files(),TRUE)
})

test_that("exportplot",{
  expect_equal("testfile-Tracks.png" %in% list.files(),TRUE)
})

test_that("Error: no input",{
  expect_error(linkFrames())
})

#-------------------------------------------------------------------------------

#TRACKFEATURES
context("trackFeatures")
tf <- trackFeatures(files = "testfile-Tracks.txt",wells = 24, fps = 25, mm = 5.4,
                    activitydist = 5, exportdata = T, exportplot = T)

test_that("Is output dataframe?",{
  expect_equal(class(tf),"data.frame")
})

test_that("dataframe columns",{
  expect_equal(all(c("plate","id","well","dist_px","dist_mm","speed_mean_pps","speed_mean_mmps","speed_max_pps",
                 "speed_max_mmps","duration_fr","duration_sec","activity","fps","mm") %in% colnames(tf)),TRUE)
})
rm(tf)

test_that("exportdata",{
  expect_equal("testfile-TrackFeatures.txt" %in% list.files(),TRUE)
})

test_that("Error: no input",{
  expect_error(trackFeatures())
})

#-------------------------------------------------------------------------------

#SPOTCOVERAGE
context("spotCoverage")
sc <- spotCoverage(files = "testfile-Tracks.txt",wells = 24, mm = 5.4, exportdata = T,
                   exportplot = T)

test_that("Is output dataframe?",{
  expect_equal(class(sc),"data.frame")
})

test_that("dataframe columns",{
  expect_equal(all(c("row","plate","well","coverage_pxsq","coverage_mmsq") %in% colnames(sc)),TRUE)
})
rm(sc)

test_that("exportdata",{
  expect_equal("testfile-SpotCoverage.txt" %in% list.files(),TRUE)
})

test_that("exportplot",{
  expect_equal("testfile-SpotCoverage.png" %in% list.files(),TRUE)
})

test_that("Error: no input",{
  expect_error(spotCoverage())
})

#-------------------------------------------------------------------------------

#SPOTALPHAHULL
context("spotAlphahull")
sa <- spotAlphahull(files = "testfile-Tracks.txt",wells = 24, mm = 5.4, alphavalue = 4,
                    exportdata = T, exportplot = T)

test_that("Is output dataframe?",{
  expect_equal(class(sa),"data.frame")
})

test_that("dataframe columns",{
  expect_equal(all(c("row","plate","well","alphahull_pxsq","alphahull_mmsq") %in% colnames(sa)),TRUE)
})
rm(sa)

test_that("exportdata",{
  expect_equal("testfile-SpotAlphahull.txt" %in% list.files(),TRUE)
})

test_that("exportplot",{
  expect_equal("testfile-SpotAlphahull.png" %in% list.files(),TRUE)
})

test_that("Error: no input",{
  expect_error(spotAlphahull())
})

#-------------------------------------------------------------------------------

#SPOTDENSITY
context("spotDensity")
spotDensity(files = "testfile-Tracks.txt",wells = 24, mm = 5.4,
                    exportdata = T, exportplot = T)

test_that("exportplot",{
  expect_equal("testfile-SpotDensity.png" %in% list.files(),TRUE)
})

test_that("Error: no input",{
  expect_error(spotDensity())
})

#-------------------------------------------------------------------------------

#SPOTRATIO
context("spotRatio")
sr <- spotRatio(files = "testfile-Tracks.txt",wells = 24, mm = 5.4,
                    exportdata = T, exportplot = T)

test_that("Is output dataframe?",{
  expect_equal(class(sr),"data.frame")
})

test_that("dataframe columns",{
  expect_equal(all(c("row","plate","well","spots_outer","spots_inner","spots_ratio") %in% colnames(sr)),TRUE)
})
rm(sa)

test_that("exportdata",{
  expect_equal("testfile-SpotRatio.txt" %in% list.files(),TRUE)
})

test_that("exportplot",{
  expect_equal("testfile-SpotRatio.png" %in% list.files(),TRUE)
})

test_that("Error: no input",{
  expect_error(spotRatio())
})

#-------------------------------------------------------------------------------

#REMOVE FILES
if(deleteoutput)
{
  file.remove("testfile-AssignedWells.txt")
  file.remove("testfile-EdgeMarks.txt")
  file.remove("testfile-EdgeMarks.png")
  file.remove("testfile-AssignedWells.png")
  file.remove("testfile-Tracks.txt")
  file.remove("testfile-Tracks.png")
  file.remove("testfile-TrackFeatures.txt")
  file.remove("testfile-SpotCoverage.txt")
  file.remove("testfile-SpotCoverage.png")
  file.remove("testfile-SpotAlphahull.txt")
  file.remove("testfile-SpotAlphahull.png")
  file.remove("testfile-SpotDensity.png")
  file.remove("testfile-SpotRatio.txt")
  file.remove("testfile-SpotRatio.png")
}

#-------------------------------------------------------------------------------

#LTRACK
context("ltrack")
ltrack(files = ll, wells = 24, mm = 5.4, markededges = T,
            framelinking = T, trackfeatures = T, spotcoverage = T,
            spotalphahull = T, spotdensity = T, spotmsd = F, spotratio = T,
                exportdata = T, exportplot = T)

test_that("Does file exist?",{
  expect_equal(file.exists("testfile-AssignedWells.txt"),TRUE)
  expect_equal(file.exists("testfile-EdgeMarks.txt"),TRUE)
  expect_equal(file.exists("testfile-EdgeMarks.png"),TRUE)
  expect_equal(file.exists("testfile-AssignedWells.png"),TRUE)
  expect_equal(file.exists("testfile-Tracks.txt"),TRUE)
  expect_equal(file.exists("testfile-Tracks.png"),TRUE)
  expect_equal(file.exists("testfile-TrackFeatures.txt"),TRUE)
  expect_equal(file.exists("testfile-SpotCoverage.txt"),TRUE)
  expect_equal(file.exists("testfile-SpotCoverage.png"),TRUE)
  expect_equal(file.exists("testfile-SpotAlphahull.txt"),TRUE)
  expect_equal(file.exists("testfile-SpotAlphahull.png"),TRUE)
  expect_equal(file.exists("testfile-SpotDensity.png"),TRUE)
  expect_equal(file.exists("testfile-SpotRatio.txt"),TRUE)
  expect_equal(file.exists("testfile-SpotRatio.png"),TRUE)
  expect_equal(file.exists("testfile-Features.txt"),TRUE)
  expect_equal(file.exists("testfile-TrackFeatures.txt"),TRUE)
})

test_that("Error: no input",{
  expect_error(ltrack())
})

#tf <- read.delim("testfile-Features.txt",header=T)
#colnames(tf)
#-------------------------------------------------------------------------------

if(deleteoutput) file.remove(list.files())
setwd(currwd)
unlink("drtrackerDemo",force=T)


