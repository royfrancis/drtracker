# drtracker 1.0.3  

`drtracker` is an R package to track movement of Zebrafish larvae in 24 or 48 well assay plates from xy data and then compute spot and track features. The features computed are well assignment of spots, minimum spanning distance of spots, convex and concave hulls of spots, linked tracks, track length, track mean speed and track max speed.  

## 1. Installation  

Download and install [R software](https://cran.rstudio.com/) appropriate for your system.  

```r
#load devtools package
library(devtools)
install_github('royfrancis/drtracker')

#load library for use
library(drtracker)
```
## 2. Considerations  

+ Code is designed for 24 or 48 well plates only.  
+ Lighting must induce minimal reflection on plates.  
+ Single larvae per well.  
+ Videos with 500 frames or more work better.  
+ Note the framerate of the video in fps (frames per second).  
+ Calibrate the size in the video to determine 1 pixel = ? mm.  
+ Uncompressed .avi video files are best.  

## 3. Workflow  
The workfow can be divided into *Detection*, *Linking* and *Feature computation*. Detection involves identifying the spots on each frame. Detection can be done using any algorithm and using any software. If you perfrom your own detection to generate xy data, you can skip the following section 3.1.  

For this demonstration, I will use ImageJ for *Detection* to generate xy coordinate data. *Linking* is the process of connecting the xy coordinates frame to frame to create tracks. And the final step is generating track statistics from tracks. *Linking* and *Feature computations* will be done in R. I am running Windows 8.1 64 bit, R 3.2.0 64 bit and Fiji ImageJ 1.50b. I am using 24 well plates.  

### 3.1 Detection

Detection is carried out in ImageJ.  
1. The video must be imported into Fiji ImageJ. If the video is not in .avi format, convert to uncompressed .avi using a tool such as VirtualDub.   
2. Save the video as a .tif sequence and/or continue to next step.  

![Plate view](vignettes/fig1.jpg)  
__Fig 1.__ *View of a typical 24 well plate with few day old single zebrafish larvae in each well.*  

3. Mark the edges of plate (4 spots) and larvae positions on frame 1. Set the brush size to something like size 8-10 black colour. If the well/wells are empty, do not mark anything. When marking the 4 edge spots, imagine connecting the 4 spots to create a rectangle. No larval position must touch that rectangle. If they do, then mark the edges further out.     

![Marked plate](vignettes/fig2.jpg)  
__Fig 2.__ *View of a plate with marking in the first frame. The edges of the wells and the positions of larvae are marked in frame 1.*  

4. Remove background. Go to `Image` > `Stacks` > `Z Project`. Set `Projection type` to `Average Intensity`. Click `OK`. Go to `Process` > `Image Calculator`. Choose the tif stack/sequence as `Image1`, the averaged single image as `Image2` and set `Operation` to `Difference`. Check `Create new window` and `32-bit (float) result`. Click `OK`. In the `Process Stack?` window, click `Yes`.  

![Z stack plate](vignettes/fig3.jpg)  
__Fig 3.__ *Blank plate after Average Z stack. Works best when larvae moves around a lot and number of frames are at least few hundred.*  

![Larvae without plate](vignettes/fig4.jpg)  
__Fig 4.__ *Larval positions without the plate after Image Difference.*  

5. Convert to binary. Go to `Image` > `Adjust` > `Threshold`. Change `Default` to `Max Entropy`. Check `Dark background` and uncheck `Stack histogram`. Click `Apply`. When `NaN Background` windows opens, uncheck `Set background pixels to NaN`. Click `OK`. In the `Convert Stack to binary` window, check `Calculate threshold for each image`, uncheck `Only convert current image`, leave rest as default, then, click `OK`.  

![Larvae thresholded](vignettes/fig5.jpg)  
__Fig 5.__ *Image after thresholding.* 

6. Remove stray spots. Go to `Process` > `Binary` > `Options`. Set `Iterations` to `2`, `Count` to `5`, check `Black background`, uncheck `Pad edges..`, set `EDM output` to `Overwrite` and `Do` to `Erode`. Click `OK`. In the `Process Stack?` window, click `Yes`. This will remove small stray spots.  

Then go to same `Options` and set `Iterations` to `10`, `Count` to `3` and `Do` to `Dilate`. Click `OK`. In the `Process Stack?` window, click `Yes`. This will remove more small stray spots and enlarge the spots.  

The first frame must have exactly 24 larvae spots and the four edge spots. Any extra spots must be removed with the brush tool by painting black.  

![Larvae thresholded](vignettes/fig6.jpg)  
__Fig 6.__ *Image after Erosion/Dilation.* 

7. Analyse spots. Go to `Analyze` > `Set measurements`. Check `Area`, `Centroid` and `Stack position`. Click `OK`.

Go to `Analyze` > `Analyze particles`. `Size` should be `0-Infinity`, `Circularity` as `0.00-1.00` and `Show` `Nothing`. Check `Display results` and uncheck `Exclude on edges`. Click `OK`. In the `Process Stack?` window, click `Yes`.  

8. Save the binary video file if required. Save the `Results` table with suitable filename as .txt. The filename will be used to refer to this file in further analyses.  

### 3.2 Linking and Feature calculation using `drtracker`  

The input file can be prepared in any way or generated from any application. The input file must be a tab-delimited dot (.) as decimal text file. The file must have a minimum of three columns named `x`, `y` and `slice`. We use the text file exported from ImageJ. Start R and load `drtracker` library. Use function `ltrack()`.  

```r
#Load the library
library(drtracker)

#Help and all arguments
?ltrack

#Usage
#Use mm and fps values for your video. Text data and plots are exported by default.
ltrack(choose.files(), wells=24, mm=5.4, fps=25)

#To not export text and plots, use
ltrack(choose.files(), wells=24, mm=5.4, fps=25, exportplot=FALSE, exportdata=FALSE)

#save results (tracks and features) to an R variable for further analyses, plotting etc.  
dframe <- ltrack(choose.files(), wells=24, mm=5.4, fps=25)
dframe <- ltrack(choose.files(), wells=24, mm=5.4, fps=25, exportplot=FALSE, exportdata=FALSE)

#turn off calculation of coverage, msd and alphahulls
ltrack(choose.files(), coverage=FALSE, msd=FALSE, alphahulls=FALSE)
```  

`drtracker` takes xy coordinate data for each frame/slice and computes several variables broadly divided to spot features and track features. Spot features are those computed from the 2D point cloud, while track features are calculated from motion tracks after linking. Spot features are __wells__, __coverage__, __alphahull__ and __msd__. Each spot is assigned to the nearest __well__. The maximum area (convex hull) covered by the spots or maximum larval movement is explained by __coverage__. The area covered by the spots/larvae taking into account gaps/holes in the point cloud (concave hull) is __alphahull__. The minimum distance covered by the larvae is explained by __msd__ (minimum spanning distance).  

Track features are __distance__, __mean speed__, __max speed__, The spots are linked from frame to frame to create tracks. Full length tracks are created for each larvae/well. The total __distance__, __mean speed__ and __max speed__ of each track is computed.  

![Feature plot](vignettes/features.jpg)  
__Fig 7.__ *Features computed by `drtracker` shown for a single well. (A) shows raw spots. (B) shows an mst. (C) shows the coverage (convex hull). (D) shows linked track. (E) and (F) shows Alpha hulls at two values of alpha, 4 and 6 respectively. Smaller values of alpha are more stringent.*  

#### Spot features 

The spots are assigned to wells and the plate layout with well numbers and spots are exported as an image if `exportplot=TRUE`.

![Spot plot](vignettes/fig7.jpg)  
__Fig 8.__ *Layout of the plate with well numbers showing spots in each well.*  

![Spot plot](vignettes/fig8.jpg)  
__Fig 9.__ *Spot plot for a 48-well plate.*  

The minimum spanning distance (msd) is computed if `msd=TRUE`. A plot showing mst lines is exported when `exportplot=TRUE`.  

![Msd plot](vignettes/fig9.jpg)  
__Fig 10.__ *Layout of the plate with spots connected by a minimum spanning tree. The total distance of this tree is taken as the minimum spanning distance.*  

The alphahulls representing area covered by larvae is computed if `alphahull=TRUE`. The default alpha value for alphahulls is set as `alphavalue=4`.  

![Alphahull plot](vignettes/fig10.jpg)  
__Fig 11.__ *Layout of the plate showing the alphahull (concave hull) polygon boundaries and spots are light grey. The area within the polygons are used as a measure of area covered by the larvae.*  

![Track coverage](vignettes/fig11.jpg)  
__Fig 12.__ *Plate layout showing the maximum area (convex hull) covered by larval activity. Number show area in pixels square.*  

#### Track features

The tracks as lines are plotted if `exportplot=TRUE`. If `exportdata=TRUE`, then the spots, tracks and all raw data is exported as a tab-delimited text file.  

![Tracks plot](vignettes/fig12.jpg)  
__Fig 13.__ *Layout of the plate with well numbers showing tracks in each well. The track length in pixels is shown below each track. The tracks are also coloured by track length.*  

The total distance for each track is plotted and exported when `exportplot=TRUE`.

![Tracks distances](vignettes/fig13.jpg)  
__Fig 14.__ *Barplot showing total distance covered by each larvae in pixels.*  

The spot and track features are exported as a text file when `exportdata=T`. 


### 4. Algorithm  
The function `ltrack()` accepts xy coordinates for each spot along with the slice/frame number.  

__Well assignment__  
The number of wells are defined. Every spot is allocated to one of the wells using nearest neighbour search (function `nn2()` from package `RANN`).  

__Convex hull__  
For the area covered by larval activity (coverage), the convex hull is computed from spot data per well and the area of the polygon is calculated. The base function `chull()` is used.  

__Concave hull__  
For the area covered by larval activity taking into account holes and gaps in spots (alphahull), the alpha hull is computed and the area of the resulting polygon is calculated. The alpha hull is computed using function `ahull()` from package `alphahull`.  

__Minimum spanning distance__  
For a measure of minimum distance covered by the larvae, the msd is calculated as the total distance in the minimum spanning tree. The mst is calculated using the function `dino.mst()` from package `fossil`.  

__Frame linking__  
Each spot on frame 1 is assigned an id. Then, each spot is connected from one frame to the next frame. A spot is selected and the algorithm searches for a spot in the next frame in the same well using one of three approaches: *single*, *nearest* or *duplicated*. If a single spot was identified in the next frame (in same well), then *single* is assigned to the point. If more than one point was found in the next frame (same well), then the nearest point is selected and assigned *nearest*. If no point was found in the next frame (same well), then the previous point is duplicated. Once a spot is defined in the next frame, the same id is assigned to that spot. This is iterated to the end of all frames.  

The total distance covered by each id is the sum of distance covered per frame. The speed is calculated as the distance moved per second. A second is defined by framerate. For example, in a 25 fps video, the distance covered every 25 frames is computed and stored. The mean value of all such stored values is the mean speed. The max value of all such stored values is the max speed.  

2015 Roy M Francis | roy.m.francis@outlook.com
