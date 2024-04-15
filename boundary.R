#######################################################################
## purpose: Creating images for removing speckles
## https://dahtah.github.io/imager/imager.html
## https://cran.r-project.org/web/packages/imager/vignettes/gettingstarted.html
#######################################################################
library(imager)
library(x3ptools)
library(ggplot2)
library(dplyr)
library(stats)
library(raster)
library(terra)

x3pA = "scans/FAU 163/Bullet A/LAPD - 163 - Bullet A - Land 6 - Sneox1 - 20x - auto light left image + 20% - threshold 2 - resolution 4 - Carley McConnell.x3p"
maskA = "masks/LAPD-FAU163-BA-L6-R02.png"

x3pB = "scans/FAU 163/Bullet B/LAPD - 163 - Bullet B - Land 1 - Sneox1 - 20x - auto light left image + 20% - threshold 2 - resolution 4 - Carley McConnell.x3p"
maskB = "masks/LAPD-FAU163-BB-L1-R01.png"

# Not so well
x3pC = "scans/FAU 163/Bullet C/LAPD - 163 - Bullet C - Land 2 - Sneox1 - 20x - auto light left image + 20% - threshold 2 - resolution 4 - Carley McConnell.x3p"
maskC = "masks/LAPD-FAU163-BC-L2-R01.png" 

x3pD = "scans/FAU 163/Bullet D/LAPD - 163 - Bullet D - Land 3 - Sneox1 - 20x - auto light left image + 20% - threshold 2 - resolution 4 - Carley McConnell.x3p"
maskD = "masks/LAPD-FAU163-BD-L3-R01.png"

land <- x3p_read(x3pC)
mask <- png::readPNG(maskC)
overlay <- x3p_add_mask(land, mask = as.raster(mask))
down <- overlay %>% x3p_sample(m=20)
down %>% x3p_image()  
rgl::rglwidget()

boundaryW = imfill(x=down$header.info$sizeX, y=down$header.info$sizeY, z=1, val=1)
boundaryW[is.na(down$surface.matrix)]=0
plot(boundaryW)

boundaryB = imfill(x=down$header.info$sizeX, y=down$header.info$sizeY, z=1, val=0)
boundaryB[is.na(down$surface.matrix)]=1
plot(boundaryB)


outRegion(boundaryW, px= 4)$px
outRegion(boundaryW)$px


# threshold not working, because what we have is binarized values
threshold(boundaryW, "15%") %>% plot
boundaryW %>% as.data.frame %>% head() 

# blur with guassian filter to vary the values
# blur and threshold, not working because the boundaries are changes
boundaryW %>% isoblur(1) %>% head() %>% as.data.frame() %>% threshold("10%")
boundaryW %>% isoblur(1,neumann = FALSE)  %>% threshold("5%") %>% plot
boundaryW %>% isoblur(1) %>% head() %>% as.data.frame() %>% threshold("10%")
boundaryW %>% isoblur(1,neumann = FALSE)  %>% threshold("15%") %>% plot

# Use blob, not working also because the altered boundaries
# Could be helpful combining with morphological operations?
# https://dahtah.github.io/imager/morphology.html
df = with(imhessian(boundaryW), (xx*yy-xy^2)) %>% threshold('99%') %>% label()%>%as.data.frame() %>% filter(value>0)
df %>% group_by(value) %>% summarise(nx=n())

blur = boundaryW %>% isoblur(1) 
plot(blur)
thresh2 = with(imhessian(blur), (xx*yy-xy^2)) %>% threshold('99%') %>% label()
plot(thresh2)
df2 = with(imhessian(blur), (xx*yy-xy^2)) %>% threshold('99%') %>% label()%>%as.data.frame() %>% filter(value>0)
df2 %>% group_by(value) %>% summarise(n=n())
center2 = df2 %>% group_by(value) %>% summarise(mx=mean(x), my=mean(y))
plot(thresh2)
with(center2, points(mx,my, col='red'))

# blob with scale, I don't quite get how this is done
library(purrr)
hessdet <- function(im,scale=1) isoblur(im,scale) %>% imhessian %$% { scale^2*(xx*yy - xy^2) }
dat <- map_df(2:4,function(scale) hessdet(boundaryW,scale) %>% as.data.frame %>% mutate(scale=scale))
scales <- seq(2,20,l=10)

d.max <- map_il(scales,function(scale) hessdet(boundaryW,scale)) %>% parmax
plot(d.max,main="Point-wise maximum across scales")


# Back to morphological operations
plot(boundaryW)
# grow as dilation, increase white dot
boundaryW %>% grow(2) %>% plot
# shrink remove isolated white dots
boundaryW %>% shrink(2) %>% plot

# dilate and erode to remove while holes
# Remove holes, not so much for separating regions
boundaryW %>% clean(2) %>% plot
boundaryW %>% shrink(2) %>% grow(2) %>% plot

# helps creating boundaries to the regions while not hurting the initially defined outside regions
boundaryW %>% grow(2) %>% shrink(2) %>% plot
boundaryW %>% fill(2) %>% plot


outRegion <- function(cimg, maxY=1.5, px = NA){
  # A helper function to determine the number of pixels of highlight regions
  numElement <- function(element){
    max(lengths(element))
  }
  
  if (is.na(px) == TRUE){ # automatically find the lines
    px = 1
    find = 0
    while (find == 0 & px <=5) {
      highlight = boundaryW %>% fill(px) %>% highlight
      numHighlights = sapply(highlight, numElement)
      # The top 2 number of pixels of highlight regions
      numMax = sort(numHighlights, decreasing = TRUE)[1:3]
      
      if (sum(numMax > dim(boundaryW)[2])==2){# top 2 has all y axis covered
        if (sum(numMax[1:2] <= maxY*dim(cimg)[2]) <= 2){
          ms = order(numHighlights, decreasing = TRUE)[1:2]
          find = 1
          h1 = highlight[[ms[1]]]
          h2 = highlight[[ms[2]]]
        }
        else{
          find = 0}
      }
      else{
        find = 0
      }
      px = px+1
    }
  }
  
  else{
    find=1
    highlight = boundaryW %>% fill(px) %>% highlight
    numHighlights = sapply(highlight, numElement)
    # The top 2 number of pixels of highlight regions
    numMax = sort(numHighlights, decreasing = TRUE)[1:3]
    ms = order(numHighlights, decreasing = TRUE)[1:2]
    h1 = highlight[[ms[1]]]
    h2 = highlight[[ms[2]]]
  }
  
  if(find==1){
    plot(cimg)
    with(h1, lines(x,y, col='red', lwd=2))
    with(h2, lines(x,y, col='blue', lwd=2))
    return(list(px=px, h1=data.frame(h1$x, h1$y), h2=data.frame(h2$x, h2$y)))}
  else{
    return("Not available for px<=5")
  }
}


# get outlines
plot(boundaryW)
highlight = boundaryW %>% fill(3) %>% highlight

# Getting the subsets of highlighted regions
numElement <- function(element){
  max(lengths(element))
}
numHighlights = sapply(highlight, numElement) # the maximum has 
# The largest number of pixels of highlight regions
numMax = sort(numHighlights, decreasing = TRUE)[1:3] 
#  LIMIT THE NUMMAX HERE LARGER THAN NUM OF Y BUT LESS THAN A NUMBER
# AS A CONDITION?
numMax > dim(boundaryW)[2]
numMax < 80
ms = order(numHighlights, decreasing = TRUE)[1:2]
# two data sets from highlighted region helps seperating the image
h1 = highlight[[ms[1]]]
h2 = highlight[[ms[2]]]
plot(boundaryW)
with(h1, lines(x,y, col='green', lwd=2))
with(h2, lines(x,y, col='green', lwd=2))
