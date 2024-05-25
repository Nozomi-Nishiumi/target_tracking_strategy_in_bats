#  Created by nozomi nishiumi(n.oz.o@hotmail.co.jp) on 2022/12/01.
if( dev.cur() > 1 ) dev.off()
rm( list = ls( envir = globalenv() ), envir = globalenv() )

targetPackages <- c("vctrs",
                    "ggplot2", 
                    "MASS",
                    "mgcv",
                    "rgl",
                    "pracma",
                    "lme4",
                    "lmerTest",
                    "tcltk",
                    "MuMIn",
                    "car",
                    "merTools",
                    "performance",
                    "ggExtra",
                    "SimDesign",
                    "multcompView",
                    "multcomp",
                    "here"
) 
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)


rm(list = ls(all.names = TRUE))

red=hsv(355/360, 49/100, 99/100, alpha=1)
green=hsv(100/360, 100/100, 89/100, alpha=1)
blue=hsv(222/360, 72/100, 100/100, alpha=1)

red2=hsv(355/360, 59/100, 80/100, alpha=0.2)
green2=hsv(100/360, 100/100, 70/100, alpha=0.2)
blue2=hsv(222/360, 72/100, 100/100, alpha=0.1)

vectoring<-function(HV){
  h=HV[,1]/180*pi
  v=HV[,2]/180*pi
  x=cos(v)
  y=sin(v)*sin(h)
  z=sin(v)*cos(h)
  E=data.frame(x,y,z)
  return(E)
}

direction2<-function(target, base=c(0,0,0)){
  R=target-base
  colnames(R)=c("x","y","z")
  h=atan2(R[,2],R[,3])
  v=atan2(sqrt(R[,2]^2+R[,3]^2),R[,1])
  dist=sqrt(R[,1]^2+R[,2]^2+R[,3]^2)
  
  
  x=cos(v)
  y=sin(v)*sin(h)
  z=sin(v)*cos(h)
  E=data.frame(x,y,z)
  
  H=h*180/pi
  V=v*180/pi
  H=normalizeRotation(H)
  
  
  LeadAngle=atan2(z,x)*180/pi
  E2=rotationPolar(E,atan2(z,x),2)
  SideAngle=atan2(E2$y,E2$x)*180/pi
  Longi=atan2(z,sqrt(x^2+y^2))*180/pi
  
  directionData<-list(data.frame(H,V,LeadAngle,SideAngle,Longi),data.frame(x=E$x,y=E$y,z=E$z),data.frame(x=R$x,y=R$y,z=R$z))  
  names(directionData)=c("hv","xyz","xyz_R")
  return(directionData)
}

cnvrtVector2LongiLatiAngle<-function(v){
  x=v[,1]
  y=v[,2]
  z=v[,3]
  
  Longi=atan2(z,sqrt(x^2+y^2))*180/pi
  Latit=atan2(y,x)*180/pi
  
  return(data.frame(Longi,Latit))
}

cnvrtAngle2Vector<-function(angle){
  
  H=angle[,1]
  V=angle[,2]
  x=cos(V*(pi/180))*cos(H*(-pi/180))
  y=cos(V*(pi/180))*sin(H*(-pi/180))
  z=sin(V*(pi/180))
  
  
  aZ=atan2(y,x)
  rotatedZ=data.frame(cos(aZ-aZ)*sqrt(x^2+y^2),sin(aZ-aZ)*sqrt(x^2+y^2),z)
  colnames(rotatedZ)=c("x","y","z")
  
  
  aY=atan2(rotatedZ$z,rotatedZ$x)
  rotatedY=data.frame(cos(aY-aY)*sqrt(rotatedZ$x^2+rotatedZ$z^2),rotatedZ$y,sin(aY-aY)*sqrt(rotatedZ$x^2+rotatedZ$z^2))
  colnames(rotatedY)=c("x","y","z")
  
  aY_90=aY-pi/2
  
  a<-data.frame(H,V,x,y,z,aZ,aY,aY_90)
  return(a)
}

angleForAlign<-function(vector){
  
  
  x=vector$x
  y=vector$y
  z=vector$z
  
  
  aZ=atan2(y,x)
  rotatedZ=data.frame(cos(aZ-aZ)*sqrt(x^2+y^2),sin(aZ-aZ)*sqrt(x^2+y^2),z)
  colnames(rotatedZ)=c("x","y","z")
  
  
  aY=atan2(rotatedZ$z,rotatedZ$x)
  rotatedY=data.frame(cos(aY-aY)*sqrt(rotatedZ$x^2+rotatedZ$z^2),rotatedZ$y,sin(aY-aY)*sqrt(rotatedZ$x^2+rotatedZ$z^2))
  colnames(rotatedY)=c("x","y","z")
  
  
  a<-data.frame(aZ,aY)
  colnames(a)=c("Z","Y")
  return(a)
}

rotationPolar<-function(vector,angle,xyz){
  
  x=vector$x
  y=vector$y
  z=vector$z
  
  aX=atan2(y,z)
  aY=atan2(z,x)
  aZ=atan2(y,x)
  
  if(xyz==1){
    rotated=data.frame(x,sin(aX-angle)*sqrt(z^2+y^2),cos(aX-angle)*sqrt(z^2+y^2))
  }
  if(xyz==2){
    rotated=data.frame(cos(aY-angle)*sqrt(z^2+x^2),y,sin(aY-angle)*sqrt(z^2+x^2))
  }
  if(xyz==3){
    rotated=data.frame(cos(aZ-angle)*sqrt(x^2+y^2),sin(aZ-angle)*sqrt(x^2+y^2),z)
  }
  
  colnames(rotated)=c("x","y","z")
  
  return(rotated)
}

normalizeRotation<-function(x){
  while(length(x[which(x>=360)])>0){
    x[which(x>=360)]=x[which(x>=360)]-360
  }
  
  while(length(x[which(x<0)])>0){
    
    x[which(x<0)]=x[which(x<0)]+360
  }
  
  return(x)
}

drawSphere<-function(r=1,size=180){
  
  a=pi*size/180
  theta=seq(0,2*pi+0.1,0.1)
  x=r*cos(theta)
  y=r*sin(theta)
  z=0*theta
  circle=data.frame(x,y,z)
  plot3d(circle,type="n",axes=F,
         xlab=c(""),ylab=c(""),zlab=c(""),
         xlim=c(-1.5*r,1.5*r),ylim=c(-1.5*r,1.5*r),zlim=c(-1.5*r,1.5*r))
  aspect3d(1,1,1)
  bearing=seq(0,pi,pi/8)
  for(i in 1:length(bearing)){
    longitude.x=circle$x
    longitude.y=circle$y*cos(bearing[i])
    longitude.z=circle$y*sin(bearing[i])
    longitude=data.frame(longitude.x,longitude.y,longitude.z)
    rgl.linestrips(longitude,color="black",alpha=0.3,lwd=1.0)
  }
  
  circle2=data.frame(z,x,y)
  colnames(circle2)=c("x","y","z")
  
  deviation=seq(0,a,a/4)
  for(i in 1:length(deviation)){
    latitude.r=r*sin(deviation[i])
    latitude.x=circle2$x+r*cos(deviation[i])
    latitude.y=circle2$y*latitude.r/r
    latitude.z=circle2$z*latitude.r/r
    latitude=data.frame(latitude.x,latitude.y,latitude.z)
    rgl.linestrips(latitude,color="black",alpha=0.3,lwd=1.0)
  }
  
  text.longitude=seq(0,2*pi-pi/8,pi/2)
  info.longitude=c("0","-90","180/-180","90")
  deviation2=deviation-a/8
  for(i in 1:length(text.longitude)){
    text.r=r*1.1
    text.r_longitude=r*sin(deviation2[3])*1.1
    text.x=text.r*cos(deviation2[3])
    text.y=text.r_longitude*sin(text.longitude[i])
    text.z=text.r_longitude*cos(text.longitude[i])
    text3d(text.x, text.y,text.z, texts=info.longitude[i], color="black")
  }
  
  info.latitude=seq(0,size,size/4)
  for(i in 2:length(deviation)){
    text.r=r*1.1
    text.r_latitude=r*sin(deviation[i])
    text.x=text.r*cos(deviation[i])
    text.y=text.r_latitude*sin(pi/16)
    text.z=text.r_latitude*cos(pi/16)
    text3d(text.x, text.y,text.z, texts=info.latitude[i], color="black")
  }
  
  deviation=seq(pi,a,-pi/4)
  for(i in 1:length(deviation)){
    latitude.r=r*sin(deviation[i])
    latitude.x=circle2$x+r*cos(deviation[i])
    latitude.y=circle2$y*latitude.r/r
    latitude.z=circle2$z*latitude.r/r
    latitude=data.frame(latitude.x,latitude.y,latitude.z)
    rgl.linestrips(latitude,color="black",alpha=0.3,lwd=1.0)
  }
  
  labels=c("x","y","z")
  x=c(r*1.2,0,0)
  y=c(0,r*1.2,0)
  z=c(0,0,r*1.2)
  axes=data.frame(x,y,z)
  rgl.texts(axes, text=labels, color="black")
}

scaling_polar<-function(vector,inputangle,targetangle){
  
  x=vector$x
  y=vector$y
  z=vector$z
  
  aX=atan2(y,z)
  rotated=rotationPolar(vector,aX,1)
  
  one=pi/180
  calib=one*targetangle/inputangle
  
  angle=atan2(rotated$z,rotated$x)
  scaled<-data.frame(cos(angle*calib)*sqrt(rotated$z^2+rotated$x^2),rotated$y,sin(angle*calib)*sqrt(rotated$z^2+rotated$x^2))
  colnames(scaled)=c("x","y","z")
  
  Re_rotated=rotationPolar(scaled,-aX,1)
  return(Re_rotated)
}

Polar_line2<-function(data1,data2,data3,size,limit,xlab,ylab,pointON,pathON,densityON, textON,backON,title,polar=1,hist="", area=NULL){
  
  if(!is.null(data1)){
    longitude1=cnvrtVector2LongiLatiAngle(vectoring(data1))[,1]
    latitude1=-cnvrtVector2LongiLatiAngle(vectoring(data1))[,2]
    dataset1=na.omit(as.data.frame(cbind(latitude1,longitude1)))
  }else{dataset1=NULL}
  
  if(!is.null(data2)){
    longitude2=cnvrtVector2LongiLatiAngle(vectoring(data2))[,1]
    latitude2=-cnvrtVector2LongiLatiAngle(vectoring(data2))[,2]
    dataset2=na.omit(as.data.frame(cbind(latitude2,longitude2)))
  }else{dataset2=NULL}
  
  if(!is.null(data3)){
    longitude3=cnvrtVector2LongiLatiAngle(vectoring(data3))[,1]
    latitude3=-cnvrtVector2LongiLatiAngle(vectoring(data3))[,2]
    dataset3=na.omit(as.data.frame(cbind(latitude3,longitude3)))
  }else{dataset3=NULL}
  
  
  
  draw<-ggplot() 
  if(densityON=="densityON_b"){
    draw<-draw+
      stat_density2d(data=dataset2[,1:2], aes(x=dataset2[,1:2][,1], y=dataset2[,1:2][,2], alpha=..density..), 
                     fill=rgb(1,0.4,0.4,alpha=0.8), geom="raster", contour = "false")
  }
  
  theta=seq(0,2*pi, by=pi/100)
  circlex=cos(theta)*limit
  circley=sin(theta)*limit
  
  if(polar){
    draw<-draw+
      geom_polygon(aes(x=circlex, y=circley), colour="black", fill=NA,size = 0.4)+
      geom_polygon(aes(x=circlex*2/4, y=circley*2/4), colour="black", fill=NA,size = 0.1)+
      geom_line(aes(x=c(-limit,limit),y=c(0,0)), colour="black",size = 0.4)+
      geom_line(aes(x=c(0,0),y=c(-limit,limit)), colour="black",size = 0.4)
  }
  
  if(pointON=="pointON_a"){
    draw<-draw+
      geom_point(data=dataset1[,1:2], aes(x=dataset1[,1:2][,1], y=dataset1[,1:2][,2]),
                 shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = rgb(0, 0, 0, alpha=0.5))+
      geom_point(data=dataset2[,1:2], aes(x=dataset2[,1:2][,1], y=dataset2[,1:2][,2]),
                 shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = rgb(0, 0, 0, alpha=0.5))+
      geom_point(data=dataset3[,1:2], aes(x=dataset3[,1:2][,1], y=dataset3[,1:2][,2]),
                 shape = 21, size = size*0.8, stroke = 1, colour = rgb(1, 0.35, 0.4, alpha=1), fill = rgb(1, 0.35, 0.4, alpha=0.5))
    
  }
  
  if(pointON=="pointON_b"){
    draw<-draw+
      geom_point(data=dataset1[,1:2], aes(x=dataset1[,1], y=dataset1[,2]),
                 shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = rgb(0.5, 1, 0.5, alpha=0.6))+
      geom_point(data=dataset2[,1:2], aes(x=dataset2[,1], y=dataset2[,2]),
                 shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = hsv(4/360, 0/100, 80/100, alpha=0.6))+
      geom_point(data=dataset3[,1:2], aes(x=dataset3[,1], y=dataset3[,2]),
                 shape = 21, size = size, stroke = 0.3, colour = rgb(1, 0.35, 0.4, alpha=1), fill = rgb(1, 0.35, 0.4, alpha=0.5))
    
  }
  
  if(pointON=="pointON_c"){
    draw<-draw+
      geom_point(data=dataset2[,1:2], aes(x=dataset2[,1:2][,1], y=dataset2[,1:2][,2]),
                 shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.2), fill = rgb(0, 0, 0, alpha=0.5))
  }
  
  if(pathON=="pathON"){
    draw<-draw+
      geom_path(data=dataset1[,1:2], aes(x=dataset1[,1:2][,1], y=dataset1[,1:2][,2]),
                lwd = 1,col = rgb(0.5, 1, 0.5,alpha=0.3))+
      geom_path(data=dataset2[,1:2], aes(x=dataset2[,1:2][,1], y=dataset2[,1:2][,2]),
                lwd = 1,col = rgb(0, 0, 0,alpha=0.3))+
      geom_path(data=dataset3[,1:2], aes(x=dataset3[,1:2][,1], y=dataset3[,1:2][,2]),
                lwd = 1,col = rgb(1, 0.35, 0.4, alpha=0.3))
  }
  
  if(densityON=="densityON_a"){
    draw<-draw+
      geom_density2d(data=dataset2[,1:2], aes(x=dataset2[,1:2][,1], y=dataset2[,1:2][,2]),col = rgb(1, 0,0),size=0.5)
  }
  
  
  
  if(textON=="textON"){
    sample1=seq(1,nrow(dataset1),6)
    sample2=seq(1,nrow(dataset2),6)
    sample3=seq(1,nrow(dataset3),6)
    
    omit1=dataset1[sample1,]
    omit2=dataset2[sample2,]
    omit3=dataset3[sample3,]
    
    print((dataset1))
    
    text1<-geom_text(data=omit1[,1:2],
                     aes(x=omit1[,1], y=omit1[,2], label=sample1-1),cex=7,fontface = "bold")
    
    text2<-geom_text(data=omit2[,1:2],
                     aes(x=omit2[,1], y=omit2[,2], label=sample2-1),cex =7,fontface = "bold",col=hsv(4/360, 65/100, 10/100, alpha=1))
    
    text3<-geom_text(data=omit3[,1:2],
                     aes(x=omit3[,1], y=omit3[,2], label=sample3-1),cex =7,fontface = "bold",col = hsv(222/360, 60/100, 10/100, alpha=1))
    
    text1_circle<-geom_point(data=omit1[,1:2],
                             aes(x=omit1[,1], y=omit1[,2], label=sample1-1),
                             shape = 21, size = size, stroke = 1, colour = hsv(0, 0, 30/100, alpha=1), fill = rgb(1, 0.35, 0.4, alpha=0))
    text2_circle<-geom_point(data=omit2[,1:2],
                             aes(x=omit2[,1], y=omit2[,2], label=sample2-1),
                             shape = 21, size = size, stroke = 1, colour = hsv(0, 0, 30/100, alpha=1), fill = hsv(222/360, 60/100, 100/100, alpha=0))
    text3_circle<-geom_point(data=omit3[,1:2],
                             aes(x=omit3[,1], y=omit3[,2], label=sample3-1),
                             shape = 21, size = size, stroke = 1, colour = hsv(0, 0, 30/100, alpha=1), fill = hsv(222/360, 60/100, 100/100, alpha=0))
    
    draw<-draw+text2_circle+text3_circle
  }
  
  if(backON=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA), 
      axis.text= element_blank(),
      plot.title = element_blank(),
      axis.title = element_blank()
    )
  }
  
  
  
  if(pointON=="pointON_c"){
    draw<-draw+
      geom_point(data=dataset1[,1:2], aes(x=dataset1[,1:2][,1], y=dataset1[,1:2][,2]),
                 shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = rgb(0, 0, 0, alpha=0.5))+
      geom_point(data=dataset3[,1:2], aes(x=dataset3[,1:2][,1], y=dataset3[,1:2][,2]),
                 shape = 21, size = size*0.8, stroke = 1, colour = rgb(1, 0.35, 0.4, alpha=1), fill = rgb(1, 0.35, 0.4, alpha=0.5))
    
  }
  
  if(!is.null(area)){
    draw<-draw+xlim(-area,area)+ylim(-area,area)+coord_fixed()
  }else{
    draw<-draw+coord_fixed()
  }
  draw<-draw+xlab(xlab)+ylab(ylab)
  
  if(hist=="histON"){
    
    draw<-ggMarginal(
      draw,
      binwidth = 0.5,
      type = "histogram",
      margins = "both",
      size = 7,
      col=hsv(0,0,0.5),
      fill=hsv(0,0,0.5)
    )
  }
  
  
  print(draw)
  if(!is.null(title)){
    ggsave(title, draw, bg = "transparent",dpi=dpi1,width=10,height=10)
  }
  return(draw)
}

Polar_line3<-function(data1,data2,data3,size,limit,xlab,ylab,pointON,pathON,densityON, textON,backON,title,polar=1,hist="", area=NULL){
  
  if(!is.null(data1)){
    longitude1=cnvrtVector2LongiLatiAngle(vectoring(data1))[,1]
    latitude1=-cnvrtVector2LongiLatiAngle(vectoring(data1))[,2]
    dataset1=na.omit(as.data.frame(cbind(latitude1,longitude1)))
  }else{dataset1=NULL}
  
  if(!is.null(data2)){
    longitude2=cnvrtVector2LongiLatiAngle(vectoring(data2))[,1]
    latitude2=-cnvrtVector2LongiLatiAngle(vectoring(data2))[,2]
    dataset2=na.omit(as.data.frame(cbind(latitude2,longitude2)))
  }else{dataset2=NULL}
  
  if(!is.null(data3)){
    longitude3=cnvrtVector2LongiLatiAngle(vectoring(data3))[,1]
    latitude3=-cnvrtVector2LongiLatiAngle(vectoring(data3))[,2]
    dataset3=na.omit(as.data.frame(cbind(latitude3,longitude3)))
  }else{dataset3=NULL}
  
  
  
  draw<-ggplot() 
  theta=seq(0,2*pi, by=pi/100)
  circlex=cos(theta)*limit
  circley=sin(theta)*limit
  
  if(polar){
    draw<-draw+
      geom_polygon(aes(x=circlex, y=circley), colour="black", fill=NA,size = 0.4)+
      geom_polygon(aes(x=circlex*2/4, y=circley*2/4), colour="black", fill=NA,size = 0.1)+
      geom_line(aes(x=c(-limit,limit),y=c(0,0)), colour="black",size = 0.4)+
      geom_line(aes(x=c(0,0),y=c(-limit,limit)), colour="black",size = 0.4)
  }
  
  if(pointON=="pointON_b"){
    draw<-draw+
      geom_point(data=dataset2[,1:2], aes(x=dataset2[,1], y=dataset2[,2]),
                 shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = hsv(4/360, 0/100, 80/100, alpha=0.6))+
      geom_point(data=dataset3[,1:2], aes(x=dataset3[,1], y=dataset3[,2]),
                 shape = 21, size = size, stroke = 0.3, colour = rgb(1, 0.35, 0.4, alpha=1), fill = rgb(1, 0.35, 0.4, alpha=0.5))
    
  }
  
  if(is.null(data1)){
    draw<-draw+
      geom_path(data=dataset3[,1:2], aes(x=dataset3[,1:2][,1], y=dataset3[,1:2][,2]),
                lwd = 1,col = rgb(1, 0.35, 0.4, alpha=0.3))+
      geom_path(data=dataset2[,1:2], aes(x=dataset2[,1:2][,1], y=dataset2[,1:2][,2]),
                lwd = 1,col = rgb(0, 0, 0,alpha=0.3))
  }else{
    draw<-draw+
      geom_path(data=dataset1[,1:2], aes(x=dataset1[,1:2][,1], y=dataset1[,1:2][,2]),
                lwd = 1,col = rgb(1, 0.35, 0.4, alpha=0.3))+
      geom_path(data=dataset2[,1:2], aes(x=dataset2[,1:2][,1], y=dataset2[,1:2][,2]),
                lwd = 1,col = rgb(0, 0, 0,alpha=0.3))
  }
  
  if(backON=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA), 
      axis.text= element_blank(),
      plot.title = element_blank(),
      axis.title = element_blank()
    )
  }
  
  
  
  
  if(!is.null(area)){
    draw<-draw+xlim(-area,area)+ylim(-area,area)+coord_fixed()
  }else{
    draw<-draw+coord_fixed()
  }
  draw<-draw+xlab(xlab)+ylab(ylab)
  
  
  print(draw)
  ggsave(title, draw, bg = "transparent",dpi=dpi1,width=10,height=10)
}

olthoPlot<-function(v1,v2,p1,p2,angle,size,xlab,ylab,textON,backOFF,title){
  point1=p1
  point2=p2
  
  draw<-ggplot()
  pulse_vector<-geom_spoke(data=p1[,1:2],aes(x=p1[,1],y=p1[,2],angle = angle*pi/180, radius = 0.200),
                           lwd = 0.5, col = rgb(0, 0, 0, alpha=0.7))
  point_a<-geom_point(data=point1[,1:2], aes(x=point1[,1:2][,1], y=point1[,1:2][,2]),
                      shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = hsv(222/360, 60/100, 100/100, alpha=0.6))
  point_b<-geom_point(data=point2[,1:2], aes(x=point2[,1:2][,1], y=point2[,1:2][,2]),
                      shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = rgb(1, 0.35, 0.4, alpha=0.5))
  
  point_a_init<-geom_point(aes(x=p1[1,1], y=p1[1,2]),
                           shape = 22, size = size, stroke = 1, colour = rgb(0, 0, 0), fill = hsv(222/360, 60/100, 100/100, alpha=1))
  point_a_end<-geom_point(aes(x=p1[nrow(p1),1], y=p1[nrow(p1),2]),
                          shape = 23, size = size, stroke = 1, colour = rgb(0, 0, 0), fill = hsv(222/360, 60/100, 100/100, alpha=1))
  
  point_b_init<-geom_point(aes(x=p2[1,1], y=p2[1,2]),
                           shape = 22, size = size, stroke = 1, colour = rgb(0, 0, 0), fill = rgb(0.5,0.5,0.5, alpha=1))
  point_b_end<-geom_point(aes(x=p2[nrow(p2),1], y=p2[nrow(p2),2]),
                          shape = 23, size = size, stroke = 1, colour = rgb(0, 0, 0), fill = rgb(0.5,0.5,0.5, alpha=1))
  
  
  path_a<-geom_path(data=v1[,1:2], aes(x=v1[,1:2][,1], y=v1[,1:2][,2]),
                    lwd = 1,col = hsv(222/360, 60/100, 100/100,alpha=0.3))
  path_b<-geom_path(data=v2[,1:2], aes(x=v2[,1:2][,1], y=v2[,1:2][,2]),
                    lwd = 1,col = rgb(1, 0.35, 0.4,alpha=0.3))
  
  
  
  if(backOFF=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA), 
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text= element_blank()
    )
  }
  
  xmin=min(na.omit(c(v1[,1],v2[,1])))
  xmax=max(na.omit(c(v1[,1],v2[,1])))
  ymin=min(na.omit(c(v1[,2],v2[,2])))
  ymax=max(na.omit(c(v1[,2],v2[,2])))
  
  xrange=xmax-xmin
  yrange=ymax-ymin
  aspect_ratio=yrange/xrange
  
  xlimitation<-xlim(xmin-xrange*0,xmax+xrange*0.0)
  ylimitation<-ylim(ymin-yrange*0.8,ymax+yrange*0.8)
  
  draw<-draw+
    point_a+point_b+
    path_a+path_b+
    # point_a_init+point_b_init+
    # point_a_end+point_b_end+
    pulse_vector
  
  
  if(textON=="textON"){
    
    sample1=seq(1,nrow(p1),6)
    sample2=seq(1,nrow(p2),6)
    
    omit1=p1[sample1,]
    omit2=p2[sample2,]
    
    
    text1<-geom_text(data=omit1[,1:2],
                     aes(x=omit1[,1], y=omit1[,2]+0.150, label=sample1-1),cex=3,fontface = "bold")
    
    text2<-geom_text(data=omit2[,1:2],
                     aes(x=omit2[,1], y=omit2[,2]+0.150, label=sample2-1),cex =3,fontface = "bold",colour=rgb(0, 0, 0.6))
    
    
    text1_circle<-geom_point(data=omit1[,1:2],
                             aes(x=omit1[,1], y=omit1[,2], label=sample1-1),
                             shape = 21, size = size, stroke = 1, colour = rgb(0, 0, 0, alpha=0.6), fill = rgb(1, 0.35, 0.4, alpha=0))
    text2_circle<-geom_point(data=omit2[,1:2],
                             aes(x=omit2[,1], y=omit2[,2], label=sample2-1),
                             shape = 21, size = size, stroke = 1, colour = rgb(0, 0, 0, alpha=0.6), fill = hsv(222/360, 60/100, 100/100, alpha=0))
    
    
    # text1_init<-geom_text(aes(x=dataset1[1,1], y=dataset1[1,2], label="i"),cex =6,fontface = "bold")
    # text2_init<-geom_text(aes(x=dataset2[1,1], y=dataset2[1,2], label="i"),cex =6,fontface = "bold")
    # 
    # text1_end<-geom_text(aes(x=dataset1[nrow(dataset1),1], y=dataset1[nrow(dataset1),2], label="e"),cex =6,fontface = "bold")
    # text2_end<-geom_text(aes(x=dataset2[nrow(dataset2),1], y=dataset2[nrow(dataset2),2], label="e"),cex =6,fontface = "bold")
    draw<-draw+
      # text1+
      # text2+
      text1_circle+
      text2_circle
  }
  
  draw<-draw+
    xlimitation+ylimitation+coord_fixed()
  draw<-draw+xlab(xlab)+ylab(ylab)
  print(draw)
  # theme(aspect.ratio=aspect_ratio)
  ggsave(title, draw, bg = "transparent",dpi=dpi1,width=10,height=10) #width and height were previously 5
  
  
}

olthoPlot_sim<-function(v1,v2,v3,p1,p2,p3,size,xlab, ylab,textON,backOFF,title=""){
  point1=p1
  point2=p2
  point3=p3
  
  
  draw<-ggplot()
  
  point_a<-geom_point(data=point1[,1:2], aes(x=point1[,1:2][,1], y=point1[,1:2][,2]),
                      shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = hsv(222/360, 60/100, 100/100, alpha=0.6))
  point_b<-geom_point(data=point2[,1:2], aes(x=point2[,1:2][,1], y=point2[,1:2][,2]),
                      shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = rgb(1, 0.35, 0.4, alpha=0.5))
  point_c<-geom_point(data=point3[,1:2], aes(x=point3[,1:2][,1], y=point3[,1:2][,2]),
                      shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = hsv(100/360, 0/100, 30/100, alpha=0.6))
  
  
  path_a<-geom_path(data=v1[,1:2], aes(x=v1[,1:2][,1], y=v1[,1:2][,2]),
                    lwd = 1,col = hsv(222/360, 60/100, 100/100,alpha=0.3))
  path_b<-geom_path(data=v2[,1:2], aes(x=v2[,1:2][,1], y=v2[,1:2][,2]),
                    lwd = 1,col = rgb(1, 0.35, 0.4,alpha=0.3))
  path_c<-geom_path(data=v3[,1:2], aes(x=v3[,1:2][,1], y=v3[,1:2][,2]),
                    lwd = 1,col = hsv(100/360, 0/100, 20/100,alpha=0.3))
  
  
  
  if(backOFF=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA), 
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text= element_blank()
    )
  }
  
  xmin=min(na.omit(c(v1[,1],v2[,1],v3[,1])))
  xmax=max(na.omit(c(v1[,1],v2[,1],v3[,1])))
  ymin=min(na.omit(c(v1[,2],v2[,2],v3[,2])))
  ymax=max(na.omit(c(v1[,2],v2[,2],v3[,2])))
  
  xrange=xmax-xmin
  yrange=ymax-ymin
  aspect_ratio=yrange/xrange
  
  xlimitation<-xlim(xmin-xrange*0,xmax+xrange*0.0)
  ylimitation<-ylim(ymin-yrange*0.8,ymax+yrange*0.8)
  
  draw<-draw+
    point_a+point_b+point_c+
    path_a+path_b+path_c
  
  
  if(textON=="textON"){
    
    sample1=seq(1,nrow(p1),6)
    sample2=seq(1,nrow(p2),6)
    sample3=seq(1,nrow(p3),6)
    
    omit1=p1[sample1,]
    omit2=p2[sample2,]
    omit3=p3[sample3,]
    
    
    text1<-geom_text(data=omit1[,1:2],
                     aes(x=omit1[,1], y=omit1[,2]+150, label=sample1-1),cex=3,fontface = "bold")
    
    text2<-geom_text(data=omit2[,1:2],
                     aes(x=omit2[,1], y=omit2[,2]+150, label=sample2-1),cex =3,fontface = "bold",colour=rgb(0, 0, 0.6))
    
    
    text1_circle<-geom_point(data=omit1[,1:2],
                             aes(x=omit1[,1], y=omit1[,2], label=sample1-1),
                             shape = 21, size = size, stroke = 1, colour = rgb(0, 0, 0, alpha=0.6), fill = rgb(1, 0.35, 0.4, alpha=0))
    text2_circle<-geom_point(data=omit2[,1:2],
                             aes(x=omit2[,1], y=omit2[,2], label=sample2-1),
                             shape = 21, size = size, stroke = 1, colour = rgb(0, 0, 0, alpha=0.6), fill = hsv(222/360, 60/100, 100/100, alpha=0))
    text3_circle<-geom_point(data=omit3[,1:2],
                             aes(x=omit3[,1], y=omit3[,2], label=sample3-1),
                             shape = 21, size = size, stroke = 1, colour = rgb(0, 0, 0, alpha=0.6), fill = hsv(222/360, 60/100, 100/100, alpha=0))
    
    # text1_init<-geom_text(aes(x=dataset1[1,1], y=dataset1[1,2], label="i"),cex =6,fontface = "bold")
    # text2_init<-geom_text(aes(x=dataset2[1,1], y=dataset2[1,2], label="i"),cex =6,fontface = "bold")
    #
    # text1_end<-geom_text(aes(x=dataset1[nrow(dataset1),1], y=dataset1[nrow(dataset1),2], label="e"),cex =6,fontface = "bold")
    # text2_end<-geom_text(aes(x=dataset2[nrow(dataset2),1], y=dataset2[nrow(dataset2),2], label="e"),cex =6,fontface = "bold")
    draw<-draw+
      # text1+
      # text2+
      text1_circle+
      text2_circle+
      text3_circle
  }
  
  draw<-draw+
    # coord_cartesian(xlim=xlimitation,ylim=ylimitation)+
    xlimitation+ylimitation+
    coord_fixed()
  draw<-draw+xlab(xlab)+ylab(ylab)
  print(draw)
  # theme(aspect.ratio=aspect_ratio)
  if(title!=""){
    ggsave(title, draw, bg = "transparent",dpi=dpi1,width=10,height=10)#width and height were previously 5
  }
  
  
}

drawDirection3d<-function(STD,OBJ,REF,xlab,ylab,name,phase=c(1:nrow(REF)),limit=0){
  dir.create("graphics")
  dir.create("graphics/pulse_direction")
  if(limit==0){
    if(max(na.omit(REF[,2]))<1.1){
      size=3.0
    }else{
      size=30
    }
  }else{
    size=limit
  }
  
  S_hv=STD[phase,]
  O_hv=OBJ[phase,]
  R_hv=REF[phase,]
  
  S_xyz=vectoring(S_hv)
  O_xyz=vectoring(O_hv)
  R_xyz=vectoring(R_hv)
  S_xyz=na.omit(S_xyz)
  O_xyz=na.omit(O_xyz)
  R_xyz=na.omit(R_xyz)
  
  
  Polar_line2(S_hv,O_hv,R_hv,8,size,xlab,ylab,"pointON_b","pathON","densityOFF","textOFF","backON",paste("graphics/pulse_direction/",name,".png", sep=""),0)
  
  if(draw_RGLfigures){
    
    if(rgl.cur()!=0){
      rgl.close()
    }
    drawSphere(1,size)
    par3d(windowRect = c(0, 0, 700, 700))
    rollx=rotationMatrix(-pi/2, 1,0,0)
    rolly=rotationMatrix(0, 0,1,0)
    rollz=rotationMatrix(pi/2, 0,0,1)
    rgl.viewpoint(userMatrix=rollx%*%rolly%*%rollz, fov=0, zoom=0.006*size, scale=par3d("scale"), interactive=TRUE)
    legend3d("topright", legend = c(head3,'Pulse','Target'), pch = 16, col = c(rgb(1,1,1),rgb(0.5,0.5,0.5),rgb(1,0.5,0.5)), cex=1, inset=c(0.02))
    
    rgl.points(0,0,0,col="black")
    rgl.points(S_xyz,col=rgb(0,0,0),alpha=0.4,size=8)
    rgl.points(O_xyz,col=rgb(0,0,0),alpha=0.4,size=8)
    rgl.points(R_xyz,col=rgb(1,0.3,0.3),alpha=0.6,size=8)
    rgl.linestrips(S_xyz,col=rgb(0,0,0),alpha=0.4,size=8)
    rgl.linestrips(O_xyz,col=rgb(0,0,0),alpha=0.4,size=8)
    rgl.linestrips(R_xyz,col=rgb(1,0.3,0.3),alpha=0.6,size=8)
    
    setwd("./graphics/pulse_direction/")
    HTML <- rglwidget( width=1000,height=600)
    filepath=paste0(name,"3D")
    fullfilepath=paste0(filepath,".html")
    htmlwidgets::saveWidget(HTML, fullfilepath,selfcontained = T)
    setwd("../../")
    rgl.close()
  }
}

drawflatcircle<-function(x,y,z,surface=3,rad=10,r=0,g=0,b=0,alpha=1){
  
  theta=seq(0,2*pi,0.5)
  axes=data.frame(theta,theta,theta)
  axes[,surface]=0
  axes[,-surface][,1]=rad*cos(theta)
  axes[,-surface][,2]=rad*sin(theta)
  
  for(i in 1:length(x)){
    circle=data.frame(axes[,1]+x[i],axes[,2]+y[i],axes[,3]+z[i])
    
    polygon3d(circle,color=rgb(r,g,b,alpha=alpha))
  }
}

drawScene<-function(pulse,position){
  dir.create("graphics")
  dir.create("graphics/scene")  
  
  olthoPlot(data.frame(position$x_bx,position$x_by),
            data.frame(position$x_mx,position$x_my),
            data.frame(pulse$x_bx,pulse$x_by),
            data.frame(pulse$x_mx,pulse$x_my),
            -pulse$peakdirectionH,
            2,"x","y","textON","backON","graphics/scene/scene2Dx_y.png")
  
  olthoPlot(data.frame(position$x_bx,position$x_bz),
            data.frame(position$x_mx,position$x_mz),
            data.frame(pulse$x_bx,pulse$x_bz),
            data.frame(pulse$x_mx,pulse$x_mz),
            pulse$peakdirectionV,
            2,"x","z","textON","backON","graphics/scene/scene2Dx_z.png")
  
  if(draw_RGLfigures){
    if(rgl.cur()!=0){
      rgl.close()
    }
    range_x=c(-0.500,5.500)
    range_y=c(-1.500,1.500)
    range_z=c(-0.500,1.500)
    plot3d(position$x_bx,position$x_by,position$x_bz,type="n",
           xlim=range_x,ylim=range_y,zlim=range_z,xlab="x",ylab="y",zlab="z")
    par3d(windowRect = c(0, 0, 700, 700))
    aspect3d(6,3,2)
    lines3d(position$x_bx,position$x_by,position$x_bz,col=rgb(0.5, 0.5, 1, alpha=0.6))
    points3d(pulse$x_bx, pulse$x_by,pulse$x_bz,col=rgb(0.5, 0.5, 1, alpha=0.6))
    lines3d(position$x_mx,position$x_my,position$x_mz,col= rgb(1, 0.5, 0.5, alpha=0.6))
    points3d(pulse$x_mx, pulse$x_my,pulse$x_mz,col= rgb(1, 0.5, 0.5, alpha=0.6))
    plsvec_x=cos(-pulse$peakdirectionH*pi/180)
    plsvec_y=sin(-pulse$peakdirectionH*pi/180)
    plsvec_z=sin(pulse$peakdirectionV*pi/180)
    gain=0.3
    for(h in 1:nrow(pulse)){
      x=c(pulse$x_bx[h],pulse$x_bx[h]+plsvec_x[h]*gain)
      y=c(pulse$x_by[h],pulse$x_by[h]+plsvec_y[h]*gain)
      z=c(pulse$x_bz[h],pulse$x_bz[h]+plsvec_z[h]*gain)
      lines3d(x,y,z,color=rgb(0.5, 0.5, 0.5),alpha=0.6,lwd=2.0)
    }
    view3d(theta=0, phi=-90, fov=60, zoom=0.8, scale=par3d("scale"), interactive=TRUE)
    legend3d("topright", legend = c(head3,'Bat','Moth'), pch = 16, col = c(rgb(1,1,1),rgb(0.5, 0.5, 1),rgb(1, 0.5, 0.5)), cex=1, inset=c(0.02))
    
    setwd("./graphics/scene/")
    HTML <- rglwidget( width=1000,height=600)
    filepath="scene3D"
    fullfilepath=paste0(filepath,".html")
    htmlwidgets::saveWidget(HTML, fullfilepath,selfcontained = T)
    setwd("../../")
    rgl.close()
  }
} 

fragmentation<-function(data){
  if(nrow(data)==1){return (list(data))}
  index=as.numeric(rownames(data))  
  j=1
  after<-numeric()
  blanklength<-numeric()
  for(i in 1:(length(index)-1)){
    if(index[i+1]-index[i]>1){
      after[j]=i
      j=j+1
    }
  }
  
  blank=c(0,after,length(index))
  for(i in 1:(length(blank)-1)){
    blanklength[i]=blank[i+1]-blank[i]
  }
  start=(blank+1)[-length(blank)]
  end=start+blanklength-1
  frag_num=length(start)
  result=list()
  for(i in 1:frag_num){
    frag=data[start[i]:end[i],]
    result=c(result,list(frag))
  }
  
  return(result)
}

continuous<-function(data){
  if(length(data)==1){return (list(data))}
  bool=0
  result<-logical()
  for(i in 1:(length(data))){
    if(data[i]&bool==0){
      result=c(result,TRUE)
    }else{
      bool=bool+1
      result=c(result,FALSE)
    }
  }
  
  return(result)
}

remove_value<-function(data){
  for(form in 1:length(data)){
    for(element in 1:length(data[[form]])){
      for(coord in 1:length(data[[form]][[element]])){
        
        if(nrow(data[[form]][[element]][[coord]])>0){
          a=data[[form]][[element]][[coord]]
          data[[form]][[element]][[coord]]=a[-1:-length(a),]
        }
        
      }
    }
  }
  return(data)
}

add_NA_head_vec<-function(data,idealnum){
  
  diff=idealnum-length(data)
  if(diff<0){
    print("error:add_NA_head_vec")
    exit()
  }
  
  for(i in rep(0,length=diff)){
    data=c(NA,data)
  }
  return(data)
}

add_NA_tail_vec<-function(data,idealnum){
  
  diff=idealnum-length(data)
  if(diff<0){
    print("error:add_NA_head_vec")
    exit()
  }
  
  for(i in rep(0,length=diff)){
    data=c(data,NA)
  }
  return(data)
}

add_NA_head2<-function(data,idealnum){
  
  diff=idealnum-nrow(data)
  if(diff<0){
    print("error:add_NA_head2")
    exit()
  }
  add=data.frame(matrix(rep(NA,ncol(data)), nrow=1))
  colnames(add)=colnames(data)
  for(i in rep(0,length=diff)){
    data=rbind(add,data)
  }
  return(data)
}


add_NA_head<-function(data,idealnum){
  
  for(form in 1:length(data)){
    for(element in 1:length(data[[form]])){
      for(coord in 1:length(data[[form]][[element]])){
        data[[form]][[element]][[coord]]=add_NA_head2(data[[form]][[element]][[coord]],idealnum)        
        
        
      }
    }
  }
  return(data)
}

add_NA_tail2<-function(data,idealnum){
  
  diff=idealnum-nrow(data)
  if(diff<0){
    print("error:add_NA_tail2")
    exit()
  }
  add=data.frame(matrix(rep(NA,ncol(data)), nrow=1))
  colnames(add)=colnames(data)
  for(i in rep(0,length=diff)){
    data=rbind(data,add)
  }
  return(data)
}

add_NA_tail<-function(data,idealnum){
  
  for(form in 1:length(data)){
    for(element in 1:length(data[[form]])){
      for(coord in 1:length(data[[form]][[element]])){
        data[[form]][[element]][[coord]]=add_NA_tail2(data[[form]][[element]][[coord]],idealnum)        
      }
    }
  }
  return(data)
}



length_formatting<-function(data){
  if(nrow(data$r2$STD$hv)>0){
    data$r2$STD$hv=rbind(NA,data$r2$STD$hv)
  }else{
    data$r2$STD$hv=rbind(data$r2$STD$hv,NA)[-1,]
  }
  
  if(nrow(data$r2$OBJ$hv)>0){
    data$r2$OBJ$hv=rbind(NA,data$r2$OBJ$hv)
  }else{
    data$r2$OBJ$hv=rbind(data$r2$OBJ$hv,NA)[-1,]
  }
  
  if(nrow(data$r2$REF$hv)>0){
    data$r2$REF$hv=rbind(NA,data$r2$REF$hv)
  }else{
    data$r2$REF$hv=rbind(data$r2$REF$hv,NA)[-1,]
  }
  
  
  if(nrow(data$s$STD$hv)>0){
    data$s$STD$hv=rbind(NA,data$s$STD$hv)
  }else{
    data$s$STD$hv=rbind(data$s$STD$hv,NA)[-1,]
  }
  
  if(nrow(data$s$OBJ$hv)>0){
    data$s$OBJ$hv=rbind(NA,data$s$OBJ$hv)
  }else{
    data$s$OBJ$hv=rbind(data$s$OBJ$hv,NA)[-1,]
  }
  
  if(nrow(data$s$REF$hv)>0){
    data$s$REF$hv=rbind(NA,data$s$REF$hv)
  }else{
    data$s$REF$hv=rbind(data$s$REF$hv,NA)[-1,]
  }
  
  return(data)
}

defineNAList<-function(){
  hv=data.frame(H=NA,V=NA,LeadAngle=NA,SideAngle=NA,Longi=NA)
  xyz=data.frame(x=NA,y=NA,z=NA)
  coordType <-list(hv=hv,xyz=xyz) #hv,xyz
  directionType <-list(STD=coordType,OBJ=coordType,REF=coordType) #STD,OBJ,REF
  nallList<-list(v=directionType,r1=directionType,r2=directionType,s=directionType)#v,r1,r2,s
  return (nallList)
}

defineNullList2<-function(){
  nullList=defineNAList()
  nullList=remove_value(nullList)
  return (nullList)
}

defineNullList<-function(){
  hv=data.frame(H=NULL,V=NULL,LeadAngle=NULL)
  xyz=data.frame(x=NULL,y=NULL,z=NULL)
  coordType <-list(hv=hv,xyz=xyz) #hv,xyz
  directionType <-list(STD=coordType,OBJ=coordType,REF=coordType) #STD,OBJ,REF
  nullList<-list(v=directionType,r1=directionType,r2=directionType,s=directionType)#v,r1,r2,s
  return (nullList)
}

rbind_List<-function(nullList,add){
  for(i in 1:length(nullList)){
    for(j in 1:length(nullList[[i]])){
      for(k in 1:length(nullList[[i]][[j]])){
        nullList[[i]][[j]][[k]]<-rbind(nullList[[i]][[j]][[k]],add[[i]][[j]][[k]])
      }
    }
  }
  return(nullList)
}

rotationALL_polar<-function(vectorSTD,vectorOBJ, vectorREF,limit1,limit2,limit3,graphics=0){
  
  
  if(nrow(vectorSTD)==0
     ||nrow(vectorOBJ)==0
     ||nrow(vectorREF)==0){
    return(defineNullList2())
    
  }
  
  
  colnames(vectorSTD)=c("x","y","z")
  colnames(vectorOBJ)=c("x","y","z")
  colnames(vectorREF)=c("x","y","z")
  
  vT_polar=direction2(vectorSTD)
  vS_polar=direction2(vectorOBJ)
  vF_polar=direction2(vectorREF)
  
  
  #rotating so as to align vectorSTD on x-axis, vectorREF on x-z plane
  alignSTD=angleForAlign(vectorSTD)
  rotatedSTD=rotationPolar(vectorSTD,alignSTD$Z,3)#rotating around z-axis
  rotatedSTD=rotationPolar(rotatedSTD,alignSTD$Y,2)#rotating around y-axis
  rotatedOBJ=rotationPolar(vectorOBJ,alignSTD$Z,3)
  rotatedOBJ=rotationPolar(rotatedOBJ,alignSTD$Y,2)
  rotatedREF=rotationPolar(vectorREF,alignSTD$Z,3)
  rotatedREF=rotationPolar(rotatedREF,alignSTD$Y,2)
  
  rT_polar=direction2(rotatedSTD)
  rS_polar=direction2(rotatedOBJ)
  rF_polar=direction2(rotatedREF)
  
  #align rotatedREF on x-z plane
  aX=atan2(rotatedREF$y,rotatedREF$z)
  rotatedSTD2=rotationPolar(rotatedSTD,aX,1)
  rotatedOBJ2=rotationPolar(rotatedOBJ,aX,1)
  rotatedREF2=rotationPolar(rotatedREF,aX,1)
  
  rT2_polar=direction2(rotatedSTD2)
  rS2_polar=direction2(rotatedOBJ2)
  rF2_polar=direction2(rotatedREF2)
  
  #option: angular scaling so as to direct rotatedREF toward 1-0-1
  aY_REF=atan2(rotatedREF2$z,rotatedREF2$x)
  
  scaledSTD=scaling_polar(rotatedSTD2,aY_REF,1)
  scaledOBJ=scaling_polar(rotatedOBJ2,aY_REF,1)
  scaledREF=scaling_polar(rotatedREF2,aY_REF,1)
  
  sT_polar=direction2(scaledSTD)
  sS_polar=direction2(scaledOBJ)
  sF_polar=direction2(scaledREF)
  
  if(graphics){
    if(rgl.cur()!=0){
      rgl.close()
    }
    drawSphere()
    plot3d(rT2_polar[[2]],xlab="x",ylab="y",zlab="z",xlim = c(-1.2,1.2),ylim = c(-1.2,1.2),zlim = c(-1.2,1.2),add = T,col="blue")
    plot3d(rS2_polar[[2]],xlab="x",ylab="y",zlab="z",xlim = c(-1.2,1.2),ylim = c(-1.2,1.2),zlim = c(-1.2,1.2),add = T,col="red")
    plot3d(rF2_polar[[2]],xlab="x",ylab="y",zlab="z",xlim = c(-1.2,1.2),ylim = c(-1.2,1.2),zlim = c(-1.2,1.2),add = T,col="green")
    bgplot3d({plot.new(); title(format(aY_REF[1]*180/pi,digits = 1))})
  }
  
  
  v=list(vT_polar,vS_polar,vF_polar)
  names(v)=c("STD","OBJ","REF")
  r1=list(rT_polar,rS_polar,rF_polar)
  names(r1)=c("STD","OBJ","REF")
  r2=list(rT2_polar,rS2_polar,rF2_polar)
  names(r2)=c("STD","OBJ","REF")
  s=list(sT_polar,sS_polar,sF_polar)
  names(s)=c("STD","OBJ","REF")
  
  return(list(v=v,r1=r1,r2=r2,s=s))
}

LRT<-function(testmodel, controlmodels){
  result=data.frame(rep(0,1))[,-1]
  for(i in 1:length(controlmodels)){
    anova<-anova(testmodel, controlmodels[[i]], test = "Chisq")
    res=(anova$"Pr(>Chisq)")[2]
    result=cbind(result, as.data.frame(res))
    colnames(result)[ncol(result)]=c(paste0("model",i-1))
  }
  rownames(result)=c("p in ANOVA(LRT)")
  return(result)
}

model_selection_addNLS<-function(x,y,ID,Ex,poly){
  
  d<-data.frame(x,y,ID,Ex)
  d_omit<-as.data.frame(d[is.finite(rowSums(d)),])
  colnames(d_omit)=c("x","y","ID","Ex")
  x=d_omit$x
  y=d_omit$y
  ID=d_omit$ID
  Ex=d_omit$Ex
  
  modellist1=generate.MM_with_intercept(x,y,ID,Ex,poly)
  modellist2=generate.MM_NLS(x,y,ID,Ex,poly)
  allmodels=c(modellist1,modellist2[[1]])
  index_aic=AICselection(allmodels)
  bestmodel=c(allmodels)[[index_aic]]
  
  num_model=length(allmodels)
  
  nullmodel0=generate.nullmodel.MM(x,y,ID,Ex,0)
  nullmodel1=generate.nullmodel.MM(x,y,ID,Ex,1)
  coef_a=numeric()
  coef_b=numeric()
  aics=numeric()
  lrt0s=numeric()
  lrt1s=numeric()
  modtypes=character()
  
  for(i in 1:num_model){
    a=fixed.effects(allmodels[[i]])[2]
    b=fixed.effects(allmodels[[i]])[1]
    aic=AIC(allmodels[[i]])
    lrt0=LRT(allmodels[[i]],list(nullmodel0,nullmodel1))[,1]
    lrt1=LRT(allmodels[[i]],list(nullmodel0,nullmodel1))[,2]
    if(i!=num_model){
      mod.type=summary(allmodels[[i]])[15]
      mod=substr(mod.type,regexpr("y ~", mod.type),regexpr("(1 | Ex)", mod.type)-4)
    }else{
      mod="y ~ a/x + b, (a>=0)" 
      tmp=a
      a=b
      b=tmp
    }
    
    
    coef_a=c(coef_a,a)  
    coef_b=c(coef_b,b)  
    aics=c(aics,aic)  
    lrt0s=c(lrt0s,lrt0)  
    lrt1s=c(lrt1s,lrt1)  
    modtypes=c(modtypes,mod)
  }
  
  result=cbind("model"=modtypes,"a"=coef_a,"b"=coef_b,"AIC"=aics,"p (y=0)"=lrt0s,"p (y=b)"=lrt1s)
  rownames(result)=c(1:nrow(result))
  
  aiclabel=rep("",nrow(result))
  aiclabel[index_aic]="<-lowestAIC"
  result=cbind(result,aiclabel)
  back<-list(bestmodel,result,modellist2[[2]])
  return(back)
}

model_selection<-function(x,y,ID,Ex,poly,tgtpoly=-1){
  
  if(is.null(x)){
    x=c(1:length(y))
  }
  
  d<-data.frame(x,y,ID,Ex)
  d_omit<-as.data.frame(d[is.finite(rowSums(d)),])
  colnames(d_omit)=c("x","y","ID","Ex")
  x=d_omit$x
  y=d_omit$y
  ID=d_omit$ID
  Ex=d_omit$Ex
  
  modellist1=generate.MM_with_intercept(x,y,ID,Ex,poly)
  index_aic=AICselection(c(modellist1))
  bestmodel=c(modellist1)[[index_aic]]
  
  num_model=length(modellist1)
  
  
  if(tgtpoly!=-1){
    bestmodel=modellist1[[tgtpoly]]
  }
  
  nullmodel0=generate.nullmodel.MM(x,y,Ex,ID,0)
  nullmodel1=generate.nullmodel.MM(x,y,Ex,ID,1)
  coef_a=numeric()
  coef_b=numeric()
  aics=numeric()
  lrt0s=numeric()
  lrt1s=numeric()
  modtypes=character()
  
  for(i in 1:num_model){
    a=fixed.effects(modellist1[[i]])[2]
    b=fixed.effects(modellist1[[i]])[1]
    aic=AIC(modellist1[[i]])
    lrt0=LRT(modellist1[[i]],list(nullmodel0,nullmodel1))[,1]
    lrt1=LRT(modellist1[[i]],list(nullmodel0,nullmodel1))[,2]
    mod.type=summary(modellist1[[i]])[15]
    mod=substr(mod.type,regexpr("y ~", mod.type),regexpr("(1 | Ex)", mod.type)-4)
    
    coef_a=c(coef_a,a)  
    coef_b=c(coef_b,b)  
    aics=c(aics,aic)  
    lrt0s=c(lrt0s,lrt0)  
    lrt1s=c(lrt1s,lrt1)  
    modtypes=c(modtypes,mod)
  }
  
  result=cbind("model"=modtypes,"a"=coef_a,"b"=coef_b,"AIC"=aics,"p (y=0)"=lrt0s,"p (y=b)"=lrt1s)
  rownames(result)=c(1:nrow(result))
  aiclabel=rep("",nrow(result))
  aiclabel[index_aic]="<-lowestAIC"
  result=cbind(result,aiclabel)
  back<-list(bestmodel,result)
  
  return(back)
}

ggplot_drawMM<-function(x,y,ID,Ex,poly,xlim,ylim,color,xlab,ylab,title, xyline=0,transparent=0, drawing=1,tgtpoly=-1){
  d<-data.frame(x,y,ID,Ex)
  d_omit<-as.data.frame(d[is.finite(rowSums(d)),])
  colnames(d_omit)=c("x","y","ID","Ex")
  x=d_omit$x
  y=d_omit$y
  ID=d_omit$ID
  Ex=d_omit$Ex
  
  modeling_result=model_selection(x,y,ID,Ex,poly,tgtpoly)
  best.model=modeling_result[[1]]
  print(best.model)
  print("modellist:")
  print(modeling_result[[2]])
  write.csv(x = modeling_result[[2]], file = paste0(title,".csv"))
  
  predx=(seq(min(xlim),max(xlim),length=1000))
  newdata_mix=data.frame("x"=predx, "ID"=(rep(min(ID):max(ID),length=1000)),"Ex"=(rep(min(Ex):max(Ex),length=1000)))
  predy=predictInterval(merMod = best.model, newdata = newdata_mix,
                        level = 0.95, n.sims = 1000,
                        stat = "median", type="linear.prediction",
                        include.resid.var=F,
                        which = "fixed") 
  bestreg=data.frame(predx,predy[,1])
  confint_up=data.frame(predx,predy[,2])
  confint_down=data.frame(predx,predy[,3])
  
  if(drawing==1){
    draw<-ggplot() 
    
    colnames(confint_up)=colnames(confint_down)=c("x","y")
    confint<-rbind(confint_up,confint_down[rev(1:nrow(confint_down)),])
    inarea=which(confint[,1]>=xlim[1]&confint[,1]<=xlim[2]&confint[,2]>=ylim[1]&confint[,2]<=ylim[2])
    confint=confint[inarea,]
    draw<-draw+
      geom_polygon(data=confint,aes(x=confint[,1],y=confint[,2]),fill = rgb(0,0,0, alpha=0.1))
    
    
    draw<-draw+
      geom_point(data=d_omit, aes(x=x, y=y),
                 shape = 21, size = 3.5, stroke = 0.3, colour = rgb(0, 0, 0), fill = rgb(0,0,0, alpha=0.5))
    # draw<-draw+
    #   geom_hline(yintercept = 0,size=0.5,linetype="solid", colour = rgb(0,0,0),alpha=0.7)
    # 
    # hline=data.frame(xlim,c(0,0))
    # draw<-draw+geom_line(data=hline,aes(x=hline[,1],y=hline[,2]),size=0.5,linetype="solid", colour = rgb(0,0,0),alpha=0.7)
    # 
    if(xyline!=0){
      a=xyline-1
      xy=data.frame(xlim,a*xlim)
      
      # draw<-draw+
      #   geom_abline(slope = xyline-1, intercept = 0,linetype="solid",size=1.5, colour = hsv(222/360, 60/100, 100/100, alpha=0.9))
      draw<-draw+geom_line(data=xy,aes(x=xy[,1],y=xy[,2]),size=3, colour = hsv(222/360, 60/100, 100/100, alpha=1))
      
    }
    
    draw<-draw+geom_line(data=bestreg,aes(x=bestreg[,1],y=bestreg[,2]),colour = color, size=3)
    
    
    if(transparent){
      draw<-draw+theme(
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        axis.text= element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank()
      )
    }
    
    
    draw<-draw+xlim(xlim)+ylim(ylim)+xlab(xlab)+ylab(ylab)+ggtitle(title)
    print(draw)
    ggsave(paste0(title,".png"), draw, bg = "transparent",dpi=dpi2,width=11.3,height=10)
  }
  return(best.model)
  
  
}

regression<-function(mod,xlim,ylim){
  predx=(seq(min(xlim),max(xlim),length=1000))
  frag="lmerModLmerTest"==summary(mod[[1]])$"objClass"[1]
  if(frag){
    rf=summary(mod[[1]])$ngrps
    newdata_mix=data.frame("x"=predx, "ID"=(rep(1:max(rf[2]),length=1000)),"Ex"=(rep(1:max(rf[1]),length=1000)))
    predy=predictInterval(merMod = mod[[1]], newdata = newdata_mix,
                          level = 0.95, n.sims = 1000,
                          stat = "median", type="linear.prediction",
                          include.resid.var=F,
                          which = "fixed") 
    
  }else{
    predy=predictCI(predx,mod[[1]],mod[[3]],1000)
  }
  reg=data.frame(x=predx,y=predy[,1])
  confint_up=data.frame(predx,predy[,2])
  confint_down=data.frame(predx,predy[,3])
  
  colnames(confint_up)=colnames(confint_down)=c("x","y")
  confint<-rbind(confint_up,confint_down[rev(1:nrow(confint_down)),])
  inarea=which(confint[,1]>=xlim[1]&confint[,1]<=xlim[2]&confint[,2]>=ylim[1]&confint[,2]<=ylim[2])
  confint=confint[inarea,]
  
  return(list(reg,confint))
}


ggplot_test<-function(draw,x,y,xlim,ylim,xlab,ylab,title=NULL,mod=NULL,color=NA,xyline=0,back="backON",transparent=NULL){
  
  ci=list()
  regline=list()
  point=list()
  refline=list()
  c1=c2=rgb(0.3,0.3,0.3,alpha = 0.5)
  
  
  if(!is.null(mod)){
    for(i in 1:length(mod)){
      regdata=regression(mod[[i]],xlim,ylim)
      confint=geom_polygon(data=regdata[[2]],aes(x,y),fill = rgb(0,0,0, alpha=0.1))
      reg=geom_line(data=regdata[[1]],aes(x,y),colour = color, size=3)
      ci=c(ci,list(confint))
      regline=c(regline,list(reg))
    }
  }
  
  if(!is.null(x)&!is.null(y)){
    d<-data.frame(x,y)
    d_omit<-as.data.frame(d[is.finite(rowSums(d)),])
    colnames(d_omit)=c("x","y")
    
    p=geom_point(data=d_omit, aes(x, y),
                 shape = 21, size = 3.5, stroke = 0.0, colour = c1, fill = c2)
    point=c(point,p)
  }
  
  if(xyline!=0){
    a=xyline-1
    xy=data.frame(xlim,a*xlim)
    ref=geom_line(data=xy,aes(x=xy[,1],y=xy[,2]), colour = hsv(222/360, 60/100, 100/100, alpha=1),size=3)
    refline=c(refline,list(ref))
  }
  
  
  
  if(back=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA),
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text= element_blank()
    )
  }
  
  draw<-draw+ci+point+refline+regline
  if(!is.null(xlim))
    draw<-draw+xlim(xlim)
  if(!is.null(ylim))
    draw<-draw+ylim(ylim)
  
  draw<-draw+xlab(xlab)+ylab(ylab)
  if(!is.null(title)){
    draw<-draw+ggtitle(title)
  }
  return(draw)
}

multiplot<-function(mod,xlim,ylim,xlab,ylab,title="",color=NULL,back="backON"){
  
  draw<-ggplot() 
  regline=list()
  
  if(!is.null(mod)){
    for(i in 1:length(mod)){
      regdata=regression(mod[[i]],c(xlim[1]-15,xlim[2]+15),c(ylim[1]-15,ylim[2]+15))
      reg=geom_line(data=regdata[[1]],aes(x,y),colour = color[[i]], size=2)
      regline=c(regline,list(reg))
    }
  }
  
  
  if(back=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA),
      plot.title = element_blank(),
      # axis.title = element_blank(),
      axis.text= element_blank()
    )
  }
  
  draw<-draw+regline
  # if(!is.null(xlim))
  #   draw<-draw+xlim(xlim)
  # if(!is.null(ylim))
  #   draw<-draw+ylim(ylim)
  draw<-draw+coord_cartesian(xlim=xlim,ylim=ylim)
  
  draw<-draw+xlab(xlab)+ylab(ylab)
  draw<-draw+ggtitle(title)
  return(draw)
}


ggplotter<-function(x,y,xlim,ylim,xlab,ylab,title,back="backON", p=0,save=0){
  x=as.numeric(x)
  y=as.numeric(y)
  data=data.frame(x,y)
  colnames(data)=c("x","y")
  
  draw<-ggplot() 
  draw<-draw+geom_line(data=data,aes(x,y), size=1.5)
  
  if(back=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA),
      axis.text= element_blank(),
      plot.title = element_blank(),
      axis.title = element_blank()
    )
  }
  draw<-draw+coord_cartesian(xlim=xlim,ylim=ylim)+scale_x_continuous(breaks=seq(xlim[1],xlim[2],diff(x_lim)/4))+
    xlab(xlab)+ylab(ylab)+ggtitle(title)
  
  if(save){
    ggsave(paste0(title,".png"), draw, bg = "transparent",dpi=dpi2,width=10,height=5)
  }
}

ggmultiplotter<-function(x,ylist,xlim,ylim,xlab,ylab,title,back="backON", p=0,save=0){
  x=as.numeric(x)
  
  line_list=list()
  for(i in 1:length(ylist)){
    
    y=as.numeric(ylist[[i]])
    data=data.frame(x,y)
    colnames(data)=c("x","y")
    
    line<-geom_line(data=data,aes(x,y), colour = hsv((i-1)/length(ylist)*1.2, 0.7, 0.7), size=1.5)
    line_list<-c(line_list,list(line))
  }
  
  draw<-ggplot() 
  for(i in 1:length(line_list)){
    draw<-draw+line_list[[i]]
  }
  
  if(back=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA),
      axis.text= element_blank(),
      plot.title = element_blank(),
      axis.title = element_blank()
    )
  }
  draw<-draw+coord_cartesian(xlim=xlim,ylim=ylim)+scale_x_continuous(breaks=seq(xlim[1],xlim[2],diff(x_lim)/4))+
    xlab(xlab)+ylab(ylab)+ggtitle(title)
  
  if(save){
    ggsave(paste0(title,".png"), draw, bg = "transparent",dpi=dpi2,width=10,height=5)
  }
}

coefficient_fix<-function(best.model){
  coeff_fix=as.data.frame(t(as.data.frame(fixef(best.model))))
  coeff=numeric()
  
  if(is.null(coeff_fix$'(Intercept)')){
    value=0
  }else{
    value=coeff_fix$'(Intercept)'
  }
  coeff=c(coeff,value)
  
  if(is.null(coeff_fix$'x')){
    value=0
  }else{
    value=coeff_fix$'x'
  }
  coeff=c(coeff,value)
  
  if(is.null(coeff_fix$'I(x^2)')){
    value=0
  }else{
    value=coeff_fix$'I(x^2)'
  }
  coeff=c(coeff,value)
  
  if(is.null(coeff_fix$'I(x^3)')){
    value=0
  }else{
    value=coeff_fix$'I(x^3)'
  }
  coeff=c(coeff,value)
  
  return(coeff)
}

generate.nullmodel.MM<-function(x,y,ID,Ex,intercept){
  d=data.frame(x,y,ID,Ex)
  d_omit<-as.data.frame(d[is.finite(rowSums(d)),])
  switch (as.character(intercept),
          "0" =  lmer(y ~ 0 + (1|ID)+ (1|Ex), data=d_omit)->null.model,
          "1" =  lmer(y ~ 1 + (1|ID)+ (1|Ex), data=d_omit)->null.model,
  )
  return(null.model)
}

AICselection<-function(modellist){
  AIC_list=numeric()
  print("AIC List")
  for(i in 1:length(modellist)){
    AIC_list=c(AIC_list,AIC(modellist[[i]]))
    
    print(paste(round(AIC(modellist[[i]]),3),": ",
                as.character(formula(modellist[[i]])[2]), ", ",
                as.character(formula(modellist[[i]])[3])))
    
  }
  
  index=which(AIC_list==min(AIC_list))
  lowestAIC.model=modellist[[index]]
  return(index)
}

generate.MM_with_intercept<-function(x,y,ID,Ex, poly){
  
  d=data.frame(x,y,ID,Ex)
  d_omit<-as.data.frame(d[is.finite(rowSums(d)),])
  options(na.action="na.fail")
  switch (as.character(poly),
          "0" =  lmer(y ~ 1 + (1|Ex) + (1|ID), data=d_omit)->model,
          "1" =  lmer(y ~ x + (1|Ex) + (1|ID), data=d_omit)->model,
          "2" =  lmer(y ~ x+I(x^2) + (1|Ex) + (1|ID), data=d_omit)->model,
          "3" =  lmer(y ~ x+I(x^2)+I(x^3) + (1|Ex) + (1|ID), data=d_omit)->model,
          stop("Acceptable poly value is Only 0, 1, 2 or 3")
  )
  
  all.model<-dredge(model,rank="AIC")
  sort.model<-all.model[order(as.numeric(rownames(all.model))),]
  print(sort.model)
  
  modellist<-get.models(sort.model,subset = 1:nrow(all.model))
  options(na.action="na.omit")
  return(modellist)
}

generate.MM_without_intercept<-function(x,y,ID, poly){
  d=data.frame(x,y,ID)
  d_omit<-as.data.frame(d[is.finite(rowSums(d)),])
  options(na.action="na.fail")
  switch (as.character(poly),
          "0" =  lmer(y ~ 0 + (1|ID), data=d_omit)->model,
          "1" =  lmer(y ~ 0 + x + (1|ID), data=d_omit)->model,
          "2" =  lmer(y ~ 0 + x+I(x^2) + (1|ID), data=d_omit)->model,
          "3" =  lmer(y ~ 0 + x+I(x^2)+I(x^3) + (1|ID), data=d_omit)->model,
          stop("Acceptable poly value is Only 0, 1, 2 or 3")
  )
  
  all.model<-dredge(model,rank="AIC")
  sort.model<-all.model[order(as.numeric(rownames(all.model))),]
  print(sort.model)
  
  modellist<-get.models(all.model,subset = 1:nrow(all.model))
  options(na.action="na.omit")
  return(modellist)
}

generate.MM_NLS<-function(x,y,ID,Ex,poly){
  d=data.frame(x,y,ID,Ex)
  d_omit<-as.data.frame(d[is.finite(rowSums(d)),])
  
  f1 <- deriv(~a/x+b, namevec = c("a", "b"), function.arg = c("x","a","b"))
  tart.vec <- list(a = 0,b = 0)
  model_1 <- nlmer(y~f1(x,a,b)~(a+b|ID)+(a+b|Ex), start = c(a = 0, b = 0), data=d_omit, na.action = "na.omit",
                   control=nlmerControl(optimizer = "nloptwrap",
                                        optCtrl = list(
                                          maxit=10000,
                                          iter.max=10000,  
                                          eval.max=10000,
                                          lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0),
                                          upper = c(Inf,Inf,Inf,Inf,Inf,Inf,Inf))))
  #upper/lower index
  #c(b,randomeffect,randomeffect,randomeffect,randomeffect,randomeffect,a)
  
  # model_1 <- update(model_1, control=nlmerControl(optimizer="nloptwrap"), # Default optimizer
  #                   nAGQ=0L)
  
  print(fixed.effects(model_1))
  return(list(model_1,f1))
}

predictCI<-function(predx,model,f,n){
  Mean <- fixef(model)
  Var <- vcov(model)
  
  ###### Sampling from a probability distribution of parameters ######
  set.seed(1)
  samp <- rmvnorm(n, Mean, as.matrix(Var)) 
  
  ###### Repeating the model building and predictions ######
  samp_pred_array <- NULL
  for (i in 1:n) {
    a=numeric()
    for(j in 1:ncol(samp)){
      a <-c(a, samp[i,j]) 
    }
    
    if(length(a)<2){a=c(a,0)}
    samp_pred <- mapply(f,predx, a[1],a[2]) 
    samp_pred_array <- rbind(samp_pred_array, samp_pred)
  }
  
  if(length(Mean)<2){Mean=c(Mean,0)}
  regression <- mapply(f,predx, Mean[1],Mean[2]) 
  
  CI <- apply(samp_pred_array, MARGIN=2, function(x) quantile(x, probs=c(0.025, 0.975)))
  return(data.frame(regression,t(CI)))
}

multi_distribution_comp<-function(x,y,ID,Ex){
  comp0_bool=data.frame(rep(0,length(y)))[,-1]
  comp1_bool=data.frame(rep(0,length(y)))[,-1]
  comp0_value=data.frame(rep(0,length(y)))[,-1]
  comp1_value=data.frame(rep(0,length(y)))[,-1]
  for(s in 1:length(y)){
    anova0_bool=numeric()
    anova1_bool=numeric()
    anova0_value=numeric()
    anova1_value=numeric()
    for(t in 1:length(y)){
      print(paste(paste("s:",s),paste("t:",t)))
      
      if(t==s){
        anova0_bool=c(anova0_bool,NA)
        anova1_bool=c(anova1_bool,NA)
        anova0_value=c(anova0_value,NA)
        anova1_value=c(anova1_value,NA)
      }else{
        Y=y[[t]]-y[[s]]
        anovares=model_selection(x,Y,ID,Ex,1)[[2]]
        
        #test
        siglevel=0.05/(length(y)*(length(y)-1)/2)
        print(paste("alpha=",siglevel))
        anova0_bool=c(anova0_bool,anovares[,1]<siglevel)
        anova1_bool=c(anova1_bool,anovares[,2]<siglevel)
        anova0_value=c(anova0_value,anovares[,1])
        anova1_value=c(anova1_value,anovares[,2])
      }
      
    }
    comp0_bool=cbind(comp0_bool, as.data.frame(anova0_bool))
    comp1_bool=cbind(comp1_bool, as.data.frame(anova1_bool))
    comp0_value=cbind(comp0_value, as.data.frame(anova0_value))
    comp1_value=cbind(comp1_value, as.data.frame(anova1_value))
    colnames(comp0_bool)[length(comp0_bool)]=c(as.character(s))
    colnames(comp1_bool)[length(comp1_bool)]=c(as.character(s))
    colnames(comp0_value)[length(comp0_value)]=c(as.character(s))
    colnames(comp1_value)[length(comp1_value)]=c(as.character(s))
    
  }
  
  print("multicomp: anova0")
  print(comp0_value)
  print(comp0_bool)
  print("multicomp: anova1")
  print(comp1_value)
  print(comp1_bool)
  
  return (list(simtest_bool=comp0_bool,
               simtest_value=comp0_value,
               cortest_bool=comp1_bool,
               cortest_value=comp1_value))
}

multi_distribution_comp2<-function(data,ID,Ex){
  comp0_bool=data.frame(rep(0,length(data)))[,-1]
  comp0_value=data.frame(rep(0,length(data)))[,-1]
  for(s in 1:length(data)){
    anova0_bool=numeric()
    anova0_value=numeric()
    for(t in 1:length(data)){
      print(paste(paste("s:",s),paste("t:",t)))
      
      if(t==s){
        anova0_bool=c(anova0_bool,NA)
        anova0_value=c(anova0_value,NA)
      }else{
        diff=data[[t]][,2]-data[[s]][,2]
        anovares=model_selection(AllSessions$av_mag,diff,ID,Ex,1)[[2]]
        #test
        siglevel=0.05/(length(data)*(length(data)-1)/2)
        print(paste("alpha=",siglevel))
        anova0_bool=c(anova0_bool,anovares[,1]<siglevel)
        anova0_value=c(anova0_value,anovares[,1])
      }
      
    }
    comp0_bool=cbind(comp0_bool, as.data.frame(anova0_bool))
    comp0_value=cbind(comp0_value, as.data.frame(anova0_value))
    colnames(comp0_bool)[length(comp0_bool)]=c(as.character(s))
    colnames(comp0_value)[length(comp0_value)]=c(as.character(s))
  }
  
  print("multicomp: anova0")
  print(comp0_value)
  print(comp0_bool)
  
  
  return (list(simtest_bool=comp0_bool,
               simtest_value=comp0_value))
}

multi_plot_modeling_comp<-function(x,y,ID,Ex,title,xlimit=NA,ylimit=NA){
  res=multi_distribution_comp(x,y,ID,Ex)
  multi_plot_modeling(x,y,ID,Ex,title,xlimit,ylimit,"backON")
  return(res)  
}

multi_plot_modeling_comp2<-function(data,ID,Ex,title,xlimit=NA,ylimit=NA){
  res=multi_distribution_comp2(data,ID,Ex)
  multi_plot_modeling2(data,ID,Ex,title,xlimit,ylimit,"backON")
  return(res)  
}

multi_plot_modeling<-function(x,y,ID,Ex,title,xlimit=NA,ylimit=NA,backOFF="backOFF"){
  
  
  yvalue=numeric()
  model_list=list()
  
  sim_anova0=numeric()
  sim_anova1=numeric()
  for(i in 1:length(y)){
    print(as.character(i))
    
    modeling_result=model_selection(x,y[[i]],ID,Ex,1)
    best.model=modeling_result[[1]]
    print(best.model)
    print(fixed.effects(best.model))
    print("LRT:")
    print(modeling_result[[2]])
    anovares=modeling_result[[2]]
    sim_anova0=c(sim_anova0,anovares[,1]<0.05)  
    sim_anova1=c(sim_anova1,anovares[,2]<0.05)      
    
    
    
    predx=(seq(min(x),max(x),length=1000))
    frag="lmerModLmerTest"==summary(best.model)$"objClass"[1]
    if(frag){
      newdata_mix=data.frame("x"=predx, "ID"=(rep(min(ID):max(ID),length=1000)),"Ex"=(rep(min(Ex):max(Ex),length=1000)))
      predy=predictInterval(merMod = best.model, newdata = newdata_mix,
                            level = 0.95, n.sims = 1000,
                            stat = "median", type="linear.prediction",
                            include.resid.var=F,
                            which = "fixed") 
      
    }else{
      predy=predictCI(predx,best.model,modeling_result[[3]],1000)
    }
    
    
    
    
    pred=data.frame(predx,predy[,1])
    model_list=c(model_list, list(pred))
    yvalue=c(yvalue,pred[,2])
  }
  
  print("eachReg: anova0")
  print(sim_anova0)
  print("eachReg: anova1")
  print(sim_anova1)
  
  
  
  xvalue=c(0,800)
  
  if(all(!is.na(xlimit))){
    xvalue=xlimit 
  }
  if(all(!is.na(ylimit))){
    yvalue=ylimit 
  }
  
  col=hsv(0,0,0)
  col1=hsv(0, 60/100, 0/4)
  col2=hsv(0, 60/100, 1/4)
  col3=hsv(0, 60/100, 2/4)
  col4=hsv(0, 60/100, 3/4)
  col5=hsv(0.5, 60/100, 0/4)
  col6=hsv(0.5, 60/100, 1/4)
  col7=hsv(0.5, 60/100, 2/4)
  col8=hsv(0.5, 60/100, 3/4)
  
  colset=c(col,col,col,col,
           col,col,col,col5)
  # colset[i]=rgb(1, 0.35, 0.4)
  
  
  p<-ggplot()
  p<-p+geom_line(data=model_list[[1]],aes(x=model_list[[1]][,1],y=model_list[[1]][,2]),colour = col1,size=1)
  p<-p+geom_line(data=model_list[[2]],aes(x=model_list[[2]][,1],y=model_list[[2]][,2]),colour = col2,size=1)
  p<-p+geom_line(data=model_list[[3]],aes(x=model_list[[3]][,1],y=model_list[[3]][,2]),colour = col3,size=1)
  p<-p+geom_line(data=model_list[[4]],aes(x=model_list[[4]][,1],y=model_list[[4]][,2]),colour = col4,size=1)
  p<-p+geom_line(data=model_list[[5]],aes(x=model_list[[5]][,1],y=model_list[[5]][,2]),colour = col5,size=1)
  p<-p+geom_line(data=model_list[[6]],aes(x=model_list[[6]][,1],y=model_list[[6]][,2]),colour = col6,size=1)
  p<-p+geom_line(data=model_list[[7]],aes(x=model_list[[7]][,1],y=model_list[[7]][,2]),colour = col7,size=1)
  p<-p+geom_line(data=model_list[[8]],aes(x=model_list[[8]][,1],y=model_list[[8]][,2]),colour = col8,size=1)
  
  num=c(1:8)
  labx=numeric()
  laby=numeric()
  for(i in 1:8){
    length=nrow(model_list[[i]])
    labx=c(labx,model_list[[i]][length,1])
    laby=c(laby,model_list[[i]][length,2])
  }
  
  label_posit<-data.frame(labx+3,laby)
  # text<-geom_text(data=label_posit,
  #                  aes(x=label_posit[,1], y=label_posit[,2], label=num),cex=7,fontface = "bold")
  # 
  # p<-p+text
  if(backOFF=="backOFF"){
    p<-p+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA),
      axis.text= element_blank(),
      plot.title = element_blank()
      # axis.ticks= element_blank()#
    )
  }
  # yrange=abs(max(yvalue)-min(yvalue))
  # p<-p+geom_hline(yintercept = 0,size=0.3, colour = rgb(0,0,0),alpha=0.5)+
  #   geom_hline(yintercept = -1,size=0.3, colour = rgb(0,0,0),alpha=0.5)+
  #   geom_vline(xintercept = 0,size=0.5, colour = rgb(0,0,0))
  p<-p + scale_y_continuous(breaks=c(0,-1))
  p<-p+xlim(xvalue[1],xvalue[2])+ylim(min(yvalue),max(yvalue))
  p<-p + xlab("Angular Velocity")+ylab("Sighting Error")
  print(p)
  
  if(title!=""){
    ggsave(title, p, bg = "transparent",dpi=dpi2,width=11.5,height=10)
  }
  
  
}

gather_pred_res<-function(data,ID,Ex,predx,src,lab){
  
  for(i in 1:length(data)){
    print(as.character(i))
    mod=model_selection(data[[i]][,1],data[[i]][,2],ID,Ex,1)
    best.model=mod[[1]]
    
    
    
    
    frag="lmerModLmerTest"==summary(best.model)$"objClass"[1]
    if(frag){
      newdata_mix=data.frame("x"=predx, "ID"=(rep(min(ID):max(ID),length=length(predx))),"Ex"=(rep(min(Ex):max(Ex),length=length(predx))))
      predy=predictInterval(merMod = best.model, newdata = newdata_mix,
                            level = 0.95, n.sims = 1000,
                            stat = "median", type="linear.prediction",
                            include.resid.var=F,
                            which = "fixed") 
      
    }else{
      predy=predictCI(predx,best.model,modeling_result[[3]],length(predx))
    }
    
    
    
    
    pred=data.frame(predx,predy[,1])
    name=paste0(lab,as.character(i))
    add=as.data.frame(pred[,2])
    colnames(add)=name
    src=cbind(src,add)
  }
  
  return(src)
  
}

multi_plot_modeling2<-function(data,ID,Ex,xlimit=NA,ylimit=NA,xlab,ylab,
                               size=2,colset=NULL,backOFF="backOFF",
                               title=NULL,F_dir=".",T_dir="."){
  
  
  yvalue=numeric()
  model_list=list()
  
  sim_anova0=numeric()
  sim_anova1=numeric()
  for(i in 1:length(data)){
    print(as.character(i))
    mod=model_selection(data[[i]][,1],data[[i]][,2],ID,Ex,1)
    model_list=c(model_list, list(mod))
    
    if(!is.null(T_dir)){
      write.csv(mod[[2]],paste0(T_dir,"/",title,"_",i,".csv"))
    }
  }
  
  if(is.null(colset)){
    vlist=c(1/4,1.5/4,2.5/4,3.5/4)
    for(i in 1:length(data)){
      h=1-1/((i-1)%/%4+1)
      v=vlist[(i-1)%%4+1]
      colset=c(colset,list(hsv(h, 60/100, v)))
    }
  }
  fig=multiplot(model_list,xlimit,ylimit,xlab,ylab,title,colset,backOFF)
  
  if(!is.null(F_dir)){
    ggsave(paste0(F_dir,"/",title,".png"), fig, bg = "transparent",dpi=dpi2,width=10,height=10)
  }
  return(fig)
  
}

linear_pred<-function(step1_set,step2_set,deltatime_set,leadtime_set){
  lengthcheck=all(sapply(list(nrow(step1_set),
                              nrow(step2_set),
                              length(deltatime_set),
                              length(leadtime_set)), 
                         function(x) x == length(leadtime_set)))
  if(!lengthcheck){
    print("length not matched")
    break()
  }
  
  d=na.omit(data.frame(step1_set=step1_set,step2_set=step2_set,deltatime_set,leadtime_set))
  step1_set=data.frame(d$step1_set.x,d$step1_set.y,d$step1_set.z)
  step2_set=data.frame(d$step2_set.x,d$step2_set.y,d$step2_set.z)
  deltatime_set=d$deltatime_set
  leadtime_set=d$leadtime_set
  
  result=data.frame(matrix(rep(NA,3), nrow=1))[numeric(0), ]
  colnames(result)=c("x","y","z")
  
  if(nrow(step1_set)==0){
    return(result)
  }
  
  
  
  for(i in 1:nrow(step1_set)){
    
    step1=as.numeric(step1_set[i,])
    step2=as.numeric(step2_set[i,])
    deltatime=deltatime_set[i]
    leadtime=leadtime_set[i]
    
    x_delta=step2-step1
    v=x_delta/deltatime
    pred_d=step2+v*leadtime
    result=rbind(result,t(as.data.frame(pred_d)))
  }
  colnames(result)=c("x","y","z")
  rownames(result)=c(1:nrow(step1_set))
  
  mag=sqrt(result$x^2+result$y^2+result$z^2)
  return(result) 
}

absTGTV_pred_core<-function(step1,step2,TGTV,deltatime,lead_gain,lead_offset=0){
  step1=as.numeric(step1)
  step2=as.numeric(step2)
  TGTV=as.numeric(TGTV)
  deltatime=as.numeric(deltatime)
  lead_gain=as.numeric(lead_gain)
  lead_offset=as.numeric(lead_offset)
  
  n=cross(step2-TGTV,step2)
  av=angle_2vec(step1,step2)*pi/180/deltatime
  leadang=av*lead_gain+lead_offset
  step3=rotationMatrix(leadang,n[1],n[2],n[3])[1:3,1:3]%*%step2
  step3_mag=sqrt(step3[1,1]^2+step3[2,1]^2+step3[3,1]^2)
  step3_dir=step3/step3_mag
  return(step3_dir)
}


angular_pred_core<-function(step1,step2,deltatime,lead_gain,lead_offset=0){
  step1=as.numeric(step1)
  step2=as.numeric(step2)
  deltatime=as.numeric(deltatime)
  lead_gain=as.numeric(lead_gain)
  lead_offset=as.numeric(lead_offset)
  
  n=cross(step1,step2)
  av=angle_2vec(step1,step2)*pi/180/deltatime
  leadang=av*lead_gain+lead_offset
  step3=rotationMatrix(leadang,n[1],n[2],n[3])[1:3,1:3]%*%step2
  step3_mag=sqrt(step3[1,1]^2+step3[2,1]^2+step3[3,1]^2)
  step3_dir=step3/step3_mag
  return(step3_dir)
}


angular_pred<-function(step1_set,step2_set,deltatime_set,leadtime_set){
  
  lengthcheck=all(sapply(list(nrow(step1_set),
                              nrow(step2_set),
                              length(deltatime_set),
                              length(leadtime_set)), 
                         function(x) x == length(leadtime_set)))
  if(!lengthcheck){
    print("length not matched")
    break()
  }
  d=na.omit(data.frame(step1_set=step1_set,step2_set=step2_set,deltatime_set,leadtime_set))
  step1_set=data.frame(x=d$step1_set.x,y=d$step1_set.y,z=d$step1_set.z)
  step2_set=data.frame(x=d$step2_set.x,y=d$step2_set.y,z=d$step2_set.z)
  deltatime_set=d$deltatime_set
  leadtime_set=d$leadtime_set
  
  result=data.frame(matrix(rep(NA,3), nrow=1))[numeric(0), ]
  colnames(result)=c("x","y","z")
  
  if(nrow(step1_set)==0){
    return(result)
  }
  
  for(i in 1:nrow(step1_set)){
    
    step3_dir=angular_pred_core(step1_set[i,],step2_set[i,],deltatime_set[i],leadtime_set[i])
    # step1=as.numeric(step1_set[i,])
    # step2=as.numeric(step2_set[i,])
    # deltatime=deltatime_set[i]
    # leadtime=leadtime_set[i]
    # 
    # n=cross(step1,step2)
    # av=angle_2vec(step1,step2)*pi/180/deltatime
    # leadang=av*leadtime
    # step3=rotationMatrix(leadang,n[1],n[2],n[3])[1:3,1:3]%*%step2
    # step3_mag=sqrt(step3[1,1]^2+step3[2,1]^2+step3[3,1]^2)
    result=rbind(result,t(as.data.frame(step3_dir)))
  }
  rownames(result)=c(1:nrow(step1_set))
  # result=rbind(c(NA),result,c(NA))
  return(result) 
}

omega<-function(step1_set,step2_set,deltatime_set){
  
  lengthcheck=all(sapply(list(nrow(step1_set),
                              nrow(step2_set),
                              length(deltatime_set)), 
                         function(x) x == length(deltatime_set)))
  if(!lengthcheck){
    print("length not matched")
    break()
  }
  result=data.frame(matrix(rep(NA,3), nrow=1))[numeric(0), ]
  colnames(result)=c("x","y","z")
  
  if(nrow(step1_set)==0){
    return(result)
  }
  
  for(i in 1:nrow(step1_set)){
    
    step1=as.numeric(step1_set[i,])
    step2=as.numeric(step2_set[i,])
    deltatime=deltatime_set[i]
    
    n=cross(step1,step2)
    av=angle_2vec(step1,step2)*pi/180/deltatime
    w=av*n/sqrt(n[1]^2+n[2]^2+n[3]^2)
    result=rbind(result,t(as.data.frame(w)))
  }
  colnames(result)=c("x","y","z")
  rownames(result)=c(1:nrow(step1_set))
  return(result) 
}

turnAcc<-function(step1_set,step2_set,deltatime_set,basevector){
  
  lengthcheck=all(sapply(list(nrow(step1_set),
                              nrow(step2_set),
                              length(deltatime_set),
                              nrow(basevector)), 
                         function(x) x == nrow(basevector)))
  if(!lengthcheck){
    print("length not matched")
    break()
  }
  LAresult=data.frame(matrix(rep(NA,3), nrow=1))[numeric(0), ]
  AVresult=data.frame(matrix(rep(NA,3), nrow=1))[numeric(0), ]
  
  if(nrow(step1_set)==0){
    colnames(LAresult)=colnames(AVresult)=c("x","y","z")
    return(list(LAresult,AVresult))
  }
  
  for(i in 1:nrow(step1_set)){
    
    step1=as.numeric(step1_set[i,])
    step2=as.numeric(step2_set[i,])
    deltatime=deltatime_set[i]
    base=as.numeric(basevector[i,])
    
    n=cross(step1,step2)
    av=angle_2vec(step1,step2)*pi/180/deltatime
    w=av*n/sqrt(n[1]^2+n[2]^2+n[3]^2)
    LA=cross(w,base)
    LAresult=rbind(LAresult,t(as.data.frame(LA)))
    AVresult=rbind(AVresult,t(as.data.frame(w)))
  }
  colnames(LAresult)=colnames(AVresult)=c("x","y","z")
  rownames(LAresult)=rownames(AVresult)=c(1:nrow(step1_set))
  return(list(LAresult,AVresult))
}

angle_2vec<-function(a,b){
  mag_a=sqrt(a[1]^2+a[2]^2+a[3]^2)
  mag_b=sqrt(b[1]^2+b[2]^2+b[3]^2)
  rad=acos(dot(as.matrix(a),as.matrix(b))/(mag_a*mag_b))
  deg=rad*180/pi
  return(as.numeric(deg))
}

angle_2vec_DF<-function(s_a,s_b=as.data.frame(cbind(rep(1,length=nrow(s_a)),
                                                    rep(0,length=nrow(s_a)),
                                                    rep(0,length=nrow(s_a))))){
  res=vector()
  i=1
  while(i<=nrow(s_a)){
    a=s_a[i,]
    b=s_b[i,]
    if(is.na(a[1]+a[2]+a[3]+b[1]+b[2]+b[3])){
      res=c(res,NA)
      i=i+1
      next()
    }
    mag_a=sqrt(a[1]^2+a[2]^2+a[3]^2)
    mag_b=sqrt(b[1]^2+b[2]^2+b[3]^2)
    
    outer=cross(as.matrix(a),as.matrix(b))
    inner=dot(as.matrix(a),as.matrix(b))
    rad=atan2(sqrt(outer[,1]^2+outer[,2]^2+outer[,3]^2),inner)
    
    # rad=acos(dot(as.matrix(a),as.matrix(b))/(mag_a*mag_b))
    deg=rad*180/pi
    res=c(res,as.numeric(deg))
    i=i+1
  }
  return(res)
}

func_X<-function(t,func,obj){
  na.index=which(is.na(t))
  t[na.index]=0
  x <- predict(func[[obj*3+1]],t, deriv = 0)$y
  y <- predict(func[[obj*3+2]],t, deriv = 0)$y
  z <- predict(func[[obj*3+3]],t, deriv = 0)$y
  X=data.frame(x,y,z)
  X[na.index,]=NA
  return(X)
}

func_R2<-function(t,func){
  b <- func_X(t, func, 0)
  m <- func_X(t, func, 1)
  R=m-b
  colnames(R)=c("Rx","Ry","Rz")
  return(R)
} 

func_R<-function(t,func){
  bx <- predict(func[[1]],t, deriv = 0)$y
  by <- predict(func[[2]],t, deriv = 0)$y
  bz <- predict(func[[3]],t, deriv = 0)$y
  
  mx <- predict(func[[4]],t, deriv = 0)$y
  my <- predict(func[[5]],t, deriv = 0)$y
  mz <- predict(func[[6]],t, deriv = 0)$y
  
  Rx=mx-bx
  Ry=my-by
  Rz=mz-bz
  
  return(data.frame(Rx,Ry,Rz))
} 

Three_d_plot<-function(x,y,ref,limit,range_x,range_y,range_z, model1,model2,title){
  if(rgl.cur()!=0){
    rgl.close()
  }  
  model1=c(model1,0,0)
  model2=c(model2,0,0)
  theta=seq(0,2*pi, by=pi/100)
  circlex=cos(theta)*limit
  circley=sin(theta)*limit
  circlez=rep(0,length=length(theta))
  longi=x
  lati=-y
  range=which(lati>range_x[1]&lati<range_x[2]&longi>range_y[1]&longi<range_y[2]&ref>range_z[1]&ref<range_z[2])
  plot3d(-lati[range], longi[range],ref[range],box=F,axes=F,xlab="", ylab="", zlab='',alpha=0.5,xlim=range_x,ylim=range_y,zlim=range_z)
  axes3d(labels=F,tick=F,box=F)
  par3d(windowRect = c(0, 0, 500, 500))
  aspect3d(1,1,1)
  lines3d(circlex,circley,circlez,col = rgb(0,0,0),alpha=0.4)
  segments3d(c(-limit,limit,0,0),c(0,0,limit,-limit),c(0,0,0,0),col = rgb(0,0,0),alpha=0.4)
  segments3d(c(-limit,limit,0,0),c(0,0,limit,-limit),c(range_z[2],range_z[2],range_z[2],range_z[2]), col = rgb(0,0,0),alpha=0.4)
  segments3d(c(0,0),c(0,0),c(0,range_z[2]),col = rgb(0,0,0),alpha=0.4)
  
  segments3d(c(0,0),c(0,range_z[2]),c(0,range_z[2]), lwd=3,col=hsv(222/360, 60/100, 100/100, alpha=0.9))
  segments3d(c(model2[1],range_z[2]*model2[2]+model2[1]),c(model1[1],range_z[2]*model1[2]+model1[1]),c(0,range_z[2]),lwd=3,col = hsv(100/360, 100/100, 98/100, alpha=0.9))
  
  
  view3d(theta=-120, phi=30, fov=20, zoom=0.8, scale=par3d("scale"), interactive=TRUE)
  rgl.snapshot(title)
  
  #for Blender
  scale_x=range_x[2]-range_x[1]
  scale_y=range_y[2]-range_y[1]
  scale_z=range_z[2]-range_z[1]
  blender=cbind(-lati[range]/scale_x, longi[range]/scale_y,ref[range]/scale_z)
  pulse_blender=na.omit(blender)
  # write.csv(cbind(c(1:nrow(pulse_blender)),pulse_blender),
  #           "pulse_blender.csv",row.names=F, col.names = F)
  
}

smoothed_position<-function(positiondata,func, output_fps){
  
  TIME=seq(min(positiondata$time), max(positiondata$time), by = 1/output_fps)
  result=data.frame(TIME)
  for(funcindex in 1:6){
    predicted <- predict(func[[funcindex]],TIME, deriv = 0)$y
    result=cbind(result,predicted)
  }
  
  colnames(result)=c("time", "x_bx","x_by","x_bz","x_mx","x_my","x_mz")
  
  # write.csv(result,"position_GL.csv",row.names=F, col.names = F)
  return(result)
}

add_delta_AV<-function(src,LOS,time,index1,index2,name){
  deltaT=time[index1]-time[index2]
  theta=angle_2vec_DF(LOS[index1,],LOS[index2,])
  av_mag_d=theta/deltaT/180*pi
  av_mag_d=as.data.frame(av_mag_d)
  colnames(av_mag_d)=name
  src=cbind(src,av_mag_d)
  return(src)
}

CalcSightingError<-function(sightmodel,interval,beamwidth,src,func){
  totaldelay=feedbackdelay_Sound+interval
  
  baseTime=(src$time[src$curr2])
  nextTime=(baseTime+totaldelay)
  
  baseLOS=func_R2(baseTime,func)
  nextLOS=func_R2(nextTime,func)
  
  TGTshift=angle_2vec_DF(baseLOS,nextLOS)
  av=TGTshift/totaldelay
  x=TGTshift
  pulseleadangle=predict(sightmodel,newdata=as.data.frame(x))
  AngleDiff=pulseleadangle-TGTshift
  S.Error=AngleDiff/(beamwidth/2)
  
  return(data.frame(S.Error,AngleDiff,TGTshift,av))
}

CalcSightingError2<-function(sightmodel,interval,beamwidth,src,func){
  totaldelay=feedbackdelay_Sound+interval
  
  baseTime=(src$time[src$curr2])
  nextTime=(baseTime+totaldelay)
  
  baseTGT=func_X(baseTime,func,1)
  nextBAT=func_X(nextTime,func,0)
  
  baseLOS=baseTGT-nextBAT
  nextLOS=func_R2(nextTime,func)
  
  TGTshift=angle_2vec_DF(baseLOS,nextLOS)
  av=TGTshift/totaldelay
  x=TGTshift
  pulseleadangle=predict(sightmodel,newdata=as.data.frame(x))
  
  AngleDiff=pulseleadangle-TGTshift
  S.Error=AngleDiff/(beamwidth/2)
  
  return(data.frame(S.Error,AngleDiff,TGTshift,av))
}


SightSimulation<-function(src,func){
  
  #sight direction model for each session
  x=src$LOSl_LOSc.LeadAngle
  y=src$LOSl_P.LeadAngle
  pred=lm(y~x)  #bat's predictive model
  non_pred=lm(y~0)  #non-prediction model
  dir.model=list(non_pred,pred)
  
  #rate modification
  interval=src$time[src$curr]-src$time[src$curr2]
  interval_i=rep(na.omit(interval)[1],length=length(interval))
  interval_max=rep(max(na.omit(interval)),length=length(interval))
  rate=list(interval_i,interval,interval_max)
  
  #beamwidth
  bw=src$widthHV
  bw_i=rep(na.omit(src$widthHV)[1],length=length(src$widthHV))
  bw_min=rep(min(na.omit(src$widthHV)),length=length(src$widthHV))
  beamwidth=list(bw_i,bw,bw_min)
  
  for(i in 1:length(dir.model)){
    for(j in 1:length(rate)){
      for(k in 1:length(beamwidth)){
        name=paste0("sim",i-1,j-1,k-1)  
        res=CalcSightingError(dir.model[[i]],rate[[j]],beamwidth[[k]],
                              src,func)
        colnames(res)=paste0(name,".",colnames(res))
        
        src=cbind(src,res)
      }
    }
  }
  
  return(src)
  
}


SightSimulation2<-function(src,func){
  
  #sight direction model for each session
  x=src$Tl_Tc.LeadAngle
  y=src$Tl_P.LeadAngle
  pred=lm(y~x)  #bat's predictive model
  non_pred=lm(y~0)  #non-prediction model
  dir.model=list(non_pred,pred)
  
  #rate modification
  interval=src$time[src$curr]-src$time[src$curr2]
  interval_i=rep(na.omit(interval)[1],length=length(interval))
  interval_max=rep(max(na.omit(interval)),length=length(interval))
  rate=list(interval_i,interval,interval_max)
  
  #beamwidth
  bw=src$widthHV
  bw_i=rep(na.omit(src$widthHV)[1],length=length(src$widthHV))
  bw_min=rep(min(na.omit(src$widthHV)),length=length(src$widthHV))
  beamwidth=list(bw_i,bw,bw_min)
  
  for(i in 1:length(dir.model)){
    for(j in 1:length(rate)){
      for(k in 1:length(beamwidth)){
        name=paste0("simloc",i-1,j-1,k-1)  
        res=CalcSightingError2(dir.model[[i]],rate[[j]],beamwidth[[k]],
                               src,func)
        colnames(res)=paste0(name,".",colnames(res))
        
        src=cbind(src,res)
      }
    }
  }
  
  return(src)
  
}

add_TGTshift_sim2<-function(src,interval2,func){
  #prediction model for each session
  x=src$LOSl_LOSc.LeadAngle
  y=src$LOSl_P.LeadAngle
  model=lm(y~x)
  
  #rate modification
  interval=src$time[src$curr]-src$time[src$curr2]
  interval_i=na.omit(interval)[1]
  # interval_i=max(na.omit(interval))
  
  
  totaldelay_a=feedbackdelay_Sound+interval
  totaldelay_i=feedbackdelay_Sound+interval_i
  #steps=feedbackdelay_Sound%/%interval_i+as.integer(feedbackdelay_Sound==0|feedbackdelay_Sound%/%interval_i>0)
  #totaldelay_a=src$time[src$curr]-src$time[src$last]
  #totaldelay_i=interval_i*steps
  
  baseTime=(src$time[src$curr2])
  nextTime_a=(baseTime+totaldelay_a)
  nextTime_i=(baseTime+totaldelay_i)
  baseLOS=func_R2(baseTime,func)
  nextLOS_i=func_R2(nextTime_i,func)
  nextLOS_a=func_R2(nextTime_a,func)
  
  baseLOS[which(nextTime_i>0|nextTime_a>0),]=c(NA,NA,NA)
  nextLOS_i[which(nextTime_i>0|nextTime_a>0),]=c(NA,NA,NA)
  nextLOS_a[which(nextTime_i>0|nextTime_a>0),]=c(NA,NA,NA)
  
  TGTshift_i=angle_2vec_DF(baseLOS,nextLOS_i)
  TGTshift_a=angle_2vec_DF(baseLOS,nextLOS_a)
  anglediff_i=TGTshift_i-predict(model,newdata=as.data.frame(TGTshift_i))
  anglediff_a=TGTshift_a-predict(model,newdata=as.data.frame(TGTshift_a))
  
  av_i=TGTshift_i/totaldelay_i
  av_a=TGTshift_a/totaldelay_a
  
  
  baseTGT=func_X(baseTime,func,1)
  nextBAT_i=func_X(nextTime_i,func,0)
  nextBAT_a=func_X(nextTime_a,func,0)
  baseLOSloc_i=baseTGT-nextBAT_i
  baseLOSloc_a=baseTGT-nextBAT_a
  
  baseLOSloc_i[which(nextTime_i>0|nextTime_a>0),]=c(NA,NA,NA)
  baseLOSloc_a[which(nextTime_i>0|nextTime_a>0),]=c(NA,NA,NA)
  
  TGTshiftloc_i=angle_2vec_DF(baseLOSloc_i,nextLOS_i)
  TGTshiftloc_a=angle_2vec_DF(baseLOSloc_a,nextLOS_a)
  
  src=cbind(src,
            totaldelay_i,
            totaldelay_a,
            TGTshift_i,
            TGTshift_a,
            anglediff_i,
            anglediff_a,
            av_i,
            av_a,
            TGTshiftloc_i,
            TGTshiftloc_a
  )
  
  return(src)
}


initialformatting<-function(posi,pul){
  posi=na.omit(posi)
  
  if(any(posi$mx==0&posi$my==0&posi$mz==0)){
    capture_step=min(which(posi$mx==0&posi$my==0&posi$mz==0))
  }else{
    capture_step=length(posi$mx)
  }
  
  ##set time as "time to capture"
  capture_time=posi$time[capture_step]
  posi$time=posi$time-capture_time
  pul$time=pul$time-capture_time
  
  ##focusing on BEFORE capture
  posi=posi[1:(capture_step-1),]
  
  ## matchinng to position data
  pul=pul[which(pul$time>min(posi$time)),]
  pul=pul[which(pul$time<=0),]
  
  
  ##save for GL
  # setwd(subdir)
  # write.csv(pulse,"pulse_GL.csv",row.names=F, col.names = F)
  
  
  ##set initial bat position as principle point
  for(dim in 1:3){
    basepoint=posi[1,dim+1]/1000.0
    posi[,dim+1]=posi[,dim+1]/1000.0-basepoint
    posi[,dim+4]=posi[,dim+4]/1000.0-basepoint
  }
  
  #add beam width
  pul$widthH[which(pul$widthH==0)]=NA
  pul$widthV[which(pul$widthV==0)]=NA
  pul=na.omit(pul)
  pul$widthH=pul$widthH*2
  pul$widthV=pul$widthV*2
  
  widthHV=sqrt(pul$widthH^2+pul$widthV^2)
  widthHV_init=sqrt(pul$widthH^2+pul$widthV^2)[1]
  pul=cbind(pul,widthHV,widthHV_init)
  
  #clockwise positive
  pul$peakdirectionH=pul$peakdirectionH*(-1)
  pul$cobH=pul$cobH*(-1)
  
  #add pulse interval
  interval=(pul[-1,]-pul[-nrow(pul),])$time
  pul=cbind(pul, interval=c(NA,interval))
  
  return(list(posi,pul))
}


getRW<-function(bat,tgt,interval){
  R=tgt-bat
  Rmag=sqrt(R[,1]^2+R[,2]^2+R[,3]^2)
  Wmag=angle_2vec_DF(R[-1,],R[-nrow(R),])/interval
  Wmag=c(NA,Wmag)
  
  return(data.frame(Rmag,Wmag))
}

oripursuit<-function(dataset){
  interval=dataset$time[-1]-dataset$time[-nrow(dataset)]
  sns_b <- data.frame(x=dataset$x_bx,y=dataset$x_by,z=dataset$x_bz)
  sns_m<-data.frame(x=dataset$x_mx,y=dataset$x_my,z=dataset$x_mz)
  mnv_b <- data.frame(x=dataset$x_bx_L,y=dataset$x_by_L,z=dataset$x_bz_L)
  mnv_m <- data.frame(x=dataset$x_mx_L,y=dataset$x_my_L,z=dataset$x_mz_L)
  sns_RW=getRW(sns_b,sns_m,interval)
  mnv_RW=getRW(mnv_b,mnv_m,interval)
  
  
  sn=cbind(sns_b,sns_RW)
  select=sn$Rmag>dist_threshold
  sn[which(!(select)),]=NA
  
  timeselect=(dataset$time+feedbackdelay_Flight)<0
  mn=cbind(mnv_b,mnv_RW)
  select=mn$Rmag>dist_threshold
  mn[which(!(timeselect&select)),]=NA
  return(cbind(sn=sn,mn=mn))  
}

purepursuit<-function(dataset,sensingmode=0){
  
  interval=dataset$time[-1]-dataset$time[-nrow(dataset)]
  sns_b <- data.frame(x=dataset$x_bx,y=dataset$x_by,z=dataset$x_bz)
  sns_m<-data.frame(x=dataset$x_mx,y=dataset$x_my,z=dataset$x_mz)
  mnv_b <- data.frame(x=dataset$x_bx_L,y=dataset$x_by_L,z=dataset$x_bz_L)
  mnv_m <- data.frame(x=dataset$x_mx_L,y=dataset$x_my_L,z=dataset$x_mz_L)
  
  
  deltaloc=mnv_b[-1,]-mnv_b[-nrow(mnv_b),]
  deltaloc_mag=sqrt(deltaloc[,1]^2+deltaloc[,2]^2+deltaloc[,3]^2)
  v=deltaloc_mag/interval
  
  sim_result=data.frame(sn.x=sns_b$x[1],sn.y=sns_b$y[1],sn.z=sns_b$z[1],
                        mn.x=mnv_b$x[1],mn.y=mnv_b$y[1],mn.z=mnv_b$z[1])
  if(nrow(mnv_b)>1){
    for(i in 1:(nrow(mnv_b)-1)){  
      
      sensing_s=(sim_result[i,1:3])
      mnvring_s=(sim_result[i,4:6])
      
      if(sensingmode){
        sensing_s=mnvring_s
      }
      
      sim_R=sns_m[i,]-sensing_s
      sim_Rmag=as.numeric(sqrt(sim_R[1]^2+sim_R[2]^2+sim_R[3]^2))
      sim_TGTdir=sim_R/sim_Rmag
      
      next_sensing_s=mnvring_s+sim_TGTdir*(v[i]*(interval[i]-feedbackdelay_Flight))
      next_mnvring_s=mnvring_s+sim_TGTdir*(v[i]*(interval[i]-0))
      
      result<-c(as.vector(as.matrix(next_sensing_s)),
                as.vector(as.matrix(next_mnvring_s)))
      sim_result=rbind(sim_result,result)
    }
  }else{sim_result=data.frame(NA,NA,NA,NA,NA,NA)}
  
  
  sns_b=sim_result[,1:3]
  mnv_b=sim_result[,4:6]
  sns_RW=getRW(sns_b,sns_m,interval)
  mnv_RW=getRW(mnv_b,mnv_m,interval)
  
  sn=cbind(sns_b,sns_RW)
  colnames(sn)=c("x","y","z","Rmag","Wmag")
  select=sn$Rmag>dist_threshold
  sn[which(!(select)),]=NA
  
  mn=cbind(mnv_b,mnv_RW)
  colnames(mn)=c("x","y","z","Rmag","Wmag")
  timeselect=(dataset$time+feedbackdelay_Flight)<0
  select=mn$Rmag>dist_threshold
  mn[which(!(timeselect&select)),]=NA
  return(cbind(sn=sn,mn=mn))  
}

leadpursuit<-function(dataset,sensingmode=0){
  interval=dataset$time[-1]-dataset$time[-nrow(dataset)]
  sns_b <- data.frame(x=dataset$x_bx,y=dataset$x_by,z=dataset$x_bz)
  sns_m<-data.frame(x=dataset$x_mx,y=dataset$x_my,z=dataset$x_mz)
  mnv_b <- data.frame(x=dataset$x_bx_L,y=dataset$x_by_L,z=dataset$x_bz_L)
  mnv_m <- data.frame(x=dataset$x_mx_L,y=dataset$x_my_L,z=dataset$x_mz_L)
  deltaloc=mnv_b[-1,]-mnv_b[-nrow(mnv_b),]
  deltaloc_mag=sqrt(deltaloc[,1]^2+deltaloc[,2]^2+deltaloc[,3]^2)
  v=deltaloc_mag/interval
  
  
  if(nrow(mnv_b)>2){
    sim_result=data.frame(sn.x=sns_b$x[1:2],sn.y=sns_b$y[1:2],sn.z=sns_b$z[1:2],
                          mn.x=mnv_b$x[1:2],mn.y=mnv_b$y[1:2],mn.z=mnv_b$z[1:2])
    for(i in 2:(nrow(mnv_b)-1)){  
      LOS=sns_m[(i-1):i,]-sim_result[(i-1):i,1:3]
      sim_TGTdir=angular_pred_core(LOS[1,],LOS[2,],interval[i-1],0.1535,12.79/180*pi)
      
      mnvring_s=sim_result[i,4:6]
      next_sensing_s=mnvring_s+sim_TGTdir*(v[i]*(interval[i]-feedbackdelay_Flight))
      next_mnvring_s=mnvring_s+sim_TGTdir*(v[i]*(interval[i]-0))
      
      result<-c(as.vector(as.matrix(next_sensing_s)),
                as.vector(as.matrix(next_mnvring_s)))
      sim_result=rbind(sim_result,result)
    }
  }else{sim_result=data.frame(NA,NA,NA,NA,NA,NA)}
  
  sns_b=sim_result[,1:3]
  mnv_b=sim_result[,4:6]
  sns_RW=getRW(sns_b,sns_m,interval)
  mnv_RW=getRW(mnv_b,mnv_m,interval)
  
  sn=cbind(sns_b,sns_RW)
  colnames(sn)=c("x","y","z","Rmag","Wmag")
  select=sn$Rmag>dist_threshold
  sn[which(!(select)),]=NA
  
  mn=cbind(mnv_b,mnv_RW)
  colnames(mn)=c("x","y","z","Rmag","Wmag")
  timeselect=(dataset$time+feedbackdelay_Flight)<0
  select=mn$Rmag>dist_threshold
  mn[which(!(timeselect&select)),]=NA
  return(cbind(sn=sn,mn=mn))  
}

leadpursuit_abs<-function(dataset,sensingmode=0){
  interval=dataset$time[-1]-dataset$time[-nrow(dataset)]
  sns_b <- data.frame(x=dataset$x_bx,y=dataset$x_by,z=dataset$x_bz)
  sns_m<-data.frame(x=dataset$x_mx,y=dataset$x_my,z=dataset$x_mz)
  mnv_b <- data.frame(x=dataset$x_bx_L,y=dataset$x_by_L,z=dataset$x_bz_L)
  mnv_m <- data.frame(x=dataset$x_mx_L,y=dataset$x_my_L,z=dataset$x_mz_L)
  deltaloc=mnv_b[-1,]-mnv_b[-nrow(mnv_b),]
  deltaloc_mag=sqrt(deltaloc[,1]^2+deltaloc[,2]^2+deltaloc[,3]^2)
  v=deltaloc_mag/interval
  
  
  
  if(nrow(mnv_b)>2){
    sim_result=data.frame(sn.x=sns_b$x[1:2],sn.y=sns_b$y[1:2],sn.z=sns_b$z[1:2],
                          mn.x=mnv_b$x[1:2],mn.y=mnv_b$y[1:2],mn.z=mnv_b$z[1:2])
    for(i in 2:(nrow(mnv_b)-1)){  
      LOS=sns_m[(i-1):i,]-sim_result[(i-1):i,1:3]
      mothv=sns_m[i,]-sns_m[i-1,]
      
      sim_TGTdir=absTGTV_pred_core(LOS[1,],LOS[2,],mothv,interval[i-1],0.1535,12.79/180*pi)
      
      mnvring_s=sim_result[i,4:6]
      next_sensing_s=mnvring_s+sim_TGTdir*(v[i]*(interval[i]-feedbackdelay_Flight))
      next_mnvring_s=mnvring_s+sim_TGTdir*(v[i]*(interval[i]-0))
      
      result<-c(as.vector(as.matrix(next_sensing_s)),
                as.vector(as.matrix(next_mnvring_s)))
      sim_result=rbind(sim_result,result)
    }
  }else{sim_result=data.frame(NA,NA,NA,NA,NA,NA)}
  
  sns_b=sim_result[,1:3]
  mnv_b=sim_result[,4:6]
  sns_RW=getRW(sns_b,sns_m,interval)
  mnv_RW=getRW(mnv_b,mnv_m,interval)
  
  sn=cbind(sns_b,sns_RW)
  colnames(sn)=c("x","y","z","Rmag","Wmag")
  select=sn$Rmag>dist_threshold
  sn[which(!(select)),]=NA
  
  mn=cbind(mnv_b,mnv_RW)
  colnames(mn)=c("x","y","z","Rmag","Wmag")
  timeselect=(dataset$time+feedbackdelay_Flight)<0
  select=mn$Rmag>dist_threshold
  mn[which(!(timeselect&select)),]=NA
  return(cbind(sn=sn,mn=mn))  
}



FlightSimulation<-function(pulse, sensingmode=0){
  
  
  timerange=pulse$time+feedbackdelay_Flight
  ori_p=oripursuit(pulse)
  pure_p=purepursuit(pulse,sensingmode)
  lead_p=leadpursuit(pulse,sensingmode)
  labs_p=leadpursuit_abs(pulse,sensingmode)
  
  result=data.frame(ori_p=ori_p,pure_p=pure_p,lead_p=lead_p,labs_p=labs_p)
  EachFlight=result
  return(EachFlight)
}

addKinematicData<-function(pulse,func){
  location <- data.frame(matrix(rep(NA, 1), nrow=nrow(pulse)))[numeric(0), ]
  for(funcindex in 1:6){
    predicted <- predict(func[[funcindex]],pulse$time, deriv = 0)$y
    location=as.data.frame(cbind(location,predicted))
  }
  
  for(funcindex in 1:6){
    predicted <- predict(func[[funcindex]],pulse$time+feedbackdelay_Flight, deriv = 0)$y
    location=as.data.frame(cbind(location,predicted))
  }
  
  colnames(location)=c("x_bx","x_by","x_bz","x_mx","x_my","x_mz",
                       "x_bx_L","x_by_L","x_bz_L","x_mx_L","x_my_L","x_mz_L")
  
  
  velocity <- data.frame(matrix(rep(NA, 1), nrow=nrow(pulse)))[numeric(0), ]
  for(funcindex in 1:6){
    predicted <- predict(func[[funcindex]],pulse$time, deriv = 1)$y
    velocity=as.data.frame(cbind(velocity,predicted))
  }
  for(funcindex in 1:6){
    predicted <- predict(func[[funcindex]],pulse$time+feedbackdelay_Flight, deriv = 1)$y
    velocity=as.data.frame(cbind(velocity,predicted))
  }
  colnames(velocity)=c("v_bx","v_by","v_bz","v_mx","v_my","v_mz",
                       "v_bx_L","v_by_L","v_bz_L","v_mx_L","v_my_L","v_mz_L")
  
  acceleration <- data.frame(matrix(rep(NA, 1), nrow=nrow(pulse)))[numeric(0), ]
  for(funcindex in 1:6){
    predicted <- predict(func[[funcindex]],pulse$time, deriv = 2)$y
    acceleration=as.data.frame(cbind(acceleration,predicted))
  }
  colnames(acceleration)=c("a_bx","a_by","a_bz","a_mx","a_my","a_mz")
  
  
  
  pulse=cbind(pulse,location,velocity,acceleration)
  
  
  
  ## linear vector 
  R=data.frame(x=(pulse$x_mx-pulse$x_bx), 
               y=(pulse$x_my-pulse$x_by), 
               z=(pulse$x_mz-pulse$x_bz))
  R_mag=sqrt(R$x^2+R$y^2+R$z^2)
  
  R_L=data.frame(x=(pulse$x_mx_L-pulse$x_bx_L), 
                 y=(pulse$x_my_L-pulse$x_by_L), 
                 z=(pulse$x_mz_L-pulse$x_bz_L))
  
  ## linear velocity vector (instantaneous)
  Vx=velocity$v_mx-velocity$v_bx
  Vy=velocity$v_my-velocity$v_by
  Vz=velocity$v_mz-velocity$v_bz
  V=data.frame(x=Vx,y=Vy,z=Vz)
  V_mag=sqrt(V$x^2+V$y^2+V$z^2)
  
  ## angular velocity vector (instantaneous)
  av=cross(as.matrix(R),as.matrix(V))/(R_mag^2)
  av=data.frame(av)
  colnames(av)=c("x","y","z")
  av_mag=sqrt(av$x^2+av$y^2+av$z^2)*180/pi
  
  pulse=cbind(pulse,R=R,R_mag,V_mag,av_mag,R_L=R_L)
  
  TGTV=data.frame(x=velocity$v_mx,y=velocity$v_my,z=velocity$v_mz)
  absoluteTGTav=cross(as.matrix(R),as.matrix(TGTV))/(R_mag^2)
  absoluteTGTav=data.frame(absoluteTGTav)
  colnames(absoluteTGTav)=c("x","y","z")
  
  abTGTav_mag=sqrt(absoluteTGTav$x^2+absoluteTGTav$y^2+absoluteTGTav$z^2)*180/pi
  pulse=cbind(pulse,abTGTav_mag)
  
  return(pulse)
}

CalcAV.vector<-function(vector,interval){
  # lengthcheck=all(sapply(list(nrow(vector),
  #                             length(interval)),
  #                        function(x) x == length(interval)))
  # if(!lengthcheck){
  #   print("length not matched")
  #   break()
  # }
  result=data.frame(matrix(rep(NA,3), nrow=1))[numeric(0), ]
  colnames(result)=c("x","y","z")
  
  if(nrow(vector)==0){
    return(result)
  }
  
  vector1=vector[-nrow(vector),]
  vector2=vector[-1,]
  
  for(i in 1:nrow(vector1)){
    step1=as.numeric(vector1[i,])
    step2=as.numeric(vector2[i,])
    n=cross(step1,step2)
    av=angle_2vec(step1,step2)*pi/180/interval[i]
    av.vector=av*n/sqrt(n[1]^2+n[2]^2+n[3]^2)
    result=rbind(result,t(as.data.frame(av.vector)))
  }
  rownames(result)=c(1:nrow(vector1))
  
  return(result)
}

FlightAnalysis<-function(dataset){
  EachSession=data.frame()
  for(j in 1:length(dataset)){
    p=dataset[[j]]
    result=data.frame(index=c(1:nrow(p)))
    
    LOS=data.frame(x=p$R.x,y=p$R.y,z=p$R.z)
    interval=p$time[p$curr]-p$time[p$curr2]
    
    BatXL=data.frame(x=p$x_bx_L,y=p$x_by_L,z=p$x_bz_L)
    BatVDL=(BatXL[p$curr,]-BatXL[p$curr2,])/interval
    BatV=data.frame(x=p$v_bx,y=p$v_by,z=p$v_bz)
    BatVL=data.frame(x=p$v_bx_L,y=p$v_by_L,z=p$v_bz_L)
    
    MothX=data.frame(x=p$x_mx,y=p$x_my,z=p$x_mz)
    MothV=data.frame(x=p$v_mx,y=p$v_my,z=p$v_mz)
    MothVD=(MothX[p$curr,]-MothX[p$curr2,])/interval
    
    LOS_Tl=MothX-BatXL
    
    ## Bat's flight direction relative to TGT direction
    
    result=FlightDirection(LOS[p$curr2,],BatVDL[p$curr,],MothVD[p$curr,],result,"LOSl_F_V")
    result=FlightDirection(LOS_Tl[p$curr2,],BatVDL[p$curr,],MothVD[p$curr,],result,"Tl_F_V")
    
    result=FlightLA(BatVDL,LOS,interval,result,p,"vbat")
    
    result=cbind(result,TGTdir=FlightSimulation(p,0))
    result=cbind(result,TGTloc=FlightSimulation(p,1))
    
    ## storing to EachSession
    EachSession=rbind(EachSession,result)
  }
  return(EachSession)
  
}

FlightDirection<-function(non_pred, dir,ref,res, name){
  a=rotationALL_polar(non_pred,
                      dir,
                      non_pred+ref,
                      180,100,1)
  a=add_NA_head(a,nrow(res))
  b=a$r2$OBJ$hv
  c=a$s$OBJ$hv
  colnames(b)=paste0(name,".",colnames(b))
  colnames(c)=paste0(name,"_scaled",".",colnames(c))
  res=cbind(res,b,c)
  return(res)
}

FlightLA<-function(v.vector, los, intvl,res,src, name){
  #x,o,o,...,x
  turn_v.vector=turnAcc(v.vector[src$curr,],
                        v.vector[src$curr2,],
                        intvl[src$curr],
                        v.vector[src$curr2,])
  LA_v.vector=turn_v.vector[[1]]
  LA_v.vector=add_NA_head2(LA_v.vector,nrow(res))
  
  #x,o,o,...,o
  turn_los=turnAcc(los[src$curr,], 
                   los[src$curr2,], 
                   intvl[src$curr], 
                   v.vector[src$curr2,])
  LA_los=turn_los[[1]]
  LA_los=add_NA_head2(LA_los,nrow(res))
  
  v.vector=add_NA_head2(v.vector,nrow(res))
  steering_LA=rotationALL_polar(v.vector[src$curr2,],
                                LA_v.vector,
                                LA_los,
                                180,100,1)
  a=add_NA_tail(steering_LA,nrow(res))
  b=as.data.frame(a$r2$OBJ$xyz_R$z)
  c=as.data.frame(a$r2$REF$xyz_R$z)
  d=as.data.frame(sqrt(LA_v.vector$x^2+LA_v.vector$y^2+LA_v.vector$z^2))
  
  colnames(b)=c(paste0("LA.",name,".projected"))
  colnames(c)=c(paste0("LA.",name,".LOS"))
  colnames(d)=c(paste0("LA.",name,".mag"))
  
  res=cbind(res,b,c,d)
  return(res)
}

getLastDetectedTimestep<-function(p,sensing_feedback_delay){
  
  t=p$time
  
  sensed.t=numeric()
  sensed.index=numeric()
  current.index=c(1:length(t))
  curr2=current.index-1
  curr2[which(curr2<=0)]=NA
  curr3=current.index-2
  curr3[which(curr3<=0)]=NA
  
  if(length(t)<2){
    for(i in 1:length(t)){
      sensed.t[i]=NA
      sensed.index[i]=NA
    }
  }else{
    
    for(i in length(t):2){
      curr.t=t[i]
      delay.t=curr.t-sensing_feedback_delay
      
      for(j in i:2){
        
        if(t[j]>=delay.t&delay.t>t[j-1]){
          sensed.t[i]=t[j-1]
          sensed.index[i]=j-1
          break
        }else{
          sensed.t[i]=NA
          sensed.index[i]=NA
        }
      }
    }
  }
  
  
  last2=sensed.index-1
  last2[which(last2==0)]=NA
  p=cbind(p,curr=current.index,
          curr2=curr2,
          curr3=curr3,
          last=sensed.index,
          last2=last2)
  return(p)
  
}

SightAnalysis<-function(dataset){
  EachSession=data.frame()
  
  for(j in 1:length(dataset)){
    p=dataset[[j]]
    print(paste0(head3,"-",j, ", num of data:", nrow(p)))
    #basic information
    time=p$time
    interval2=p[p$curr,]$time-p[p$last,]$time
    LOS=data.frame(x=p$R.x,y=p$R.y,z=p$R.z)
    
    BatX=data.frame(x=p$x_bx,y=p$x_by,z=p$x_bz)
    TGTX=data.frame(x=p$x_mx,y=p$x_my,z=p$x_mz)
    
    LOS_Tl=TGTX[p$last,]-BatX[p$curr,]
    
    PD_hv<-data.frame(p$peakdirectionH,p$peakdirectionV)
    PD=cnvrtAngle2Vector(PD_hv)[,3:5]
    
    p=add_delta_AV(p,LOS,time,p$last,p$last2,"perceived_av_Sound")
    p=add_delta_AV(p,LOS,time,p$curr2,p$curr3,"perceived_av_Flight")
    
    ## naming rule
    ## P:pulse directon
    ## L:linear predicted direction
    ## A:angular predicted direction
    ## LOSc:Current LOS
    ## LOSl:Last LOS
    
    ## pulse and TGT direction in global coordinates
    Null_P_LOSc=rotationALL_polar(data.frame(rep(0, nrow(PD)),rep(0, nrow(PD)),rep(0, nrow(PD))),
                                  PD,
                                  LOS,
                                  180,100,1)
    p=cbind(p,"P"=Null_P_LOSc$v$OBJ$hv)
    p=cbind(p,"LOSc"=Null_P_LOSc$v$REF$hv)
    
    
    LOSl_P_LOSc=rotationALL_polar(LOS[p$last,],
                                  PD[p$curr,],
                                  LOS[p$curr,],
                                  180,100,1)
    p=cbind(p,"LOSl_P"=LOSl_P_LOSc$r2$OBJ$hv)
    p=cbind(p,"LOSl_P_scaled"=LOSl_P_LOSc$s$OBJ$hv)
    p=cbind(p,"LOSl_LOSc"=LOSl_P_LOSc$r2$REF$hv)
    p=cbind(p,"P(-1)"=LOSl_P_LOSc$v$OBJ$hv)
    p=cbind(p,"LOSc(-1)"=LOSl_P_LOSc$v$REF$hv)
    # 
    
    ## pulse direction relative to LAST TGT LOCATION (data for figEx.2)
    Tl_P_Tc=rotationALL_polar(LOS_Tl[p$curr,],
                              PD[p$curr,],
                              LOS[p$curr,],
                              180,100,1)
    Tl_P_Tc=add_NA_tail(Tl_P_Tc,nrow(p))
    p=cbind(p,"Tl_P"=Tl_P_Tc$r2$OBJ$hv)
    p=cbind(p,"Tl_P_scaled"=Tl_P_Tc$s$OBJ$hv)
    p=cbind(p,"Tl_Tc"=Tl_P_Tc$r2$REF$hv)
    
    
    ## Predicted direction (linear and angular) 
    
    
    
    
    last1_2_interval=p$time[p$last]-p$time[p$last2]
    curr_last_interval=p$time[p$curr]-p$time[p$last]
    PRED_LV=linear_pred(LOS[p$last2,],
                        LOS[p$last,],
                        last1_2_interval,
                        curr_last_interval)
    PRED_LV=add_NA_head2(PRED_LV,nrow(p))
    
    
    PRED_AV=angular_pred(LOS[p$last2,],
                         LOS[p$last,],
                         last1_2_interval,
                         curr_last_interval)
    PRED_AV=add_NA_head2(PRED_AV,nrow(p))
    
    
    P_L_A=rotationALL_polar(PD,
                            PRED_LV,
                            PRED_AV,
                            180,100,1)
    
    p=cbind(p,"P_L"=P_L_A$r2$OBJ$hv)
    p=cbind(p,"P_A"=P_L_A$r2$REF$hv)
    
    
    LOSc_L_A=rotationALL_polar(LOS,
                               PRED_LV,
                               PRED_AV,
                               180,100,1)
    
    LOSc_P_Null=rotationALL_polar(LOS,
                                  PD,
                                  data.frame(rep(0, nrow(PD)),rep(0, nrow(PD)),rep(0, nrow(PD))),
                                  180,100,1)
    
    p=cbind(p,"LOSc_L"=LOSc_L_A$r2$OBJ$hv)
    p=cbind(p,"LOSc_A"=LOSc_L_A$r2$REF$hv)
    p=cbind(p,"LOSc_P"=LOSc_P_Null$r2$OBJ$hv)
    
    
    # preparation for the simulation 
    # TGT direction shift caused by factor2 (interval) 
    # _a:actual (+), _i:fixed to initial value (-)
    p=add_TGTshift_sim2(p,interval2,motionFunc)
    p=SightSimulation(p, motionFunc)
    p=SightSimulation2(p, motionFunc)
    
    ## storing to EachSession
    EachSession=rbind(EachSession,p)
  }
  
  return(EachSession)
}

huelisting<-function(y_list){
  colis=numeric()
  for(i in 0:(length(y_list)-1)){
    colis=c(colis,i/length(y_list))
  }
  return(colis)
}

multi_paired_comp<-function(ylist,
                            id,ex,
                            ylimit=c(0,0),xlab,ylab,
                            col=NULL,backON="backOFF",
                            title="",F_dir=".",T_dir="."){
  
  y_num=length(ylist)
  col_back=character()
  n_combi=y_num*(y_num-1)/2
  if(is.null(col)){
    col=character()
    h=huelisting(ylist)
    for(i in 1:n_combi){
      col=c(col,hsv(h[i],0.65,1,alpha=1))
      col_back=c(col_back,hsv(h[i],0.65,0.8,alpha=0.1))
    }    
  }else{
    while(length(col)<n_combi){
      col=c(col,col[length(col)])
    }
    col_back=rep(hsv(0,0, 0, alpha=0.15),length(n_combi))
  }
  
  xlist=list()
  for(i in 1:y_num){
    x_=rep(i-1,length(ylist[[i]]))
    xlist=c(xlist,list(x_))
  }
  
  seglist=list()
  lmlist=list()
  plist=list()
  each_plist=list()
  meanlist=numeric()
  sizelist=c(rep(7.5,y_num-2),6,5)
  
  loop=1
  resultset=data.frame(matrix(rep(NA,5), nrow=1))[numeric(0), ]
  for(j in 1:y_num){
    for(i in 1:y_num){
      if(i>=j)next
      plot=data.frame(a.x=xlist[[i]],
                      a.y=ylist[[i]],
                      b.x=xlist[[j]],
                      b.y=ylist[[j]])
      plot<-as.data.frame(plot[is.finite(rowSums(plot)),])
      
      seg<-geom_segment(aes(x = a.x, y = a.y, xend = b.x, yend =  b.y), 
                        colour = col_back[loop],data = plot)
      seglist=c(seglist,list(seg))
      
      point_a<-geom_point(aes(x = a.x, y = a.y), 
                          colour = col_back[loop],size=1, data=plot)
      point_b<-geom_point(aes(x = b.x, y = b.y), 
                          colour = col_back[loop],size=1, data=plot)
      each_plist<-c(each_plist,point_a,point_b)
      
      
      res=paired_comparison(ylist[i],ylist[j],id,ex)
      result=cbind("pair"=paste0(i-1,"-",j-1),res)
      resultset=rbind(resultset,result)
      
      x1=xlist[[i]][1]
      x2=xlist[[j]][1]
      y1=res$b
      y2=res$b+res$a
      plotLMM=data.frame(x1=x1,y1=y1,x2=x2,y2=y2)
      plot_p=data.frame(x=c(x1,x2),y=c(y1,y2))
      
      meanlist=c(meanlist,c(y1,y2))
      lm<-geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
                       colour = col[loop],size=2.5, data=plotLMM)
      point<-geom_point(aes(x = x, y = y), 
                        colour = col[loop],size=sizelist[loop], data=plot_p)
      
      plist<-c(plist,list(point))
      lmlist<-c(lmlist,list(lm))
      loop=loop+1
    }
  }
  
  print(resultset)
  
  
  draw<-ggplot()
  
  draw<-draw+seglist
  for(i in 1:length(lmlist)){
    draw<-draw+lmlist[[i]]+plist[[i]]
  }
  
  base=mean(meanlist)
  range=abs(max(meanlist)-min(meanlist))*4
  ymin=base-range
  ymax=base+range
  
  if(ylimit[1]!=0||ylimit[2]!=0){
    ymin=ylimit[1]
    ymax=ylimit[2]
  }
  
  draw<-draw+coord_cartesian(xlim=c(-1,y_num),ylim=c(ymin,ymax))
  draw<-draw+xlab(xlab)+ylab(ylab)
  
  if(backON=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA), 
      axis.text= element_blank(),
      plot.title = element_blank(),
      axis.title = element_blank()
    )
  }
  
  if(!is.null(T_dir)){
    write.csv(resultset,paste0(T_dir,"/",title,".csv"))
  }
  if(!is.null(F_dir)){
    ggsave(paste0(F_dir,"/",title,".png"), draw, bg = "transparent",dpi=dpi2,width=11.3,height=10)
  }
  return(draw)
}

paired_comp_elements<-function(ylist,
                               id,ex,col=NULL,col_back=NULL,backON="backOFF"){
  
  y_num=length(ylist)
  n_combi=y_num*(y_num-1)/2
  if(is.null(col)){
    col=character()
    col_back=character()
    h=huelisting(ylist)
    for(i in 1:n_combi){
      col=c(col,hsv(h[i],0.65,1,alpha=1))
      col_back=c(col_back,hsv(h[i],0.65,0.8,alpha=0.1))
    }    
  }else{
    while(length(col)<n_combi){
      col=c(col,col[length(col)])
      col_back=c(col_back,col[length(col_back)])
    }
  }
  
  xlist=list()
  for(i in 1:y_num){
    x_=rep(i-1,length(ylist[[i]]))
    xlist=c(xlist,list(x_))
  }
  
  seglist=list()
  lmlist=list()
  plist=list()
  each_plist=list()
  meanlist=numeric()
  sizelist=c(rep(7.5,y_num-2),6,5)
  
  loop=1
  resultset=data.frame(matrix(rep(NA,5), nrow=1))[numeric(0), ]
  for(j in 1:y_num){
    for(i in 1:y_num){
      if(i>=j)next
      plot=data.frame(a.x=xlist[[i]],
                      a.y=ylist[[i]],
                      b.x=xlist[[j]],
                      b.y=ylist[[j]])
      plot<-as.data.frame(plot[is.finite(rowSums(plot)),])
      
      seg<-geom_segment(aes(x = a.x, y = a.y, xend = b.x, yend =  b.y), 
                        colour = col_back[loop],data = plot)
      seglist=c(seglist,list(seg))
      
      point_a<-geom_point(aes(x = a.x, y = a.y), 
                          colour = col_back[loop],size=1, data=plot)
      point_b<-geom_point(aes(x = b.x, y = b.y), 
                          colour = col_back[loop],size=1, data=plot)
      each_plist<-c(each_plist,point_a,point_b)
      
      
      res=paired_comparison(ylist[i],ylist[j],id,ex)
      result=cbind("pair"=paste0(i-1,"-",j-1),res)
      resultset=rbind(resultset,result)
      
      x1=xlist[[i]][1]
      x2=xlist[[j]][1]
      y1=res$b
      y2=res$b+res$a
      plotLMM=data.frame(x1=x1,y1=y1,x2=x2,y2=y2)
      plot_p=data.frame(x=c(x1,x2),y=c(y1,y2))
      
      meanlist=c(meanlist,c(y1,y2))
      lm<-geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
                       colour = col[loop],size=2.5, data=plotLMM)
      point<-geom_point(aes(x = x, y = y), 
                        colour = col[loop],size=sizelist[loop], data=plot_p)
      
      plist<-c(plist,list(point))
      lmlist<-c(lmlist,list(lm))
      loop=loop+1
    }
  }
  
  return(list(seglist,lmlist,plist,resultset))
}

Overlayed_paired_comp<-function(ylist1,ylist2,
                                id,ex,
                                ylimit=c(0,0),xlab,ylab1,ylab2,
                                col1=NULL,col2=NULL,col1_back=NULL,col2_back=NULL,backON="backOFF",
                                title="",F_dir=".",T_dir="."){
  
  list1=paired_comp_elements(ylist1,id,ex,col1,col1_back,backON)
  list2=paired_comp_elements(ylist2,id,ex,col2,col2_back,backON)
  
  
  draw<-ggplot()
  
  draw<-draw+list1[[1]]+list2[[1]]
  for(i in 1:length(list1[[2]])){
    draw<-draw+list1[[2]][[i]]+list1[[3]][[i]]+list2[[2]][[i]]+list2[[3]][[i]]
  }
  
  
  ymin=ylimit[1]
  ymax=ylimit[2]
  
  
  draw<-draw+coord_cartesian(xlim=c(-1,length(ylist1)),ylim=c(ymin,ymax))
  draw<-draw+xlab(xlab)+ylab(ylab1)
  
  if(backON=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA), 
      axis.text= element_blank(),
      plot.title = element_blank(),
      axis.title = element_blank()
    )
  }
  
  if(!is.null(T_dir)){
    
    lab1=as.character(eval(substitute(alist(col1))))
    lab2=as.character(eval(substitute(alist(col2))))
    
    write.csv(list1[[4]],paste0(T_dir,"/",paste0(title,"-",lab1),".csv"))
    write.csv(list2[[4]],paste0(T_dir,"/",paste0(title,"-",lab2),".csv"))
  }
  if(!is.null(F_dir)){
    ggsave(paste0(F_dir,"/",title,".png"), draw, bg = "transparent",dpi=dpi2,width=11.3,height=10)
  }
  return(draw)
}


LMM_comparison<-function(a,b,id,ex,color,xlab,ylab,backON="backOFF",ylimit=c(0,0),drawing=1,title=""){
  dataset=data.frame(a,b,id,ex) 
  dataset_omit<-as.data.frame(dataset[is.finite(rowSums(dataset)),])
  colnames(dataset_omit)=c("a","b","ID","Ex")
  Dev<-c(dataset_omit$a,dataset_omit$b)
  model_type<-c(rep("a",length(dataset_omit$a)), rep("b",length(dataset_omit$b)))
  emissionIndex<-c(as.factor(1:length(dataset_omit$a)), as.factor(1:length(dataset_omit$b)))
  ID<-c(dataset_omit$ID,dataset_omit$ID)
  Ex<-c(dataset_omit$Ex,dataset_omit$Ex)
  ER_data<-data.frame(Dev= Dev, model_type=model_type, emissionIndex=emissionIndex, ID=ID,Ex=Ex)
  
  (lmer(Dev ~ model_type + (1|emissionIndex) + (1|ID) +(1|Ex), data= ER_data)->model) 
  (lmer(Dev ~ 1 + (1|emissionIndex) + (1|ID)+(1|Ex), data= ER_data)->nullmodel) 
  
  
  print(summary(model))
  coeff_fix=as.data.frame(t(as.data.frame(fixef(model))))
  
  ratio=(coeff_fix[1]+coeff_fix[2])/coeff_fix[1]
  print("performance ratio:")
  print(ratio)
  
  lrt=anova(model,nullmodel,test = "Chisq")
  print("lrt test:")
  print(lrt)
  
  #manual calc. for wald test
  coeff=summary(model)$coefficients
  diff=coeff[2,1]/(coeff[1,1]+coeff[2,1])
  wald=1-pchisq((coeff[2,1]/coeff[2,2])^2,1) 
  print("wald test:")
  print(wald)
  
  plotdata<-data.frame(rep(0,nrow(dataset_omit)), 
                       dataset_omit$a,
                       rep(1,nrow(dataset_omit)),
                       dataset_omit$b)
  colnames(plotdata)=c("x1","a","x2","b")
  
  plotLMM<-data.frame(c(0,1),c(coeff[1,1],coeff[1,1]+coeff[2,1]))
  colnames(plotLMM)=c("x","y")
  
  draw<-ggplot()
  seg<-geom_segment(aes(x = plotdata[,1], y = plotdata[,2], xend = plotdata[,3], yend = plotdata[,4]), colour = rgb(0, 0, 0, alpha=0.15),data = plotdata)
  
  size=6
  point_a<-geom_point(data=plotLMM[,1:2], aes(x=plotLMM[,1:2][,1], y=plotLMM[,1:2][,2]),
                      shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = color)
  point_b<-geom_point(data=plotdata[,1:4], aes(x=plotdata[,3], y=plotdata[,4]),
                      shape = 21, size = size, stroke = 0.3, colour = rgb(0, 0, 0, alpha=0.6), fill = color)
  
  if(backON=="backOFF"){
    draw<-draw+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA), 
      axis.text= element_blank(),
      plot.title = element_blank(),
      axis.title = element_blank()
    )
  }
  
  
  base=coeff[1,1]+coeff[2,1]/2
  range=abs(coeff[2,1])*4
  
  ymin=base-range
  ymax=base+range
  
  if(ylimit[1]!=0||ylimit[2]!=0){
    ymin=ylimit[1]
    ymax=ylimit[2]
  }
  
  draw<-draw+seg
  draw<-draw+coord_cartesian(xlim=c(-1,2),ylim=c(ymin,ymax))
  draw<-draw+ylab(ylab)+ggtitle(title)
  lmmline<-geom_segment(aes(x = plotLMM[1,1], y = plotLMM[1,2], xend = plotLMM[2,1], yend = plotLMM[2,2]), colour = color,size=2.5, data=plotLMM)
  draw<-draw+point_a
  draw<-draw+lmmline
  if(drawing){
    print(draw)
  }
  
  if(title!=""){
    ggsave(title, draw, bg = "transparent",dpi=dpi2,width=11.3,height=10)
  }
  return(c(ratio[,1],lrt$`Pr(>Chisq)`[2],ymin,ymax,list(plotLMM)))
}

paired_comparison<-function(a,b,id,ex){
  dataset=data.frame(a,b,id,ex) 
  dataset_omit<-as.data.frame(dataset[is.finite(rowSums(dataset)),])
  colnames(dataset_omit)=c("a","b","ID","Ex")
  Dev<-c(dataset_omit$a,dataset_omit$b)
  model_type<-c(rep("a",length(dataset_omit$a)), rep("b",length(dataset_omit$b)))
  emissionIndex<-c(as.factor(1:length(dataset_omit$a)), as.factor(1:length(dataset_omit$b)))
  ID<-c(dataset_omit$ID,dataset_omit$ID)
  Ex<-c(dataset_omit$Ex,dataset_omit$Ex)
  ER_data<-data.frame(Dev= Dev, model_type=model_type, emissionIndex=emissionIndex, ID=ID,Ex=Ex)
  
  (lmer(Dev ~ model_type + (1|emissionIndex) + (1|ID) +(1|Ex), data= ER_data)->model) 
  (lmer(Dev ~ 1 + (1|emissionIndex) + (1|ID)+(1|Ex), data= ER_data)->nullmodel) 
  
  mod.type=summary(model)[15]
  a=fixed.effects(model)[2]
  b=fixed.effects(model)[1]
  lrt1=LRT(model,list(nullmodel))[,1]
  ratio=(b+a)/b
  
  result=data.frame(a,b,lrt1,ratio)
  colnames(result)=c("a","b","p(y=b)","ratio")
  rownames(result)=c(1:nrow(result))
  return(result)
}


biastest<-function(y,ID){
  d=data.frame(y,ID)
  d_omit<-as.data.frame(d[is.finite(rowSums(d)),])
  colnames(d_omit)=c("y","ID")
  lmer(y ~ 1+(1|ID),data=d_omit)->model1
  lmer(y ~ 0+(1|ID),data=d_omit)->nullmodel
  print(summary(model1))
  print(anova(model1, nullmodel))
}

Oltho_Projection<-function(from, to){
  lengthcheck=all(sapply(list(nrow(from),
                              nrow(to)), 
                         function(x) x == nrow(from)))
  if(!lengthcheck){
    print("length not matched")
    break()
  }
  result=data.frame(matrix(rep(NA,3), nrow=1))[numeric(0), ]
  colnames(result)=c("x","y","z")
  colnames(from)=c("x","y","z")
  colnames(to)=c("x","y","z")
  
  if(nrow(from)==0){
    return(result)
  }
  
  mag_from=sqrt(from$x^2+from$y^2+from$z^2)
  mag_to=sqrt(to$x^2+to$y^2+to$z^2)
  angle=angle_2vec_DF(from,to)
  result=mag_from*cos(angle*pi/180.0)*to/mag_to
  
  return(result)
}


multi_plot<-function(data,xlim,ylim,size=2,xlab,ylab,title,backOFF="backOFF",colset=NULL){
  
  
  
  linelist=list()
  for(i in 1:length(data)){
    if(is.null(colset)){
      d5=(i-1)%%4
      d4=(i-1)%%2
      d3=((i-1)%%4)%/%2
      d2=((i-1)%%8)%/%4
      d1=((i-1)%%16)%/%8
      col=hsv(d2*0.5,60/100,(d5+1)*(1/4), alpha=(2-d1)*0.5)
    }else{
      col=colset[[i]]
    }
    d0=((i-1)%%16)%/%15
    lsize=(d0+1)*2
    
    linelist=c(linelist,
               geom_line(data=data[[i]],
                         aes(x=x,y=y),
                         colour = col,
                         size=lsize))
  }
  
  
  
  p<-ggplot()
  p<-p+linelist
  
  # p<-p+text
  if(backOFF=="backOFF"){
    p<-p+theme(
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA),
      axis.text= element_blank(),
      plot.title = element_blank(),
      axis.title = element_blank())
  }
  
  p<-p+coord_cartesian(xlim=xlim,ylim=ylim)
  
  
  p<-p+xlab(xlab)+ylab(ylab)
  
  if(!is.null(title)){
    ggsave(paste0(F_dir,"/",title,".png"), p, bg = "transparent",dpi=dpi2,width=10*51/60,height=10)
  }
  return(p)
}

parameterlistNames<- function(...){
  nm1 <- as.character(eval(substitute(alist(...))))
  nm2=substr(nm1,regexpr("\\(", nm1)+1,regexpr("\\)", nm1)-1)
  nm3=strsplit(nm2,",")[[1]]
  return(nm3)
}

saveFigures<-function(...){
  
  nm1 <- as.character(eval(substitute(alist(...))))
  nm2=substr(nm1,regexpr("\\(", nm1)+1,regexpr("\\)", nm1)-1)
  names=strsplit(nm2,",")[[1]]
  
  for(i in 1:length(...)){
    ggsave(paste0(names[i],".png"), ...[[i]], bg = "transparent",dpi=dpi2,width=10,height=10)
  }
}

saveTables<-function(...){
  nm1 <- as.character(eval(substitute(alist(...))))
  nm2=substr(nm1,regexpr("\\(", nm1)+1,regexpr("\\)", nm1)-1)
  names=strsplit(nm2,",")[[1]]
  
  for(i in 1:length(...)){
    write.csv(...[[i]], file = paste0(names[i],".csv"))
  }
}



getFig_type1<-function(x,y,ID,Ex,
                       xlim,ylim,xlab,ylab,
                       col,back="backON", density=NULL,
                       title=NULL,F_dir=".",T_dir="."){
  mod_x=model_selection(NULL,x,ID,Ex,0)
  mod_y=model_selection(NULL,y,ID,Ex,0)
  
  fig=ggplot()
  fig=ggplot_test(fig,x,y,NULL,NULL,xlab,ylab,NULL,NULL,NA,0,back)
  
  if(!is.null(density)){
    if(density=="heatmap"){
      fig<-fig+
        stat_density2d(data=data.frame(x,y), aes(x=data.frame(x,y)[,1], y=data.frame(x,y)[,2], alpha=..density..), 
                       fill=rgb(1,0.4,0.4,alpha=0.8), geom="raster", contour = "false")
    }
    if(density=="contour"){
      fig<-fig+geom_density2d(data=data.frame(x,y), aes(x=data.frame(x,y)[,1:2][,1], y=data.frame(x,y)[,1:2][,2]),col = col,size=0.5)
    }
  }
  
  
  
  
  fig<-fig+scale_x_reverse(limits = c(xlim[2],xlim[1])) 
  fig<-fig+ylim(ylim)
  fig<-fig+coord_fixed()
  
  
  bin=(xlim[2]-xlim[1])/200
  fig<-ggMarginal(
    fig,
    binwidth = bin,
    type = "histogram",
    margins = "both",
    size = 7,
    col=hsv(0,0,0.5),
    fill=hsv(0,0,0.5)
  )
  
  
  if(!is.null(T_dir)){
    write.csv(mod_x[[2]],paste0(T_dir,"/",title,"_side.csv"))
    write.csv(mod_y[[2]],paste0(T_dir,"/",title,"_lead.csv"))
  }
  if(!is.null(F_dir)){
    ggsave(paste0(F_dir,"/",title,".png"), fig, bg = "transparent",dpi=dpi1,width=10,height=10)
  }
  return(fig)
}

getFig_type2<-function(x,y,ID,Ex,poly,
                       xlim,ylim,xlab,ylab,
                       col,refline=0,back="backON",
                       title=NULL,F_dir=".",T_dir="."){
  mod=model_selection(x,y,ID,Ex,poly)
  fig=ggplot()
  fig=ggplot_test(fig,x,y,xlim,ylim,xlab,ylab,title,list(mod),col,refline,back)
 
  if(!is.null(T_dir)){
    write.csv(mod[[2]],paste0(T_dir,"/",title,".csv"))
  }
  if(!is.null(F_dir)){
    ggsave(paste0(F_dir,"/",title,".png"), fig, bg = "transparent",dpi=dpi2,width=11.3,height=10)
  }
  return(fig)
}

getFig_type3<-function(x,y,ID,Ex,poly,
                       xlim,ylim,xlab,ylab,
                       col,refline=0,back="backON",
                       title=NULL,F_dir=".",T_dir=".",reverse="reverseON"){
  mod=model_selection_addNLS(x,y,ID,Ex,poly)
  
  fig=ggplot()
  fig=ggplot_test(fig,x,y,xlim,ylim,xlab,ylab,title,list(mod),col,refline,back)
  
  if(reverse=="reverseON"){
    fig<-fig+scale_x_reverse()
  }
  if(!is.null(T_dir)){
    write.csv(mod[[2]],paste0(T_dir,"/",title,".csv"))
  }
  if(!is.null(F_dir)){
    ggsave(paste0(F_dir,"/",title,".png"), fig, bg = "transparent",dpi=dpi2,width=11.3,height=10)
  }
  return(fig)
}


##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################



#simple controller
generate_results=T
draw_EachSession=T
draw_RGLfigures=T
generate_results_delay=F ##it takes long time (several hours)

dpi1=100 ##600 for print, 100 for quick view
dpi2=100 
feedbackdelay_Flight=0.120
feedbackdelay_Sound=0.120
pi=atan2(1,1)*4
G=9.80665
fps=125
dist_threshold=0.3
par(mfrow=c(1,1))
par(mar=c(3.2, 3.2, 3.2, 3.2), mgp=c(2.0, 0.7, 0))
options(warn = -1)
old.op <- options(max.print=10000)
background="backON"


initialdir=here()
setwd(initialdir)
print(initialdir)
setwd("..")
basedir=getwd()
expdir=paste0(basedir,"/experiment")
IDdir=paste0(basedir,"/basic_info")

setwd(IDdir)
Exp_Bat=read.csv("ID_table.csv",header=T)
setwd(expdir)
files  <- list.files()



AllSessions=data.frame()
AllFlight=data.frame()
AllFlight_dir=data.frame()
ExpIndex=0

for (file.name in files) {
  
  head3=substr(file.name, 1, 3)
  head1=substr(head3, 1, 1)
  if (regexpr("E", head1)  < 0) {next}
  subdir=paste0(expdir,"/",file.name)
  setwd(subdir)
  
  position=read.csv(paste0(head3,"_position.csv"), header=T)
  pulse=read.csv(paste0(head3,"_pulse.csv"), header=T)
  
  
  format=initialformatting(position,pulse)
  position=format[[1]]
  pulse=format[[2]]
  
  #########making motion functions
  motionFunc=list()
  for(funcindex in 1:6){
    func_each <- smooth.spline(position$time, position[,funcindex+1])
    motionFunc=c(motionFunc,list(func_each))
  }
  pulse=addKinematicData(pulse,motionFunc)
  pulse=cbind(pulse,forNULL=c(1:nrow(pulse)))
  ####omit low precision data due to short distance
  select=which(pulse$R_mag>dist_threshold)
  if(length(select)==0){next}
  
  ### fragmentation of data to make non-continuous sections clear 
  pulse_frag=fragmentation(pulse[select,])[1]
  
  res=data.frame()
  for(j in 1:length(pulse_frag)){
    pulse_frag[[j]]=getLastDetectedTimestep(pulse_frag[[j]],feedbackdelay_Sound)
  }
  
  ## adding parameters in each fragment to pulse_frag then storing to EachSession 
  EachSight=SightAnalysis(pulse_frag)
  EachFlight=FlightAnalysis(pulse_frag)
  
  EachSession=cbind(EachSight,EachFlight)
  
  
  
  ## adding index
  ID=as.numeric(Exp_Bat[Exp_Bat$Exp.ID==head3,]$Bat.ID)
  ExpIndex=ExpIndex+1
  EachSession=cbind(ID,ExpID=head3,ExpIndex,EachSession)
  
  ## storing each session data (EachSession) to AllSessions
  AllSessions=rbind(AllSessions,EachSession)
  
  if(draw_EachSession){
    s_Position=smoothed_position(position,motionFunc,60*5)
    R_=data.frame(x=(s_Position$x_mx-s_Position$x_bx),
                  y=(s_Position$x_my-s_Position$x_by),
                  z=(s_Position$x_mz-s_Position$x_bz))
    R_mag_=sqrt(R_$x^2+R_$y^2+R_$z^2)
    select_=which(s_Position$time>=min(pulse_frag[[1]]$time)&s_Position$time<=max(pulse_frag[[1]]$time))
    
    s_bat=s_Position[select_,2:4]
    s_moth=s_Position[select_,5:7]
    s_R=direction2(s_moth,s_bat)$hv
    pd=data.frame(EachSight$P.H,EachSight$P.V)
    td=data.frame(EachSight$LOSc.H,EachSight$LOSc.V)
    pd_v=data.frame(EachSight$`P(-1).H`,EachSight$`P(-1).V`)
    td_v=data.frame(EachSight$`LOSc(-1).H`,EachSight$`LOSc(-1).V`)
    pd_r2=data.frame(EachSight$LOSl_P.H,EachSight$LOSl_P.V)
    td_r2=data.frame(EachSight$LOSl_LOSc.H,EachSight$LOSl_LOSc.V)
    
    drawDirection3d(NULL,pd_v,td_v,"Horizontal Angle","Vertical Angle","default_coordinates(before rotation)")
    drawDirection3d(NULL,pd_r2,td_r2,"Horizontal Angle","Vertical Angle","defined_coordinates(after rotation)")
    Polar_line3(s_R,
                pd,
                td,
                8,0,
                "Horizontal Angle","Vertical Angle",
                "pointON_b","pathON","densityOFF","textOFF","backON",paste("graphics/pulse_direction/","default_coordinates_full",".png", sep=""),0)
    
    drawScene(pulse[select,],s_Position[select_,])
    
    range_ori=continuous(!is.na(EachSession$TGTdir.ori_p.mn.Rmag))
    range_pure=continuous(!is.na(EachSession$TGTdir.pure_p.mn.Rmag))
    range=which(range_ori&range_pure)
    t_range=c(EachSession$time[min(range)],EachSession$time[max(range)])
    
    select_=which(s_Position$time>=t_range[1]&s_Position$time<=t_range[2]+feedbackdelay_Flight)
    olthoPlot_sim(data.frame(s_Position[select_,]$x_bx,s_Position[select_,]$x_by),
                  data.frame(s_Position[select_,]$x_mx,s_Position[select_,]$x_my),
                  data.frame(EachSession[range,]$TGTdir.pure_p.mn.x,EachSession[range,]$TGTdir.pure_p.mn.y),
                  
                  data.frame(EachSession[range,]$x_bx_L,EachSession[range,]$x_by_L),
                  data.frame(EachSession[range,]$x_mx_L,EachSession[range,]$x_my_L),
                  data.frame(EachSession[range,]$TGTdir.pure_p.mn.x,EachSession[range,]$TGTdir.pure_p.mn.y),
                  
                  2,"x","y","textON","backON","graphics/scene/simscene2Dx_y_dir_pure.png")
    
    range_ori=continuous(!is.na(EachSession$TGTloc.ori_p.mn.Rmag))
    range_pure=continuous(!is.na(EachSession$TGTloc.pure_p.mn.Rmag))
    range=which(range_ori&range_pure)
    t_range=c(EachSession$time[min(range)],EachSession$time[max(range)])
    
    select_=which(s_Position$time>=t_range[1]&s_Position$time<=t_range[2]+feedbackdelay_Flight)
    olthoPlot_sim(data.frame(s_Position[select_,]$x_bx,s_Position[select_,]$x_by),
                  data.frame(s_Position[select_,]$x_mx,s_Position[select_,]$x_my),
                  data.frame(EachSession[range,]$TGTloc.pure_p.mn.x,EachSession[range,]$TGTloc.pure_p.mn.y),
                  data.frame(EachSession[range,]$x_bx_L,EachSession[range,]$x_by_L),
                  data.frame(EachSession[range,]$x_mx_L,EachSession[range,]$x_my_L),
                  data.frame(EachSession[range,]$TGTloc.pure_p.mn.x,EachSession[range,]$TGTloc.pure_p.mn.y),
                  
                  2,"x","y","textON","backON","graphics/scene/simscene2Dx_y_loc.png")
  }
}



if(generate_results){
  setwd(basedir)
  getwd()
  if(!dir.exists("results")){
    dir.create("results")
  }
  setwd(paste0(basedir,"/results"))
  
  dirs=c("Figures","Supp data")
  creat_dir<-function(name){
    if(!dir.exists(name)){
      dir.create(name)
    }
  }
  lapply(dirs, creat_dir)
  
  setwd("./Figures")
  dirs=c("Fig1","Fig2","Fig3","Fig4","Fig5","FigS1","FigS2","FigS3","FigS4")
  creat_dir<-function(name){
    if(!dir.exists(name)){
      dir.create(name)
    }
  }
  lapply(dirs, creat_dir)
  
  setwd("../Supp data")
  dirs=c("DataS1","DataS2","DataS3","DataS4","DataS5")
  creat_dir<-function(name){
    if(!dir.exists(name)){
      dir.create(name)
    }
  }
  lapply(dirs, creat_dir)
  
  ##load pre_generated results of delay effect (see line XX~ for detailed generating procedure)
  load(paste0(basedir,"/results/sim of delay effect/delay_result_70_170ms.Rdata"))
  colnames(delay_result)
  x_lim=c(0.07,0.17)
  head(delay_result)
  
  for(i in 1:ncol(delay_result)){
    delay_result[,i]=as.numeric(delay_result[,i])
  }
  
  F_dir=paste0(basedir,"/results/Figures/Fig2")
  T_dir=paste0(basedir,"/results/Supp data/DataS1")
  
  generate_Fig2<-function(){
    #Fig.1B
    Fig2B=getFig_type1(AllSessions$LOSl_P_scaled.SideAngle,AllSessions$LOSl_P_scaled.LeadAngle,AllSessions$ID,AllSessions$ExpIndex,
                       c(-3,3),c(-3,3),"Normalized pulse side angle","Normalized pulse lead angle",green,background,"contour",
                       "Fig2b",F_dir,NULL)

    Fig2C=multi_paired_comp(list(AllSessions$LOSl_P.V,
                                 AllSessions$LOSc_P.V),
                            AllSessions$ID,
                            AllSessions$ExpIndex,
                            c(0,40),
                            "Target Direction","Deviation angle",
                            c(green),background,
                            "Fig2c",F_dir,T_dir)
    
    Fig2D=getFig_type2(AllSessions$LOSl_LOSc.LeadAngle,AllSessions$LOSl_P.LeadAngle,AllSessions$ID,AllSessions$ExpIndex,
                       1,c(-2,32),c(-40,40),"Relative motion angle","Pulse lead angle",green,2,background,
                       "Fig2d",F_dir,T_dir)
    
    Fig2E=getFig_type2(AllSessions$LOSl_LOSc.LeadAngle,AllSessions$LOSl_P.SideAngle,AllSessions$ID,AllSessions$ExpIndex,
                       1,c(-2,32),c(-40,40),"Relative motion angle","Pulse side angle",green,1,background,
                       "Fig2e",F_dir,T_dir)
    
    Fig2F=getFig_type3(AllSessions$R_mag,1/AllSessions$interval,AllSessions$ID,AllSessions$ExpIndex,
                       1,c(0,4),c(0,80),"Distance","Sensing rate",green,0,background,
                       "Fig2f",F_dir,T_dir,"reverseON")
    
    Fig2G=getFig_type3(AllSessions$R_mag,AllSessions$widthHV,AllSessions$ID,AllSessions$ExpIndex,
                       1,c(0,4),c(25,125),"Distance","Beamwidth",green,0,background,
                       "Fig2g",F_dir,T_dir,"reverseON")
    
    Fig2J=getFig_type1(AllSessions$LOSl_F_V.SideAngle,AllSessions$LOSl_F_V.LeadAngle,AllSessions$ID,AllSessions$ExpIndex,
                       c(-75,75),c(-75,75),"Flight side angle","Flight lead angle",blue,background,NULL,
                       "Fig2j",F_dir,T_dir)
    
  }
  generate_Fig2()
  
  F_dir=paste0(basedir,"/results/Figures/Fig3")
  T_dir=paste0(basedir,"/results/Supp data/DataS2")
  
  gerenate_Fig3<-function(){
    cdr=getwd()  
    
    s1=data.frame(av=AllSessions$sim000.av,se=AllSessions$sim000.S.Error,abs_se=abs(AllSessions$sim000.S.Error),simtype=as.factor(rep(1,nrow(AllSessions))),id=AllSessions$ID,exp=AllSessions$ExpIndex)
    s2=data.frame(av=AllSessions$sim001.av,se=AllSessions$sim001.S.Error,abs_se=abs(AllSessions$sim001.S.Error),simtype=as.factor(rep(2,nrow(AllSessions))),id=AllSessions$ID,exp=AllSessions$ExpIndex)
    s3=data.frame(av=AllSessions$sim010.av,se=AllSessions$sim010.S.Error,abs_se=abs(AllSessions$sim010.S.Error),simtype=as.factor(rep(3,nrow(AllSessions))),id=AllSessions$ID,exp=AllSessions$ExpIndex)
    s4=data.frame(av=AllSessions$sim011.av,se=AllSessions$sim011.S.Error,abs_se=abs(AllSessions$sim011.S.Error),simtype=as.factor(rep(4,nrow(AllSessions))),id=AllSessions$ID,exp=AllSessions$ExpIndex)
    s5=data.frame(av=AllSessions$sim100.av,se=AllSessions$sim100.S.Error,abs_se=abs(AllSessions$sim100.S.Error),simtype=as.factor(rep(5,nrow(AllSessions))),id=AllSessions$ID,exp=AllSessions$ExpIndex)
    s6=data.frame(av=AllSessions$sim101.av,se=AllSessions$sim101.S.Error,abs_se=abs(AllSessions$sim101.S.Error),simtype=as.factor(rep(6,nrow(AllSessions))),id=AllSessions$ID,exp=AllSessions$ExpIndex)
    s7=data.frame(av=AllSessions$sim110.av,se=AllSessions$sim110.S.Error,abs_se=abs(AllSessions$sim110.S.Error),simtype=as.factor(rep(7,nrow(AllSessions))),id=AllSessions$ID,exp=AllSessions$ExpIndex)
    s8=data.frame(av=AllSessions$sim111.av,se=AllSessions$sim111.S.Error,abs_se=abs(AllSessions$sim111.S.Error),simtype=as.factor(rep(8,nrow(AllSessions))),id=AllSessions$ID,exp=AllSessions$ExpIndex)
    
    al=0.8
    multi_plot_modeling2(list(s1,s2,s3,s4,s6,s5,s7,s8),
                         AllSessions$ID,AllSessions$ExpIndex,
                         c(0,200),c(-2,0.5),
                         "Angular velocity","Sighting error", 
                         2,
                         list(hsv(0,0,0.1,alpha=al),hsv(0,0,0.2,alpha=al),hsv(0,0,0.4,alpha=al),hsv(0,0,0.5,alpha=al),
                              hsv(0.45,0.7,0.5,alpha=al),hsv(0.45,0.7,0.3,alpha=al),hsv(0.45,0.8,0.6,alpha=al),hsv(0.45,0.9,0.9)),
                         background,
                         "Fig3b",F_dir,T_dir) 
    
    
    setwd(paste0(basedir,"/experiment/E09/graphics/"))
    file.copy("scene/simscene2Dx_y_dir_pure.png",paste0(F_dir,"/Fig3d.png"),overwrite = T)
    setwd(cdr)
    
    Fig3E=multi_paired_comp(list(AllSessions$TGTdir.pure_p.sn.Wmag,
                                 AllSessions$TGTdir.ori_p.sn.Wmag),
                            AllSessions$ID,
                            AllSessions$ExpIndex,
                            c(0,200),
                            "type of flight pattern","Angular velocity",
                            c(blue),background,
                            "Fig3e",F_dir,T_dir)
    
    df_sm0000=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim0000)
    df_sm0001=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim0001)
    df_sm0010=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim0010)
    df_sm0011=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim0011)
    df_sm0100=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim0100)
    df_sm0101=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim0101)
    df_sm0110=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim0110)
    df_sm0111=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim0111)
    df_sm1000=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim1000)
    df_sm1001=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim1001)
    df_sm1010=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim1010)
    df_sm1011=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim1011)
    df_sm1100=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim1100)
    df_sm1101=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim1101)
    df_sm1110=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim1110)
    df_sm1111=data.frame(x=delay_result$feedbackdelay_Flight,y=delay_result$sim1111)
    
    
    
    al=0.6
    multi_plot(list(df_sm0000,df_sm0010,df_sm0100,df_sm0110,
                    df_sm0001,df_sm0101,df_sm0011,df_sm0111,
                    df_sm1000,df_sm1010,df_sm1100,df_sm1110,
                    df_sm1001,df_sm1011,df_sm1101,
                    
                    df_sm1111),
               x_lim,c(0,2.0),2,
               "Delay","Sighting error (absolute value)","Fig3f",
               background,
               list(hsv(0,0,0.1,alpha=al),hsv(0,0,0.2,alpha=al),hsv(0,0,0.35,alpha=al),hsv(0,0,0.4,alpha=al),
                    hsv(0.58,0.8,0.9),hsv(0.6,0.8,0.9,alpha=al),hsv(0.6,0.6,0.8,alpha=al),hsv(0.6,0.5,0.7,alpha=al),
                    hsv(0.45,0.9,0.8),hsv(0.45,0.7,0.7,alpha=al),hsv(0.45,0.7,0.7,alpha=al),hsv(0.45,0.3,0.7,alpha=al),
                    hsv(0.93,0.7,0.8, alpha = al),hsv(0.93,0.7,0.7,alpha=al),hsv(0.93,0.7,0.7,alpha=al),hsv(0.93,0.7,1))
    )
    
  }
  gerenate_Fig3()
  
  F_dir=paste0(basedir,"/results/Figures/Fig4")
  T_dir=paste0(basedir,"/results/Supp data/DataS3")
  
  gerenate_Fig4<-function(){
    
    ###Fig.3
    Fig4B_left=multi_paired_comp(list(AllSessions$P_L.V,AllSessions$P_A.V),
                                 AllSessions$ID,
                                 AllSessions$ExpIndex,
                                 c(0,30),
                                 "type of prediction model",
                                 "Fitting error",
                                 c(green),background,
                                 "Fig4b",F_dir,T_dir)
    
    
    Fig4B_right=Overlayed_paired_comp(
      list(AllSessions$P_L.V,AllSessions$P_A.V),
      list(AllSessions$LOSc_L.V,AllSessions$LOSc_A.V),
      AllSessions$ID,
      AllSessions$ExpIndex,
      c(0,30),
      "Type of prediction",
      "Fitting error to pulse",
      "Fitting error to target",
      c(green),c(red),c(green2),c(red2),
      background,
      "Fig4b",F_dir,T_dir)
    
    Fig4C=getFig_type3(AllSessions$perceived_av_Sound*180/pi,1/AllSessions$interval,AllSessions$ID,AllSessions$ExpIndex,
                       1,c(-22,352),c(0,80),"Angular velocity","Scan Rate",green,0,background,
                       "Fig4c",F_dir,T_dir,"reverseOFF")
    
    Fig4D=getFig_type3(AllSessions$perceived_av_Sound*180/pi,AllSessions$widthHV,AllSessions$ID,AllSessions$ExpIndex,
                       1,c(-22,352),c(25,120),"Angular velocity","Beamwidth",green,0,background,
                       "Fig4d",F_dir,T_dir,"reverseOFF")
    
    Fig4E=getFig_type3(AllSessions$perceived_av_Flight*180/pi,AllSessions$LOSl_F_V.LeadAngle,AllSessions$ID,AllSessions$ExpIndex,
                       1,c(-22,352),c(-50,100),"Angular velocity","Flight lead angle",blue,0,background,
                       "Fig4e",F_dir,T_dir,"reverseOFF")
    
    Fig4F=getFig_type3(AllSessions$perceived_av_Flight*180/pi,AllSessions$LOSl_F_V.SideAngle,AllSessions$ID,AllSessions$ExpIndex,
                       1,c(-22,352),c(-80,80), "Angular velocity","Flight side angle",blue,0,background,
                       "Fig4f",F_dir,T_dir,"reverseOFF")
    
  }
  gerenate_Fig4()
  
  
  F_dir=paste0(basedir,"/results/Figures/FigS1")
  
  gerenate_figS1<-function(){
    from_dir=paste0(basedir,"/experiment/E09/graphics")
    file.copy(paste0(from_dir,"/scene/scene2Dx_y.png"),paste0(F_dir,"/figS1a_upper.png"),overwrite = T)
    file.copy(paste0(from_dir,"/scene/scene2Dx_z.png"),paste0(F_dir,"/figS1a_lower.png"),overwrite = T)
    file.copy(paste0(from_dir,"/pulse_direction/default_coordinates_full.png"),paste0(F_dir,"/figS1c_left.png"),overwrite = T)
    file.copy(paste0(from_dir,"/pulse_direction/defined_coordinates(after rotation).png"),paste0(F_dir,"/figS1c_right.png"),overwrite = T)
  }
  gerenate_figS1()
  
  F_dir=paste0(basedir,"/results/Figures/FigS2")
  T_dir=paste0(basedir,"/results/Supp data/DataS4")
  
  gerenate_figS2<-function(){
    
    FigS3=getFig_type3(AllSessions$R_mag,
                       AllSessions$av_mag,
                       AllSessions$ID,AllSessions$ExpIndex,
                       1,c(-5/16,5),c(0,500),"Distance","Angular velocity",green,0,background,
                       "figS2b",F_dir,T_dir,"reverseOFF")
  }
  gerenate_figS2()
  
  
  F_dir=paste0(basedir,"/results/Figures/FigS3")
  setwd(F_dir)
  
  
  ggmultiplotter(delay_result$feedbackdelay_Sound,list(delay_result$f2C_b, delay_result$f2C_b+delay_result$f2C_a), 
                 x_lim,c(-2,20),
                 "Delay","Fiting error in fig2c",
                 "figS3a_upper",background,1,1)
  
  ggplotter(delay_result$feedbackdelay_Sound,delay_result$f2C_p, 
            x_lim,c(0,0.07),
            "Delay","P value",
            "figS3a_lower",background,1,1)
  
  
  ggplotter(delay_result$feedbackdelay_Sound,delay_result$f2D_slope, 
            x_lim,c(0,2),
            "Delay","Slope of f1d",
            "figS3b_upper",background,1,1)
  ggplotter(delay_result$feedbackdelay_Sound,delay_result$f2D_p, 
            x_lim,c(0,0.07),
            "Delay","P value",
            "figS3b_lower",background,1,1)
  
  df_sm000=data.frame(x=delay_result$feedbackdelay_Sound,y=delay_result$sim000)
  df_sm001=data.frame(x=delay_result$feedbackdelay_Sound,y=delay_result$sim001)
  df_sm010=data.frame(x=delay_result$feedbackdelay_Sound,y=delay_result$sim010)
  df_sm011=data.frame(x=delay_result$feedbackdelay_Sound,y=delay_result$sim011)
  df_sm100=data.frame(x=delay_result$feedbackdelay_Sound,y=delay_result$sim100)
  df_sm101=data.frame(x=delay_result$feedbackdelay_Sound,y=delay_result$sim101)
  df_sm110=data.frame(x=delay_result$feedbackdelay_Sound,y=delay_result$sim110)
  df_sm111=data.frame(x=delay_result$feedbackdelay_Sound,y=delay_result$sim111)
  al=0.8
  
  multi_plot(list(df_sm000,df_sm001,df_sm010,df_sm011,df_sm100,df_sm101,df_sm110,df_sm111),
             x_lim,c(0,1.5),2,
             "Delay","Sighting error (absolute value)","figS3c",
             background,
             list(hsv(0,0,0.1,alpha=al),hsv(0,0,0.2,alpha=al),hsv(0,0,0.4,alpha=al),hsv(0,0,0.5,alpha=al),
                  hsv(0.45,0.7,0.3,alpha=al),hsv(0.45,0.7,0.4,alpha=al),hsv(0.45,0.8,0.6,alpha=al),hsv(0.45,0.9,0.9))
  )
  
  
  ggmultiplotter(delay_result$feedbackdelay_Sound,list(delay_result$f4B_g_b, delay_result$f4B_g_b+delay_result$f4B_g_a), 
                 x_lim,c(-2,26),
                 "Delay","Fitting Error in f4B_green",
                 "figS3d_upper",background,1,1)
  
  ggplotter(delay_result$feedbackdelay_Sound,delay_result$f4B_g_p, 
            x_lim,c(0,0.07),
            "Delay","P value",
            "figS3d_lower",background,1,1)
  
  ggmultiplotter(delay_result$feedbackdelay_Sound,list(delay_result$f4B_r_b, delay_result$f4B_r_b+delay_result$f4B_r_a), 
                 x_lim,c(-2,26),
                 "Delay","Fitting Error in f4B_red",
                 "figS3d2_upper",background,1,1)
  
  ggplotter(delay_result$feedbackdelay_Sound,delay_result$f4B_r_p, 
            x_lim,c(0,0.07),
            "Delay","P value",
            "figS3d2_lower",background,1,1)
  
  ggplotter(delay_result$feedbackdelay_Sound,delay_result$f4C_v, 
            x_lim,c(-0.05,0.15),
            "Delay","Slope of fig4C",
            "figS3e_upper",background,1,1)
  ggplotter(delay_result$feedbackdelay_Sound,delay_result$f4C_p, 
            x_lim,c(0,0.07),
            "Delay","P value",
            "figS3e_lower",background,1,1)
  
  ggplotter(delay_result$feedbackdelay_Sound,delay_result$f4D_v, 
            x_lim,c(0.0,0.2),
            "Delay","Slope of fig4D",
            "figS3f_upper",background,1,1)
  ggplotter(delay_result$feedbackdelay_Sound,delay_result$f4D_p, 
            x_lim,c(0,0.07),
            "Delay","P value",
            "figS3f_lower",background,1,1)
  
  x_lim=c(0.07,0.17)
  ggplotter(delay_result$feedbackdelay_Flight,delay_result$f2J_lead_v, 
            x_lim,c(0,20),
            "Delay","Flight lead angle",
            "figS3g_upper",background,1,1)
  ggplotter(delay_result$feedbackdelay_Flight,delay_result$f2J_lead_p, 
            x_lim,c(0,0.07),
            "Delay","P value",
            "figS3g_lower",background,1,1)
  
  
  ggmultiplotter(delay_result$feedbackdelay_Flight,list(delay_result$f3E_b, delay_result$f3E_b+delay_result$f3E_a), 
                 x_lim,c(-20,200),
                 "Delay","A.V. in fig3E",
                 "figS3h2_upper",background,1,1)
  
  ggplotter(delay_result$feedbackdelay_Flight,delay_result$f3E_p, 
            x_lim,c(0,0.07),
            "Delay","P value",
            "figS3h_lower",background,1,1)
  
  ggplotter(delay_result$feedbackdelay_Flight,delay_result$f4E_slope, 
            x_lim,c(0,0.2),
            "Delay","Slope of fig4E",
            "figS3i_upper",background,1,1)
  ggplotter(delay_result$feedbackdelay_Flight,delay_result$f4E_p, 
            x_lim,c(0,0.07),
            "Delay","P value",
            "figS3i_lower",background,1,1)
  
  
  
  
  F_dir=paste0(basedir,"/results/Figures/FigS4")
  T_dir=paste0(basedir,"/results/Supp data/DataS5")
  gerenate_figS4<-function(){
    
    #TableS7
    figS5_2B=getFig_type1(AllSessions$Tl_P_scaled.SideAngle,AllSessions$Tl_P_scaled.LeadAngle,AllSessions$ID,AllSessions$ExpIndex,
                          c(-3,3),c(-3,3),"Normalized pulse side angle","Normalized pulse lead angle",green,background,"contour",
                          "reanalysis_2b",NULL,NULL)
    
    figS5_2C=multi_paired_comp(list(AllSessions$Tl_P.V,
                                    AllSessions$LOSc_P.V),
                               AllSessions$ID,
                               AllSessions$ExpIndex,
                               c(0,40),
                               "Target Direction","Deviation angle",
                               c(green),background,
                               "reanalysis_2c",NULL,T_dir)
    
    figS5_2D=getFig_type2(AllSessions$Tl_Tc.LeadAngle,AllSessions$Tl_P.LeadAngle,AllSessions$ID,AllSessions$ExpIndex,
                          1,c(-2,32),c(-40,40),"Relative motion angle","Pulse lead angle",green,2,background,
                          "reanalysis_2d",NULL,T_dir)
    
    figS5_2E=getFig_type2(AllSessions$Tl_Tc.LeadAngle,AllSessions$Tl_P.SideAngle,AllSessions$ID,AllSessions$ExpIndex,
                          1,c(-2,32),c(-40,40),"Relative motion angle","Pulse side angle",green,1,background,
                          "reanalysis_2e",NULL,T_dir)
    
    figS5_2J=getFig_type1(AllSessions$Tl_F_V.SideAngle,AllSessions$Tl_F_V.LeadAngle,AllSessions$ID,AllSessions$ExpIndex,
                          c(-75,75),c(-75,75),"Flight side angle","Flight lead angle",NULL,background,NULL,
                          "reanalysis_2j",NULL,T_dir)
    
    
    sl1=cbind(AllSessions$simloc000.av,AllSessions$simloc000.S.Error,id=AllSessions$ID,exp=AllSessions$ExpIndex)
    sl2=cbind(AllSessions$simloc001.av,AllSessions$simloc001.S.Error,id=AllSessions$ID,exp=AllSessions$ExpIndex)
    sl3=cbind(AllSessions$simloc010.av,AllSessions$simloc010.S.Error,id=AllSessions$ID,exp=AllSessions$ExpIndex)
    sl4=cbind(AllSessions$simloc011.av,AllSessions$simloc011.S.Error,id=AllSessions$ID,exp=AllSessions$ExpIndex)
    sl5=cbind(AllSessions$simloc100.av,AllSessions$simloc100.S.Error,id=AllSessions$ID,exp=AllSessions$ExpIndex)
    sl6=cbind(AllSessions$simloc101.av,AllSessions$simloc101.S.Error,id=AllSessions$ID,exp=AllSessions$ExpIndex)
    sl7=cbind(AllSessions$simloc110.av,AllSessions$simloc110.S.Error,id=AllSessions$ID,exp=AllSessions$ExpIndex)
    sl8=cbind(AllSessions$simloc111.av,AllSessions$simloc111.S.Error,id=AllSessions$ID,exp=AllSessions$ExpIndex)
    
    figS5_3B=multi_plot_modeling2(list(sl1,sl2,sl3,sl4,sl5,sl6,sl7,sl8),
                                  AllSessions$ID,AllSessions$ExpIndex,
                                  c(0,200),c(-2.0,0.5),
                                  "Angular velocity","Sighting error", 
                                  2,NULL,background,
                                  "reanalysis_3b",NULL,T_dir) 
    
    
    figS5_3E=multi_paired_comp(list(AllSessions$TGTloc.pure_p.sn.Wmag,
                                    AllSessions$TGTloc.ori_p.sn.Wmag),
                               AllSessions$ID,
                               AllSessions$ExpIndex,
                               c(0,200),
                               "type of flight pattern","Angular velocity",
                               c(blue),background,
                               "reanalysis_3e",NULL,T_dir)
    
    
    
    figS5_4E=getFig_type3(AllSessions$perceived_av_Flight*180/pi,AllSessions$Tl_F_V.LeadAngle,
                          AllSessions$ID,AllSessions$ExpIndex,
                          1,c(-22,352),c(-50,100),"Angular velocity","Flight lead angle",blue,0,background,
                          "reanalysis_4e",NULL,T_dir,"reverseOFF")
    
    figS5_4F=getFig_type3(AllSessions$perceived_av_Flight*180/pi,AllSessions$Tl_F_V.SideAngle,AllSessions$ID,AllSessions$ExpIndex,
                          1,c(-22,352),c(-80,80), "Angular velocity","Flight side angle",blue,0,background,
                          "reanalysis_4f",NULL,T_dir,"reverseOFF")
  }
  gerenate_figS4()
  setwd(initialdir)
  
}
setwd(initialdir)





if(!generate_results_delay){
  setwd(initialdir)
  invokeRestart("abort")
}


delay_result=data.frame()

for(k in 1:111){
  feedbackdelay_Flight=0.065+(k-1)*0.001
  feedbackdelay_Sound=0.065+(k-1)*0.001
  
  AllSessions=data.frame()
  AllFlight=data.frame()
  AllFlight_dir=data.frame()
  ExpIndex=0
  
  for (file.name in files) {
    
    head3=substr(file.name, 1, 3)
    head1=substr(head3, 1, 1)
    if (regexpr("E", head1)  < 0) {next}
    subdir=paste0(expdir,"/",file.name)
    setwd(subdir)
    
    position=read.csv(paste0(head3,"_position.csv"), header=T)
    pulse=read.csv(paste0(head3,"_pulse.csv"), header=T)
    
    
    format=initialformatting(position,pulse)
    position=format[[1]]
    pulse=format[[2]]
    
    #########making motion functions
    motionFunc=list()
    for(funcindex in 1:6){
      func_each <- smooth.spline(position$time, position[,funcindex+1])
      motionFunc=c(motionFunc,list(func_each))
    }
    pulse=addKinematicData(pulse,motionFunc)
    pulse=cbind(pulse,forNULL=c(1:nrow(pulse)))
    ####omit low precision data due to short distance
    select=which(pulse$R_mag>dist_threshold)
    if(length(select)==0){next}
    
    ### fragmentation of data to make non-continuous sections clear 
    pulse_frag=fragmentation(pulse[select,])[1]
    
    res=data.frame()
    for(j in 1:length(pulse_frag)){
      pulse_frag[[j]]=getLastDetectedTimestep(pulse_frag[[j]],feedbackdelay_Sound)
    }
    
    ## adding parameters in each fragment to pulse_frag then storing to EachSession 
    EachSight=SightAnalysis(pulse_frag)
    EachFlight=FlightAnalysis(pulse_frag)
    EachSession=cbind(EachSight,EachFlight)
    
    
    ## adding index
    ID=as.numeric(Exp_Bat[Exp_Bat$Exp.ID==head3,]$Bat.ID)
    ExpIndex=ExpIndex+1
    EachSession=cbind(ID,ExpID=head3,ExpIndex,EachSession)
    
    ## storing each session data (EachSession) to AllSessions
    AllSessions=rbind(AllSessions,EachSession)
  }
  
  
  f2C=paired_comparison(AllSessions$LOSl_P.V,
                        AllSessions$LOSc_P.V,
                        AllSessions$ID,
                        AllSessions$ExpIndex)
  f2C_a=f2C[1,1]
  f2C_b=f2C[1,2]
  f2C_v=f2C[1,4]
  f2C_p=f2C[1,3]
  
  f2D=model_selection(AllSessions$LOSl_LOSc.LeadAngle,AllSessions$LOSl_P.LeadAngle,
                      AllSessions$ID,AllSessions$ExpIndex,
                      1)
  f2D_slope=fixed.effects(f2D[[1]])[2]
  f2D_intercept=fixed.effects(f2D[[1]])[1]
  f2D_p=f2D[[2]][2,6]
  
  
  f2J_lead=model_selection(AllSessions$forNULL,AllSessions$LOSl_F_V.LeadAngle,
                           AllSessions$ID,AllSessions$ExpIndex,
                           0)
  f2J_lead_v=fixed.effects(f2J_lead[[1]])[1]
  f2J_lead_p=f2J_lead[[2]][1,5]
  
  f2J_side=model_selection(AllSessions$forNULL,AllSessions$LOSl_F_V.SideAngle,
                           AllSessions$ID,AllSessions$ExpIndex,
                           0)
  f2J_side_v=fixed.effects(f2J_side[[1]])[1]
  f2J_side_p=f2J_side[[2]][1,5]
  
  
  s1=data.frame(av=AllSessions$sim000.av,se=AllSessions$sim000.S.Error)
  s2=data.frame(av=AllSessions$sim001.av,se=AllSessions$sim001.S.Error)
  s3=data.frame(av=AllSessions$sim010.av,se=AllSessions$sim010.S.Error)
  s4=data.frame(av=AllSessions$sim011.av,se=AllSessions$sim011.S.Error)
  s5=data.frame(av=AllSessions$sim100.av,se=AllSessions$sim100.S.Error)
  s6=data.frame(av=AllSessions$sim101.av,se=AllSessions$sim101.S.Error)
  s7=data.frame(av=AllSessions$sim110.av,se=AllSessions$sim110.S.Error)
  s8=data.frame(av=AllSessions$sim111.av,se=AllSessions$sim111.S.Error)
  
  av_pure=AllSessions$TGTdir.pure_p.sn.Wmag
  av_ori=AllSessions$TGTdir.ori_p.sn.Wmag
  AllSessions=gather_pred_res(list(s1,s2,s3,s4,s5,s6,s7,s8),AllSessions$ID,AllSessions$ExpIndex,av_pure,AllSessions, "S.Error_pure")
  AllSessions=gather_pred_res(list(s1,s2,s3,s4,s5,s6,s7,s8),AllSessions$ID,AllSessions$ExpIndex,av_ori,AllSessions, "S.Error_ori")
  
  
  f3E=paired_comparison(AllSessions$TGTdir.pure_p.sn.Wmag,
                        AllSessions$TGTdir.ori_p.sn.Wmag,
                        AllSessions$ID,
                        AllSessions$ExpIndex)
  f3E_a=f3E[1,1]
  f3E_b=f3E[1,2]
  f3E_p=f3E[1,3]
  f3E_v=f3E[1,4]
  
  f4B_green=paired_comparison(AllSessions$P_L.V,
                              AllSessions$P_A.V, 
                              AllSessions$ID,
                              AllSessions$ExpIndex)
  f4B_g_a=f4B_green[1,1]
  f4B_g_b=f4B_green[1,2]
  f4B_g_v=f4B_green[1,4]
  f4B_g_p=f4B_green[1,3]
  
  
  f4B_red=paired_comparison(AllSessions$LOSc_L.V,
                            AllSessions$LOSc_A.V, 
                            AllSessions$ID,
                            AllSessions$ExpIndex)
  f4B_r_a=f4B_red[1,1]
  f4B_r_b=f4B_red[1,2]
  f4B_r_v=f4B_red[1,4]
  f4B_r_p=f4B_red[1,3]
  
  
  f4C=model_selection(AllSessions$perceived_av_Sound*180/pi,1/AllSessions$interval,
                      AllSessions$ID,AllSessions$ExpIndex,
                      1)
  f4C_v=fixed.effects(f4C[[1]])[2]
  f4C_p=f4C[[2]][2,6]
  
  f4D=model_selection(AllSessions$perceived_av_Sound*180/pi,AllSessions$widthHV,
                      AllSessions$ID,AllSessions$ExpIndex,
                      1)
  f4D_v=fixed.effects(f4D[[1]])[2]
  f4D_p=f4D[[2]][2,6]
  
  f4E=model_selection(AllSessions$perceived_av_Flight*180/pi,AllSessions$LOSl_F_V.LeadAngle,
                      AllSessions$ID,AllSessions$ExpIndex,
                      1)
  f4E_slope=fixed.effects(f4E[[1]])[2]
  f4E_intercept=fixed.effects(f4E[[1]])[1]
  f4E_p=f4E[[2]][2,6]
  
  
  delay_result=rbind(delay_result,(cbind(feedbackdelay_Sound,
                                         feedbackdelay_Flight,
                                         
                                         f2C_a,
                                         f2C_b,
                                         f2C_v,
                                         f2C_p,
                                         
                                         f2D_slope,
                                         f2D_intercept,
                                         f2D_p,
                                         
                                         f2J_lead_v,
                                         f2J_lead_p,
                                         f2J_side_v,
                                         f2J_side_p,
                                         
                                         sim000=mean(na.omit(abs(AllSessions$sim000.S.Error))),
                                         sim001=mean(na.omit(abs(AllSessions$sim001.S.Error))),
                                         sim002=mean(na.omit(abs(AllSessions$sim002.S.Error))),
                                         sim010=mean(na.omit(abs(AllSessions$sim010.S.Error))),
                                         sim011=mean(na.omit(abs(AllSessions$sim011.S.Error))),
                                         sim012=mean(na.omit(abs(AllSessions$sim012.S.Error))),
                                         sim020=mean(na.omit(abs(AllSessions$sim020.S.Error))),
                                         sim021=mean(na.omit(abs(AllSessions$sim021.S.Error))),
                                         sim022=mean(na.omit(abs(AllSessions$sim022.S.Error))),
                                         sim100=mean(na.omit(abs(AllSessions$sim100.S.Error))),
                                         sim101=mean(na.omit(abs(AllSessions$sim101.S.Error))),
                                         sim102=mean(na.omit(abs(AllSessions$sim102.S.Error))),
                                         sim110=mean(na.omit(abs(AllSessions$sim110.S.Error))),
                                         sim111=mean(na.omit(abs(AllSessions$sim111.S.Error))),
                                         sim112=mean(na.omit(abs(AllSessions$sim112.S.Error))),
                                         sim120=mean(na.omit(abs(AllSessions$sim120.S.Error))),
                                         sim121=mean(na.omit(abs(AllSessions$sim121.S.Error))),
                                         sim122=mean(na.omit(abs(AllSessions$sim122.S.Error))),
                                         
                                         sim0000=mean(na.omit(abs(AllSessions$S.Error_pure1))),
                                         sim0001=mean(na.omit(abs(AllSessions$S.Error_ori1))),
                                         sim0010=mean(na.omit(abs(AllSessions$S.Error_pure2))),
                                         sim0011=mean(na.omit(abs(AllSessions$S.Error_ori2))),
                                         sim0100=mean(na.omit(abs(AllSessions$S.Error_pure3))),
                                         sim0101=mean(na.omit(abs(AllSessions$S.Error_ori3))),
                                         sim0110=mean(na.omit(abs(AllSessions$S.Error_pure4))),
                                         sim0111=mean(na.omit(abs(AllSessions$S.Error_ori4))),
                                         sim1000=mean(na.omit(abs(AllSessions$S.Error_pure5))),
                                         sim1001=mean(na.omit(abs(AllSessions$S.Error_ori5))),
                                         sim1010=mean(na.omit(abs(AllSessions$S.Error_pure6))),
                                         sim1011=mean(na.omit(abs(AllSessions$S.Error_ori6))),
                                         sim1100=mean(na.omit(abs(AllSessions$S.Error_pure7))),
                                         sim1101=mean(na.omit(abs(AllSessions$S.Error_ori7))),
                                         sim1110=mean(na.omit(abs(AllSessions$S.Error_pure8))),
                                         sim1111=mean(na.omit(abs(AllSessions$S.Error_ori8))),
                                         
                                         f3E_a,
                                         f3E_b,
                                         f3E_v,
                                         f3E_p,
                                         
                                         f4B_g_a,
                                         f4B_g_b,
                                         f4B_g_v,
                                         f4B_g_p,
                                         
                                         f4B_r_a,
                                         f4B_r_b,
                                         f4B_r_v,
                                         f4B_r_p,
                                         
                                         f4C_v,
                                         f4C_p,
                                         f4D_v,
                                         f4D_p,
                                         
                                         f4E_slope,
                                         f4E_intercept,
                                         f4E_p)
  ))
  
  print(paste("LOOP:",k))
}

setwd(paste0(basedir,"/results/sim of delay effect"))
save(delay_result,file="delay_result_70_170ms.Rdata")

setwd(initialdir)
