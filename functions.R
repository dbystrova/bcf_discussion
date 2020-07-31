
lissage  =function(m){
  nr = dim(m)[1]
  nc = dim(m)[2]
  for (i in 1:(nr)){
    for (j in 2:(nc-1)){
      v = m[i,(j-1):(j+1)]
      m[i,j] = mean(v)
    }
  }  
  for (j in 1:(nc)){
    for (i in 2:(nr-1)){
      v = m[(i-1):(i+1),j]
      m[i,j] = mean(v)
    }
  }
  m
}

# 3D-plot
build3ds1<-function(x,y,z,par1=""){
  z<-pmin(z,5)
  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("blue", "red") )
  nbcol <- 100
  color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol <- cut(zfacet, nbcol)
  
  sigma_param<-x
  theta_param<-y
  par(mar = c(2,1,2,.1), mgp = c(1,1,0), xaxs = "i", yaxs = "i", las= 1)
  persp(
    x = sigma_param, 
    y = theta_param, 
    z = z, 
    zlim = c(min(z),max(z)), 
    col = color[facetcol], 
    theta = -25, 
    phi = 25,
    ltheta = 120, 
    ticktype = "detailed", 
    shade = 0.3,
    xlab = "", ylab = "", zlab = "", 
    d = 5, r = 10,
    cex.axis = 1, cex.lab = 1.3, nticks = 3, main =par1)
  text(.13,-.37,expression(x[1]), cex = 1.5)
  text(-.3,-.35,expression(x[2]), cex = 1.5)
}
