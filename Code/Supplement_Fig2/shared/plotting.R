
suppressPackageStartupMessages(library(grid))
labels <- c("ImmuCellAI 2.0", "BayesPrism", "CIBERSORT", "CITMIC", "DWLS", "ImmuCellAI", "MuSiC")
titles <- c("Liu's cohort (Lymph node)", "Zhang's cohort (CRC)", "Wu's cohort (BRCA)", "Liu's cohort (NPC)")
short_titles <- c("Liu (Lymph node)", "Zhang (CRC)", "Wu (BRCA)", "Liu (NPC)")
palette <- c("#B21D39", "#547FA5", "#E28C36", "#66A7A0", "#8266A5", "#599444", "#B28A50")
red <- "#8D222C"; black <- "#202020"
text_at <- function(x, y, value, size = 12, bold = FALSE, color = black, align = "left", angle = 0) {
 grid.text(value, x, y, default.units = "native", just = c(align, "bottom"), rot = angle,
           gp = gpar(fontsize = size, fontface = if (bold) "bold" else "plain", col = color, fontfamily = "Arial"))
}
rect_at <- function(x,y,w,h,fill="white",col=black,lwd=1.1)
 grid.rect(x,y,w,h,default.units="native",just=c("left","bottom"),gp=gpar(fill=fill,col=col,lwd=lwd))
line_at <- function(x,y,x2,y2,col=black,lwd=1)
 grid.lines(c(x,x2),c(y,y2),default.units="native",gp=gpar(col=col,lwd=lwd))
r_color <- function(r) {
 rgb <- colorRamp(c("#4575B4","#66C2A5","#FFFFBF","#F46D43","#9E0142"))((pmin(1,pmax(-1,r))+1)/2)
 grDevices::rgb(rgb[1],rgb[2],rgb[3],maxColorValue=255)
}
pie_at <- function(x,y,r,status,radius=11.6) {
 if(status=="unavailable"){line_at(x-3,y-6,x+3,y+6,lwd=1.15);return(invisible(NULL))}
 if(!is.finite(r)){text_at(x,y-4,"x",13,TRUE,align="center");return(invisible(NULL))}
 grid.circle(x,y,radius,default.units="native",gp=gpar(fill="white",col="#333333",lwd=1.15))
 if(abs(r)>0){
  a <- seq(pi/2,pi/2-2*pi*min(1,abs(r)),length.out=max(4,ceiling(180*abs(r))))
  grid.polygon(c(x,x+radius*cos(a)),c(y,y+radius*sin(a)),default.units="native",gp=gpar(fill=r_color(r),col="#333333",lwd=1.15))
 }
 grid.circle(x,y,radius,default.units="native",gp=gpar(fill=NA,col="#333333",lwd=1.15))
}
draw_pies <- function(panel,left,right,top,lx) {
 result <- analysis[[panel]]; d <- result$data; cells <- result$cells
 bottom <- top-44-7*29; pitch <- (right-left)/nrow(cells); k <- match(panel,LETTERS[1:4])
 rect_at(left,top-44,right-left,44,fill="#F5EBEB")
 text_at(lx,top-20,panel,23,TRUE)
 text_at((left+right)/2,top-31,titles[k],19,TRUE,red,"center")
 rect_at(left,bottom,right-left,7*29)
 for(j in seq_along(methods)){
   y <- top-44-(j-.5)*29
   text_at(left-10,y-5,labels[j],14,j==1,if(j==1)red else black,"right")
   for(i in seq_len(nrow(cells))){
     z <- d[d$method==methods[j]&d$target==cells$target[i],]
     stopifnot(nrow(z)==1)
     pie_at(left+(i-.5)*pitch,y,z$Pearson,z$status)
   }
 }
 for(i in seq_len(nrow(cells)))text_at(left+(i-.5)*pitch,bottom-12,cells$display_label[i],12.4,align="right",angle=48)
}
draw_legend <- function(top=1060) {
 x <- 1090; height <- 120
 for(i in 0:239)rect_at(x,top-height+i*height/240,17,height/240+.1,fill=r_color(-1+2*i/239),col=NA)
 text_at(x-10,top-56,"Pearson r",12,TRUE,align="right")
 for(v in c(-1,-.5,0,.5,1))text_at(x+24,top-height+(v+1)*height/2-3,sprintf("%.1f",v),10)
 text_at(1010,top-162,"Symbol key",13,TRUE,red)
 text_at(1010,top-190,"x",15,TRUE);text_at(1030,top-187,"Undefined correlation",10.5)
 line_at(1010,top-226,1019,top-212,lwd=1.4)
 text_at(1030,top-218,"Unavailable / no",10.5);text_at(1030,top-234,"matching reference",10.5)
}
draw_mean <- function(top=382) {
 d <- analysis$E
 text_at(4,top,"E",22,TRUE);text_at(45,top-17,"Pearson correlation across cell types: mean +/- SD",19,TRUE,red)
 lo <- min(-1.05,min(d$mean_Pearson-d$SD_Pearson)-.04)
 hi <- max(1.05,max(d$mean_Pearson+d$SD_Pearson)+.04)
 y <- top-254; h <- 178; w <- 230
 for(k in 1:4){
  x <- c(76,364,652,940)[k]; yp <- function(r)y+(r-lo)/(hi-lo)*h
  text_at(x+w/2,top-50,short_titles[k],14,TRUE,red,"center")
  for(v in c(-1,-.5,0,.5,1)){
   line_at(x,yp(v),x+w,yp(v),col="#DADADA",lwd=.6);text_at(x-7,yp(v)-3,format(v,trim=TRUE),11,align="right")
  }
  line_at(x,y,x,y+h);line_at(x,y,x+w,y)
  if(k==1)text_at(x-45,y+h/2,"Mean Pearson r",12,TRUE,align="center",angle=90)
  for(j in seq_along(methods)){
   row <- d[d$dataset==datasets[k]&d$method==methods[j],];xx<-x+(j-.5)*w/7
   line_at(xx,yp(row$mean_Pearson-row$SD_Pearson),xx,yp(row$mean_Pearson+row$SD_Pearson),col=palette[j])
   for(v in row$mean_Pearson+c(-1,1)*row$SD_Pearson)line_at(xx-3.5,yp(v),xx+3.5,yp(v),col=palette[j])
   yy<-yp(row$mean_Pearson)
   grid.polygon(c(xx,xx+4.6,xx,xx-4.6),c(yy+4.6,yy,yy-4.6,yy),default.units="native",gp=gpar(fill=palette[j],col=palette[j]))
   text_at(xx+3,y-10,labels[j],12,j==1,if(j==1)red else black,"right",58)
  }
 }
}
make_page <- function(height,draw) function(){
 grid.newpage();pushViewport(viewport(xscale=c(0,1200),yscale=c(0,height),gp=gpar(fontfamily="Arial")))
 draw();popViewport()
}
build_panel <- function(panel) {
 if(panel=="E")list(plot=make_page(430,function()draw_mean(410)),width=1200/72,height=430/72)
 else list(plot=make_page(410,function()draw_pies(panel,165,1173,395,18)),width=1200/72,height=410/72)
}
combine_panels <- function(plots) make_page(1448,function(){
 draw_pies("A",155,1198,1445,4);draw_pies("B",155,990,1095,4)
 draw_pies("C",155,514,745,4);draw_pies("D",755,1198,745,604)
 draw_legend(1075);draw_mean(405)
})
save_plot <- function(plot,stem,width,height){
 result<-character()
 for(ext in if(output_format=="both")c("pdf","png")else output_format){
  path<-file.path(output_dir,paste0(stem,".",ext))
  if(ext=="pdf")cairo_pdf(path,width=width,height=height,family="Arial",bg="white")
  else png(path,width=width,height=height,units="in",res=dpi,type="cairo",family="Arial",bg="white")
  tryCatch(plot(),finally=invisible(dev.off()));result<-c(result,path)
 }
 result
}
