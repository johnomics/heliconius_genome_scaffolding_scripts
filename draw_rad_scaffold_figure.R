#!/usr/bin/env Rscript

# This script is provided FOR INFORMATION ONLY and is not intended for use.
# It was written specifically for the Heliconius melpomene genome paper
# (Heliconius Genome Consortium, doi: 10.1038/nature11041) and has not
# been adapted for general use.

# Purpose: Generate input data for draw_rad_scaffold_figure.R,
#          used to produce Supplementary Figure S4.6.1

# Input  : filename for PDF output
#          files from data_for_rad_scaffold_figure.pl for one chromosome:
#            scaffold start and end positions
#            marker positions
#            cM positions
#          synteny data in tab-delimited data frame format,  eg:
# Hmel.gene       Hmel.pos        Bmor.gene       Bmor.chr        Bmor.pos
# HMEL009995-PA   29604   BGIBMGA011261-PA        chr23   17151375
# HMEL013926-PA   14393257        BGIBMGA013392-PA        chr27   3181868

# Output : PDF of figure

# Author: John Davey john.davey@ed.ac.uk
# Begun 23/10/2011

#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################

# Initialise
library(grid)
library(RColorBrewer)

args<-commandArgs(trailingOnly=T)
pdf(args[1], 5, 10)

grid.newpage()
grid.text("cM",0.08,0.98,just=c("left","bottom"),gp=gpar(fontsize=10))
grid.text("Mb",0.4,0.98,just=c("left","bottom"),gp=gpar(fontsize=10))
grid.text("Mb",0.9,0.98,just=c("left","bottom"),gp=gpar(fontsize=10))

pushViewport(viewport(x=0.5,y=0.49,width=0.95,height=0.95, gp=gpar(lineend="butt")))

# Load scaffold sizes
read.delim(args[2])->chr18.scf.pos

# Load marker positions
read.delim(args[3])->chr18.marker.chrpos.cm

# Load genetic positions and draw genetic map
read.delim(args[4])->chr18.cm.pos

# Load Bombyx synteny data
read.delim(args[5])->chr18.genes

genetic.colours<-c(rgb(244,165,130,max=255),rgb(5,113,176,max=255))

# Add alternating groups for colouring
chr18.cm.pos<-data.frame(chr18.cm.pos,Colour=rep(genetic.colours,length(chr18.cm.pos$cM)/2))

# Draw genetic map
viewport(0,0,just=c("left","bottom"),width=0.3,height=1,yscale=c(max(chr18.cm.pos$cM),min(chr18.cm.pos$cM)),name="chr18.genetic")->geneticvp
pushViewport(geneticvp)
grid.text(sprintf("%5.2f",chr18.cm.pos$cM),0.35,unit(chr18.cm.pos$cM,"native"),just=c("right","centre"),gp=gpar(fontsize=8))
grid.lines(c(0.6,0.6),unit(c(max(chr18.cm.pos$cM),min(chr18.cm.pos$cM)),"native"),gp=gpar(col=rgb(141,160,203,max=255),lwd=5,lineend="round"))
grid.polyline(
	c(rep(0.4,length(chr18.cm.pos$cM)),rep(0.8,length(chr18.cm.pos$cM))),
	unit(c(chr18.cm.pos$cM,chr18.cm.pos$cM),"native"),
	id=rep(1:length(chr18.cm.pos$cM),2),
	gp=gpar(col=genetic.colours,lwd=3,lineend="round")
)

# Create physical map viewport
popViewport()
viewport(0.35,0,just=c("left","bottom"),width=0.3,height=1,yscale=c(max(as.numeric(chr18.scf.pos$End)),1),name="chr18.physical")->physicalvp
pushViewport(physicalvp)

# Draw connections between genetic and physical map
apply (chr18.cm.pos,1,
	function(x) {
		popViewport()
		pushViewport(geneticvp)
		grid.move.to(0.8,unit(x[1],"native"))
		popViewport()
		pushViewport(physicalvp)
		midpoint = (as.numeric(x[2])+as.numeric(x[3]))/2
		grid.line.to(0,unit(midpoint,"native"),gp=gpar(col=x[4],lwd=1.5,lty="dashed"))
		grid.lines(c(0,0),unit(c(as.numeric(x[2])+20000,as.numeric(x[3])-20000),"native"),gp=gpar(col=x[4],lwd=3,lineend="round"))
	}
)


# Draw markers
paternal.markers<-merge(chr18.marker.chrpos.cm,chr18.cm.pos)
paternal.marker.num<-length(paternal.markers$Position)
grid.polyline(
	c(rep(0.075,paternal.marker.num),rep(0.125,paternal.marker.num)),
	unit(c(paternal.markers$Position,paternal.markers$Position),"native"),
	id=rep(1:paternal.marker.num,2),
	gp=gpar(col=as.character(paternal.markers$Colour),lineend="round")
)


# Draw Bombyx chromosomes
chr23.length<-23133604
chr27.length<-14467522
chr23.prop<-chr23.length/(chr23.length+chr27.length)
chr23.vp<-viewport(0.75, 1-chr23.prop+0.005, width=0.25, height=chr23.prop-0.005, just=c("left","bottom"), name="Bmor.chr23", yscale=c(1,chr23.length))
chr27.vp<-viewport(0.75, 0, width=0.25, height=1-chr23.prop-0.005, just=c("left","bottom"), name="Bmor.chr27", yscale=c(1,chr27.length))

# Draw gene connections
apply(chr18.genes[((chr18.genes$Bmor.chr=='chr23') & ((chr18.genes$Hmel.pos<843357) | (chr18.genes$Hmel.pos>1445632))),],1,
    function(x) {
        popViewport();
        pushViewport(physicalvp);
        grid.move.to(0.9,unit(x[2],"native"));
        popViewport();
        pushViewport(chr23.vp);
        grid.line.to(0.4,unit(x[5],"native"),gp=gpar(col=rgb(141,160,203,max=255)));
    }
)
apply(chr18.genes[((chr18.genes$Bmor.chr=='chr23') & ((chr18.genes$Hmel.pos>=843357) & (chr18.genes$Hmel.pos<=1445632))),],1,
    function(x) {
        popViewport()
        pushViewport(physicalvp)
        grid.move.to(0.9,unit(x[2],"native"))
        popViewport()
        pushViewport(chr23.vp)
#         grid.line.to(0.4,unit(x[5],"native"),gp=gpar(col="red"))  # BD
		grid.line.to(0.4,unit(x[5],"native"),gp=gpar(col=rgb(141,160,203,max=255)))
    }
)

apply(chr18.genes[chr18.genes$Bmor.chr=='chr27',],1,
    function(x) {
        popViewport();
        pushViewport(physicalvp);
        grid.move.to(0.9,unit(x[2],"native"));
        popViewport();
        pushViewport(chr27.vp);
        grid.line.to(0.4,unit(x[5],"native"),gp=gpar(col=rgb(141,160,203,max=255)));
    }
)

# Draw chromosome 23
popViewport()
pushViewport(chr23.vp)

# Chr 23 axes
chr23.onemil.tick<-seq(8000,chr23.length,1000000)
chr23.tenmil.tick<-seq(8000,chr23.length,10000000)
chr23.fivemil.tick<-seq(8000,chr23.length,5000000)
grid.text(sprintf("%2d", 0:23), 0.8, unit(chr23.onemil.tick, "native"), just="right", gp=gpar(fontsize=8))
grid.polyline(
	c(rep(0.6,24),rep(0.65,24)),
	unit(c(chr23.onemil.tick,chr23.onemil.tick),"native"),
	id=rep(1:24,2),
	gp=gpar(col="grey",lwd=3,lineend="round")
)
grid.polyline(
	c(rep(0.6,5),rep(0.65,5)),
	unit(c(chr23.fivemil.tick,chr23.fivemil.tick),"native"),
	id=rep(1:5,2),
	gp=gpar(col="dimgrey",lwd=4,lineend="round")
)
grid.polyline(
	c(rep(0.6,3),rep(0.65,3)),
	unit(c(chr23.tenmil.tick,chr23.tenmil.tick),"native"),
	id=rep(1:3,2),
	gp=gpar(col="black",lwd=5,lineend="round")
)


grid.rect(x=0.41,width=0.19,y=0,height=1,just=c("left","bottom"),gp=gpar(col="black",fill=rgb(186,186,186,max=255),lwd=4,lineend="round"))
grid.text("B. mori 23", 0.93,0.5,gp=gpar(font=3, fontsize=12),rot=90)



# Draw chromosome 27
popViewport()
pushViewport(chr27.vp)

# Chr 27 axes
chr27.onemil.tick<-seq(8000,chr27.length,1000000)
chr27.tenmil.tick<-seq(8000,chr27.length,10000000)
chr27.fivemil.tick<-seq(8000,chr27.length,5000000)

grid.text(sprintf("%2d", 0:14), 0.8, unit(chr27.onemil.tick, "native"), just="right", gp=gpar(fontsize=8))
grid.polyline(
	c(rep(0.6,15),rep(0.65,15)),
	unit(c(chr27.onemil.tick,chr27.onemil.tick),"native"),
	id=rep(1:15,2),
	gp=gpar(col="grey",lwd=3,lineend="round")
)
grid.polyline(
	c(rep(0.6,3),rep(0.65,3)),
	unit(c(chr27.fivemil.tick,chr27.fivemil.tick),"native"),
	id=rep(1:3,2),
	gp=gpar(col="dimgrey",lwd=4,lineend="round")
)
grid.polyline(
	c(rep(0.6,2),rep(0.65,2)),
	unit(c(chr27.tenmil.tick,chr27.tenmil.tick),"native"),
	id=rep(1:2,2),
	gp=gpar(col="black",lwd=5,lineend="round")
)
grid.rect(x=0.41,width=0.19,y=0,height=1,just=c("left","bottom"),gp=gpar(col="black",fill=rgb(186,186,186,max=255),lwd=4,lineend="round"))
grid.text("B. mori 27",  0.93, 0.5, gp=gpar(font=3, fontsize=12),rot=90)


# Draw scaffolds
popViewport()
pushViewport(physicalvp)

mb.onemil.tick<-seq(8000,16000000,1000000)
mb.tenmil.tick<-seq(8000,16000000,10000000)
mb.fivemil.tick<-seq(8000,16000000,5000000)
grid.text(sprintf("%2d", 0:15), 0.29, unit(mb.onemil.tick, "native"), just="right", gp=gpar(fontsize=8))
grid.polyline(
	c(rep(0.35,16),rep(0.4,16)),
	unit(c(mb.onemil.tick,mb.onemil.tick),"native"),
	id=rep(1:16,2),
	gp=gpar(col="grey",lwd=1,lineend="round")
)
grid.polyline(
	c(rep(0.35,4),rep(0.4,4)),
	unit(c(mb.fivemil.tick,mb.fivemil.tick),"native"),
	id=rep(1:4,2),
	gp=gpar(col="darkgrey",lwd=2,lineend="round")
)
grid.polyline(
	c(rep(0.35,2),rep(0.4,2)),
	unit(c(mb.tenmil.tick,mb.tenmil.tick),"native"),
	id=rep(1:2,2),
	gp=gpar(col="dimgrey",lwd=3,lineend="round")
)


chr.width=0.49
chr.start=0.4
grid.rect(
	y=unit(c(chr18.scf.pos$Start),"native"),
	height=unit(c(chr18.scf.pos$End-chr18.scf.pos$Start),"native"),
	x=rep(chr.start,length(chr18.scf.pos$Scaffold)),
	width=rep(chr.width,length(chr18.scf.pos$Scaffold)),
	just=c("left","bottom"),
    gp=gpar(fill=rgb(141,160,203,max=255),col="black")
)

# Draw BD scaffold in red (red band scaffold)
grid.rect(y=unit(844157,"native"),height=unit(602275,"native"),x=chr.start,width=chr.width, just=c("left","bottom"), gp=gpar(fill="red",col="black",lwd=3,lineend="round"))


ordered.scfs<-chr18.scf.pos[chr18.scf.pos$Link=="black",]
grid.polyline(
    y=unit(c(ordered.scfs$Start,ordered.scfs$End), "native"),
    x=c(rep(chr.start,length(ordered.scfs$Scaffold)),rep(chr.start,length(ordered.scfs$Scaffold))),
    id=rep(1:length(ordered.scfs$Scaffold),2),
    gp=gpar(lwd=4,col=as.character(ordered.scfs$Link),lineend="round")
)
grid.polyline(
    y=unit(c(ordered.scfs$Start,ordered.scfs$End), "native"),
    x=c(rep(chr.start+chr.width,length(ordered.scfs$Scaffold)),rep(chr.start+chr.width,length(ordered.scfs$Scaffold))),
    id=rep(1:length(ordered.scfs$Scaffold),2),
    gp=gpar(lwd=4,col=as.character(ordered.scfs$Link),lineend="round")
)

dev.off()
