#!/usr/bin/Rscript --no-save --no-restore  
## novoQplot.R
## Colin Hercus and Zayed Albertyn
## SYNOPSIS: Plot quality values 
## Usage: novoqplot.R  [input qcal.csv file]  [output pdf filename]
## e.g.  novoQplot.R   617.qcal.csv 617.qcal.pdf

#need to install the ggplot2,reshape and gridExtra packages as below:
#install.packages("gridExtra")
#install.packages("reshape")
#install.packages("ggplot2")

# Changes:
#     6Jun11 Colin - Fix problems when not enough data to calculate some values (NaNs & Infs)
print("Load Libraries")
library(gridExtra)
library(reshape)
library(ggplot2)
print("Libraries Loaded")
inputfile = "qcal.tsv"
outfile="qcal.pdf"
args=commandArgs(trailingOnly = TRUE)
inputfile= args[1]
outfile= args[2]

print(paste("Reading in input file",inputfile))
print("Read Data")
data = readLines(inputfile)
command = ""
if(regexpr("^# ", data[1]) != -1) { 
	command = sub("# (.*)","\\1",data[1]) 
	data = data[-1]  
}
print(command)
qcal <- read.csv(textConnection(data), header = TRUE, sep = "\t")
qcal$Side = factor(qcal$Side)

nside2 <- table(qcal$Side)["1"]

colourSpace = names(qcal)[2] != "Index"
if (colourSpace) {
	colnames(qcal)[2] <- "Index"
    colnames(qcal)[3] <- "Quality"
    colnames(qcal)[4] <- "Base"
	qcal$Base <- factor(qcal$Base)
	levels(qcal$Base) <- c("R","Y","G","B")
    colnames(qcal)[5] <- "calls"
    colnames(qcal)[6] <- "mismatch"
    Base="Colour"
} else {
	qcal$match <- ifelse(qcal$Base == "A", qcal$A.s, ifelse(qcal$Base == "C", qcal$C.s, ifelse(qcal$Base == "G", qcal$G.s, qcal$T.s)))
	qcal$calls <-  qcal$A.s + qcal$C.s + qcal$G.s + qcal$T.s + 0.000001
	qcal$mismatch <- qcal$calls - qcal$match
	Base = "Base"
}

title.size<-10

mxqcalled <- ((max(qcal$Quality[qcal$calls + qcal$NMs > 0]) + 9) %/% 10) * 10

qcal$emismatch <- (qcal$calls + qcal$NMs) * 10^(-qcal$Quality/10)
qcal$quality <- as.integer(0.5 -10 * log((1+qcal$mismatch)/(qcal$calls+10^(qcal$Quality/10)),10))

mxqactual <- ((max(qcal[(qcal$calls > 3 * 10^(qcal$Quality/10) | qcal$mismatch > 2),]$quality) + 9) %/% 10) * 10
if (is.nan(mxqactual)) mxqactual <- mxqcalled
print(mxqactual)
print("Generating plots")

fontsize <- 6
my_title= "Read1 Called vs Empirical Quality"
titleobj<- opts(title = my_title,plot.title=theme_text(size=title.size))
print(my_title)
#, gpar(col = "blue", lty = "solid", lwd = 3, fontsize = 7)
p <- ggplot(qcal[qcal$Side==0 & (qcal$calls > 3 * 10^(qcal$Quality/10) | qcal$mismatch > 2),], aes(x=Quality, y=quality,colour=Base))
ca1 <- p + geom_jitter(size=0.5, position=position_jitter(height=1)) + scale_colour_manual(name=Base, values = c("red", "orange", "green","blue")) + titleobj
ca1 <- ca1 + scale_y_continuous("Empirical Quality", limits=c(0, mxqactual)) 
ca1 <- ca1 + scale_x_continuous("As Called Quality", limits=c(0, mxqcalled))
ca1 <- ca1 + opts(legend.position=c(0.2,0.7))



mqcal <- melt(qcal,id=c("Side","Quality","Index","Base"))

# Select one read side and cast as sum over Index & base
#Actual Qualities across the reads
ib <- cast(mqcal[mqcal$Side==0,] , Index + Base ~ variable, sum)

ib$quality <- -10 * log(ib$mismatch/ib$calls, 10)
ib$equality <- -10 * log(ib$emismatch/(ib$calls + ib$NMs), 10)

mxqactual2 <- ((max(ib[is.finite(ib$quality) & ib$mismatch > 2,]$quality) + 9) %/% 10) * 10
#if (is.nan(mxqactual2)) mxqactual2 <- mxqcalled
print(mxqactual2)

my_title= paste("Read 1 - Empirical Quality by",Base,"Position")
print(my_title)
titleobj<- opts(title = my_title,plot.title=theme_text(size=title.size))
p <- ggplot(ib[is.finite(ib$quality) & ib$mismatch > 2,], aes(x=Index, y=quality,colour=Base))
cb1<-  p + geom_point(size=1) + scale_colour_manual(values = c("red", "orange", "green","blue")) + titleobj  
cb1 <- cb1 + opts(legend.position="none")
cb1 <- cb1 + scale_y_continuous("Empirical Quality", limits=c(0, mxqactual2)) 
cb1 <- cb1 + scale_x_continuous("Position in Read")

my_title= paste("Read 1 - Called Quality by",Base,"Position")
print(my_title)
titleobj<- opts(title = my_title,plot.title=theme_text(size=title.size))
p <- ggplot(ib, aes(x=Index, y=equality,colour=Base))
ab1 <- p + geom_point(size=1) + scale_colour_manual(values = c("red", "orange", "green","blue")) + titleobj
ab1 <- ab1 + opts(legend.position="none")
ab1 <- ab1 + scale_y_continuous("As Called Quality", limits=c(0, mxqcalled)) 
ab1 <- ab1 + scale_x_continuous("Position in Read")

my_title= paste("Read 1 -",Base,"Calls by",Base,"Position")
print(my_title)
titleobj<- opts(title = my_title,plot.title=theme_text(size=title.size))
p <- ggplot(ib, aes(x=Index, y=calls+NMs,colour=Base))
lb1<- p + geom_point(size=1) + scale_colour_manual(values = c("red", "orange", "green","blue")) + titleobj
lb1 <- lb1 + opts(legend.position="none")
lb1 <- lb1 + scale_y_continuous("Number of Calls") 
lb1 <- lb1 + scale_x_continuous("Position in Read")

if( is.na(nside2) ) {
    print("Single End Reads")
    print(paste("Writing output to",outfile))
    #Open PDF output
    pdf(outfile,paper="a4",pointsize=12, width=0, height=0)
    # The theme
    old_theme <- theme_update(axis.text.x=theme_text(size=fontsize),axis.text.y=theme_text(size=fontsize),
                          axis.title.x=theme_text(size=fontsize),axis.title.y=theme_text(size=fontsize,angle=90),
                          strip.text.x=theme_text(size=fontsize),strip.text.y=theme_text(size=fontsize),
                          legend.text=theme_text(size=fontsize),legend.title=theme_text(size=fontsize),
    legend.key.size = unit(0.4, "lines"), legend.key = theme_rect(colour = "transparent")) 
    #Arrange the plots
    grid.arrange(ca1,cb1,ab1,lb1, ncol=2,main=paste("Novoalign Quality Calibration Plot -", inputfile), sub=textGrob(paste(command,"\n",date()), gp=gpar(font=2,fontsize=fontsize)))
    dev.off()    
    q()
};
#  Side 2
print("Paired End Reads")
my_title= "Read2 Called vs Empirical Quality"
print(my_title)
titleobj<- opts(title = my_title,plot.title=theme_text(size=title.size))
p <- ggplot(qcal[qcal$Side==1 & (qcal$calls > 3 * 10^(qcal$Quality/10) | qcal$mismatch > 2),], aes(x=Quality, y=quality,colour=Base))
ca2 <-p + geom_jitter(size=0.5, position=position_jitter(height=1))  + scale_colour_manual(values = c("red", "orange", "green","blue")) + titleobj
ca2 <- ca2 + scale_y_continuous("Empirical Quality", limits=c(0, mxqactual)) 
ca2 <- ca2 + scale_x_continuous("As Called Quality", limits=c(0, mxqcalled))
ca2 <- ca2 + opts(legend.position="none")

print("Pivot Tables")
# Select one read side and cast as sum over Index & base
#Actual Qualities across the reads
ib <- cast(mqcal[mqcal$Side==1,] , Index + Base ~ variable, sum)
ib$equality <- -10 * log(ib$emismatch/(ib$calls + ib$NMs), 10)
ib$quality <- -10 * log(ib$mismatch/ib$calls, 10)

mxqactual2 <- ((max(ib[is.finite(ib$quality) & ib$mismatch > 2,]$quality) + 9) %/% 10) * 10
if (is.nan(mxqactual2)) mxqactual2 <- mxqcalled
print(mxqactual2)

my_title= paste("Read 2 - Empirical Quality by",Base,"Position")
print(my_title)
titleobj<- opts(title = my_title,plot.title=theme_text(size=title.size))
p <- ggplot(ib[is.finite(ib$quality) & ib$mismatch > 2,], aes(x=Index, y=quality,colour=Base))
cb2<-  p + geom_point(size=1) + scale_colour_manual(values = c("red", "orange", "green","blue")) + titleobj
cb2 <- cb2 + opts(legend.position="none")
cb2 <- cb2 + scale_y_continuous("Empirical Quality", limits=c(0, mxqactual2)) 
cb2 <- cb2 + scale_x_continuous("Position in Read")

my_title= paste("Read 2 - Called Quality by",Base,"Position")
print(my_title)
titleobj<- opts(title = my_title,plot.title=theme_text(size=title.size))
p <- ggplot(ib, aes(x=Index, y=equality,colour=Base))
ab2 <- p + geom_point(size=1) + scale_colour_manual(values = c("red", "orange", "green","blue")) + titleobj
ab2 <- ab2 + opts(legend.position="none")
ab2 <- ab2 + scale_y_continuous("As Called Quality", limits=c(0, mxqcalled)) 
ab2 <- ab2 + scale_x_continuous("Position in Read")

my_title= paste("Read 2 -",Base,"Calls by",Base,"Position")
print(my_title)
titleobj<- opts(title = my_title,plot.title=theme_text(size=title.size))
p <- ggplot(ib, aes(x=Index, y=calls+NMs,colour=Base))
lb2<- p + geom_point(size=1) + scale_colour_manual(values = c("red", "orange", "green","blue")) + titleobj
lb2 <- lb2 + opts(legend.position="none")
lb2 <- lb2 + scale_y_continuous("Number of Calls") 
lb2 <- lb2 + scale_x_continuous("Position in Read")

print(paste("Writing output to",outfile))
#Open PDF output
pdf(outfile,paper="a4",pointsize=12, width=0, height=0)
# The theme
old_theme <- theme_update(axis.text.x=theme_text(size=fontsize),axis.text.y=theme_text(size=fontsize),
                          axis.title.x=theme_text(size=fontsize),axis.title.y=theme_text(size=fontsize,angle=90),
                          strip.text.x=theme_text(size=fontsize),strip.text.y=theme_text(size=fontsize),
                          legend.text=theme_text(size=fontsize),legend.title=theme_text(size=fontsize),
legend.key.size = unit(0.4, "lines"), legend.key = theme_rect(colour = "transparent")) 
#Arrange the plots
grid.arrange(ca1,ca2,cb1,cb2,ab1,ab2,lb1,lb2, ncol=2,main=paste("Novoalign Quality Calibration Plot -", inputfile), sub=textGrob(paste(command,"\n",date()), gp=gpar(font=2,fontsize=fontsize)))
dev.off()
