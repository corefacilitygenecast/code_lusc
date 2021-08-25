#!/lustre/rde/user/rde_admin/bin/Rscript
require("getopt",quiet=T)
library(ggplot2)
spec = matrix(c(
	'help','h',0,'logical','this help',
	'snp_matrix','s',1,'character','snp matrix input',
	'group_file','g',1,'character','group information',
	'prefix_out','p','1','character','output file(pdf & png)',
	'anno','a','1','character','annotation file',
        'kegg','k','1','character','formated kegg annotation file'
),byrow=TRUE,ncol=5);
opt=getopt(spec)
#if(!is.null(opt$h) || is.null(opt$s) || is.null(opt$g) || is.null(opt$o)){
if(!is.null(opt$h) || is.null(opt$s) || is.null(opt$g) || is.null(opt$prefix_out) || is.null(opt$a) || is.null(opt$k) ){
        cat(paste(getopt(spec,usage = T)))
        q(status=1)
}

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
mat=read.table(opt$s,fill=T,header=TRUE,stringsAsFactors=FALSE,sep="\t",row.names=1,encoding='UTF-8',check.names=F)
mat<-as.data.frame(t(mat))
number=dim(mat)[2]
distance=dim(mat)[1]
number
distance
group_data<-as.data.frame(read.table(opt$g,header=1,sep="\t",encoding='UTF-8',check.names=F))
group_data<-group_data[order(group_data$Efficacy),]
match.id1<-match(rownames(mat),group_data$Tumor_Sample_Barcode,nomatch=0)
mat$Group<-group_data[match.id1,"Efficacy",drop=FALSE]
mat<-mat[order(mat$Group),]
mat.bak<-mat
mat<-mat[,1:number,drop=FALSE]
#if(number==1){
#	mat<-as.data.frame(mat)
#	colnames(mat)<-colnames(mat.bak)[1]
#	rownames(mat)<-rownames(mat.bak)
#}
mat[is.na(mat)] = ""
mat = t(as.matrix(mat))
head(mat)
pathway=read.table(opt$a,fill=T,header=TRUE,stringsAsFactors=FALSE,sep="\t",encoding='UTF-8',check.names=F)
match1<-match(pathway$gene,rownames(mat),nomatch=0)
newdata<-as.data.frame(mat[match1,,drop=F])
match.id2<-match(rownames(newdata),pathway$gene,nomatch=0)
#if(length(match.id2)>0){
newdata$KEGG<-pathway[match.id2,"pathway",drop=FALSE]
drawdata<-as.data.frame(newdata[order(newdata$KEGG),1:(dim(newdata)[2]-1)])
muttype <- c("nonsynonymous_SNV","frameshift_deletion","nonframeshift_deletion",
             "frameshift_insertion","nonframeshift_insertion","frameshift_substitution","nonframeshift_substitution","stopgain","stoploss",'UTR','splicing',"UNKNOWN")
col = RColorBrewer::brewer.pal(n = length(muttype), name = 'Paired')
col <- col[c(2, 1, 7, 4, 3, 6, 5, 8, 9, 10,11,12)]
names(col) = muttype
col
alter_fun_my = function(x, y, w, h, v) {
  n = sum(v)
  w = w * 0.95
  h = h * 0.9
  grid.rect(x, y, w*0.95, h*0.9, gp = gpar(fill = "white", col = NA))
  if(n){
	grid.rect(x, y - h*0.5 + 1:n/n*h, w*0.9, 1/n*h,gp = gpar(fill = col[names(which(v))], col = NA), just = "top")	
  }
}
#group_data<-as.data.frame(read.table(opt$g,header=1))
color_n = length(unique(group_data$Efficacy))

#if(color_n == 1){
#    col_nw = c("red")
#}else if(color_n == 2){
#    col_nw = c("red","blue")
#}else{
#    col_nw = RColorBrewer::brewer.pal(n=color_n, name="Set3")
#}
# if(color_n==1){
#        col_nw = "#6699CC"
#  }else if (color_n == 2){
#      col_nw = c("#6699CC","#FF9933")
 # }else if (color_n == 3){
#       col_nw = c("#6699CC","#FF9933","#006633")
#  }else if (color_n == 4){
#       col_nw = c("#6699CC","#FF9933","#006633","#CC3366")
#  }else if (color_n == 5){
#       col_nw = c("#6699CC","#FF9933","#006633","#CC3366","#FFFF66")
#  }else{
#      col_nw = RColorBrewer::brewer.pal(n=color_n, name="Set3")
#  }
 if(color_n==1){
        col_nw = "darkorange"
  }else if (color_n == 2){
      col_nw = c("#003399","#CC0033")
  }else{
      col_nw = RColorBrewer::brewer.pal(n=color_n, name="Set3")
  }

#color_nr=length(unique(newdata$KEGG))
# if(color_nr==1){
#        col_rw = "darkorange"
#  }else if (color_nr == 2){
#      col_rw = c("A6CEE3","#1F78B4")
#  }else{
#      col_rw = RColorBrewer::brewer.pal(n=color_nr, name="Set3")
#  }


group.assign<-setNames(col_nw,unique(group_data$Efficacy))
match.id<-match(colnames(mat),group_data$Tumor_Sample_Barcode,nomatch=0)
df<-group_data[match.id,"Efficacy",drop=FALSE]
ha<-HeatmapAnnotation(df=df,col=list(Efficacy=group.assign),annotation_height = unit(1, "cm"))

group1=unique(group_data$Efficacy)[1]
group2=unique(group_data$Efficacy)[2]
anno<-group_data
data1<-anno[which(anno$Efficacy==group1),]
data2<-anno[which(anno$Efficacy==group2),]
data1 <-data1[order(data1$Mutation_counts,decreasing=T),]
data2 <-data2[order(data2$Mutation_counts,decreasing=T),]
new_data <-rbind(data1,data2)

sample_order= as.character(new_data$Tumor_Sample_Barcode)
rownames(anno)<-anno$Tumor_Sample_Barcode
match.id<-match(colnames(drawdata),anno$Tumor_Sample_Barcode,nomatch=0)
temp_df<-anno[match.id,]
head(temp_df)
rownames(temp_df) = temp_df[,1]
#temp_df <-temp_df[order(temp_df$Smoking_status),]
#top_df<-
top_df<-temp_df[,c(3,5,6,7)]
#head(top_df)
col_Sex=c("#E5201F","#007FC4")
assign_Sex=setNames(col_Sex,unique(top_df$Sex))
col_Smoking=c("#E96B00","#721E81")
assign_Smoking=setNames(col_Smoking,unique(top_df$Smoking))
#col_Age =c("#0099CC","#99CC33")
#col_Age =c("#E7392A","#9CD5E3")
col_Age=pal_npg("nrc",alpha=0.8)(2)
assign_Age=setNames(col_Age,unique(top_df$Age))
top_df_1 <-temp_df[sample_order,]
top_df_2 <-top_df_1[,c(8,9)]
#top_df_2 <-top_df_1[1]
#head(top_df_2)
#top_df_3 <-tmep_df[9]
#top_ha = HeatmapAnnotation(df=top_df[sample_order,],col=list(Sex=assign_Sex,Smoking=assign_Smoking,Age=assign_Age), PD_L1_expression= anno_points(anno$PD_L1_expression,axis = TRUE),Mutation_counts = anno_barplot(anno$Mutation_counts,axis = TRUE),height=unit(c(2,2,2,2,2), "cm"),show_annotation_name=TRUE,gp = gpar(col = "white"),annotation_name_gp = gpar(fontsize=4.5))
bot_ha = HeatmapAnnotation(df=top_df[sample_order,],col=list(Efficacy=group.assign,Sex=assign_Sex,Smoking=assign_Smoking,Age=assign_Age),height=unit(c(1,1,1,1),"cm"),show_annotation_name=TRUE,gp = gpar(col = "white"),annotation_name_gp = gpar(fontsize=10))

#top_ha = HeatmapAnnotation(df=top_df_2[sample_order,],PD_L1_expression= anno_points(anno$PD_L1_expression,axis = TRUE),height = unit(2,"cm"),show_annotation_name=TRUE,gp = gpar(col = "white"),annotation_name_gp = gpar(fontsize=4.5))



#top_ha= top_ha_1+top_ha_2




#match2<-match(rownames(drawdata),pathway$gene,nomatch=0)
#rowcol<-setNames(col_rw,unique(newdata$KEGG))
#df_row<-pathway[match2,"p",drop=FALSE]
#ha_r<-rowAnnotation(df=df_row,col=list(pathway=rowcol), width = unit(0.5, "cm"))

#p<-oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],show_column_names=TRUE,column_order=NULL,axis_gp = gpar(fontsize = 4),row_names_gp= gpar(fontsize=7),column_names_gp= gpar(fontsize=7),
#    alter_fun = alter_fun_my, col = col,heatmap_legend_param = list(title = "Mutation type"))
out_pdf<-paste(opt$prefix_out,".complexheatmap.pdf",sep="")
out_png<-paste(opt$prefix_out,".complexheatmap.png",sep="")
distance 
number
if (distance>120){dis=30}
if (distance>70 && distance <= 120){dis=15}
if (distance>60 && distance <= 70){dis=30}
if (distance>50 && distance <= 60){dis=27}
if (distance>40 && distance <= 50){dis=24}
if (distance>30 && distance <= 40){dis=21}
if (distance>20 && distance <= 30){dis=18}
if (distance>10 && distance <= 20){dis=15}
if (distance>5 &&  distance <= 10){dis=12}
if (distance <= 5){dis=7}
if (number >1000){
        pdf(out_pdf,width=dis,height=50,onefile = FALSE)}
if (number >= 690 && number < 1000){
        pdf(out_pdf,width=dis,height=45,onefile = FALSE)}
if (number >= 670 && number < 690){
        pdf(out_pdf,width=dis,height=40,onefile = FALSE)}
if (number >= 650 && number < 670){
        pdf(out_pdf,width=dis,height=39,onefile = FALSE)}
if (number >= 630 && number < 650){
        pdf(out_pdf,width=dis,height=38,onefile = FALSE)}
if (number >= 610 && number < 630){
        pdf(out_pdf,width=dis,height=37,onefile = FALSE)}
if (number >= 590 && number < 610){
        pdf(out_pdf,width=dis,height=36,onefile = FALSE)}
if (number >= 570 && number < 590){
        pdf(out_pdf,width=dis,height=35,onefile = FALSE)}
if (number >= 550 && number < 570){
        pdf(out_pdf,width=dis,height=34,onefile = FALSE)}
if (number >= 530 && number < 550){
        pdf(out_pdf,width=dis,height=33,onefile = FALSE)}
if (number >= 510 && number < 530){
        pdf(out_pdf,width=dis,height=32,onefile = FALSE)}
if (number >= 490 && number < 510){
        pdf(out_pdf,width=dis,height=31,onefile = FALSE)}
if (number >= 490 && number < 510){
        pdf(out_pdf,width=dis,height=30,onefile = FALSE)}
if (number >= 470 && number < 490){
        pdf(out_pdf,width=dis,height=29,onefile = FALSE)}
if (number >= 450 && number < 470){
        pdf(out_pdf,width=dis,height=28,onefile = FALSE)}
if (number >= 430 && number < 450){
        pdf(out_pdf,width=dis,height=27,onefile = FALSE)}
if (number >= 410 && number < 430){
        pdf(out_pdf,width=dis,height=26,onefile = FALSE)}
if (number >= 390 && number < 410){
        pdf(out_pdf,width=dis,height=25,onefile = FALSE)}
if (number >= 370 && number < 390){
        pdf(out_pdf,width=dis,height=24,onefile = FALSE)}
if (number >= 350 && number < 370){
        pdf(out_pdf,width=dis,height=23,onefile = FALSE)}
if (number >= 330 && number < 350){
        pdf(out_pdf,width=dis,height=22,onefile = FALSE)}
if (number >= 310 && number < 330){
        pdf(out_pdf,width=dis,height=21,onefile = FALSE)}
if (number >= 290 && number < 310){
        pdf(out_pdf,width=dis,height=20,onefile = FALSE)}
if (number >= 270 && number < 290){
        pdf(out_pdf,width=dis,height=19,onefile = FALSE)}
if (number >= 250 && number < 270){
        pdf(out_pdf,width=dis,height=18,onefile = FALSE)}
if (number >= 230 && number < 250){
        pdf(out_pdf,width=dis,height=17,onefile = FALSE)}
if (number >= 210 && number < 230){
        pdf(out_pdf,width=dis,height=16,onefile = FALSE)}
if (number >= 190 && number < 210){
        pdf(out_pdf,width=dis,height=15,onefile = FALSE)}
if (number >= 170 && number < 190){
        pdf(out_pdf,width=dis,height=14,onefile = FALSE)}
if (number >= 160 && number < 170){
        pdf(out_pdf,width=dis,height=13,onefile = FALSE)}
if (number >= 150 && number < 160){
	pdf(out_pdf,width=dis,height=12,onefile = FALSE)}
if (number >= 140 && number < 150){
        pdf(out_pdf,width=dis,height=11.5,onefile = FALSE)}
if (number >= 130 && number < 140){
        pdf(out_pdf,width=dis,height=11,onefile = FALSE)}
if (number >= 120 && number < 130){
        pdf(out_pdf,width=dis,height=15,onefile = FALSE)}
if (number >= 110 && number < 120){
        pdf(out_pdf,width=dis,height=10.5,onefile = FALSE)}
if (number >= 100 && number < 110){
        pdf(out_pdf,width=dis,height=10,onefile = FALSE)}
if (number >= 90 && number < 100){
        pdf(out_pdf,width=dis,height=9.5,onefile = FALSE)}
if (number >= 80 && number < 90){
        pdf(out_pdf,width=dis,height=9,onefile = FALSE)}
if (number >= 70 && number < 80){
        pdf(out_pdf,width=dis,height=8.5,onefile = FALSE)}
if (number >= 60 && number < 70){
        pdf(out_pdf,width=dis,height=8.0,onefile = FALSE)}
if (number >= 50 && number < 60){
        pdf(out_pdf,width=dis,height=7.5,onefile = FALSE)}
if (number >= 40 && number < 50){
        pdf(out_pdf,width=dis,height=7,onefile = FALSE)}
if (number >= 30 && number < 40){
        pdf(out_pdf,width=dis,height=6,onefile = FALSE)}
if (number >= 20 && number < 30){
        pdf(out_pdf,width=dis,height=7,onefile = FALSE)}
if (number >= 10 && number < 20){
        pdf(out_pdf,width=dis,height=4,onefile = FALSE)}
if (number >= 5 &&  number <10){
        pdf(out_pdf,width=dis,height=3,onefile = FALSE)}
if (number <5){
        pdf(out_pdf,width=dis,height=2.7,onefile = FALSE)}


if (number <= 10){
    pr=8
    pc=5.5}
if (number > 10 && number <= 30){
    pr=7.5
    pc=6}
if (number > 30 && number <= 50){
    pr=7
    pc=6.5}
if (number > 50 && number <= 70){
    pr=6.5
    pc=6.5}
if (number > 70 && number <= 90){
    pr=6.5
    pc=6}
if (number > 90 && number <= 110){
    pr=6
    pc=7.5}
if (number > 110 && number <= 210){
    pr=9
    pc=8}
if (number > 210 && number <= 310){
    pr=5
    pc=8.5}
if (number > 310 && number <= 410){
    pr=4.5
    pc=9}
if (number > 410 && number <= 510 ){
    pr=4
    pc=9.5}
if (number > 510 && number <= 610){
    pr=3.5
    pc=10}
if (number > 610 ){
    pr=3
    pc=12}

p<-oncoPrint(drawdata[,sample_order], get_type = function(x) strsplit(x, ";")[[1]],border=TRUE,show_column_names=FALSE,column_title = "Clinical_A Heatmap",column_names_gp = gpar(fontsize=pc),pct_gp = gpar(fontsize=pr),row_names_gp = gpar(fontsize=pr),alter_fun = alter_fun_my, col = col,heatmap_legend_param = list(title = "Mutation type",bottom_annotation_height=1),remove_empty_columns=FALSE)
#p=p+theme(text=element_text(family="Arial"))+geom_text(family="Arial")


#ord <- column_order(p)
#names(ord)[2] <- "matrix"
#or <- ord$matrix
#or <- ord[[1]]
#re_order <- data.frame(or_freq = colnames(drawdata)[or], Group = NA, stringsAsFactors = F)
#for (i in 1:nrow(re_order)){
#  re_order$Group[i] <- as.character(group_data$Efficacy[which(colnames(drawdata) == re_order$or_freq[i])])
#}
#print(group_data)
#matchid<-match(re_order$or_freq,group_data$Tumor_Sample_Barcode,nomatch=0)
#re_order$Group<-group_data[matchid,"Efficacy"]
#re_order <- rbind(subset(re_order, Group == "NR"),
#                  subset(re_order, Group == "R"),
#                  subset(re_order, is.na(Group)))
#re_order <- re_order[order(re_order$Group, decreasing = F),]
#dim(re_order)
#stop()

#print(pathway)
#ord <- row_order(p)
#or <- ord[[1]]
#row.re_order <- data.frame(or_freq = rownames(drawdata)[or], Group = NA, stringsAsFactors = F)
#for (i in 1:nrow(row.re_order)){
#  row.re_order$Group[i] <- as.character(pathway$p[which(pathway$gene == row.re_order$or_freq[i])])
#}
#row.re_order <- row.re_order[order(row.re_order$Group, decreasing = F),]

p<-oncoPrint( drawdata[,sample_order], get_type = function(x) strsplit(x, ";")[[1]],border=TRUE,show_column_names=FALSE,column_order=sample_order,column_title = "Clinical_A Heatmap",column_names_gp = gpar(fontsize=pc),pct_gp = gpar(fontsize=8),row_names_gp = gpar(fontsize=8,fontface="italic"),alter_fun = alter_fun_my, col = col,heatmap_legend_param = list(title = "Mutation type",bottom_annotation_height=1),remove_empty_columns=FALSE,top_annotation = HeatmapAnnotation(Mutation_counts = anno_barplot(top_df_2$Mutation_counts,axis = TRUE),PD_L1_expression= anno_points(top_df_2$PD_L1_expression,axis = TRUE),height = unit(c(3,3),"cm"),show_annotation_name = TRUE,annotation_name_gp = gpar(fontsize=10)),bottom_annotation=bot_ha)

#p<-oncoPrint(drawdata[,sample_order], get_type = function(x) strsplit(x, ";")[[1]],show_column_names=TRUE,column_order=sample_order,column_title = "CNV Heatmap",axis_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize=pc),pct_gp = gpar(fontsize=pr),row_names_gp = gpar(fontsize=pr),alter_fun = alter_fun_my, col = col,heatmap_legend_param = list(title = "Mutation type"),remove_empty_columns=FALSE,show_row_barplot = FALSE,top_annotation = HeatmapAnnotation(df=top_df_2,gp = gpar(fontsize=3,col = "white"),annotation_name_gp = gpar(fontsize = 6),OS = anno_points(anno$OS,axis = TRUE),height = unit(2,"cm"),show_annotation_name  = TRUE),bottom_annotation=bot_ha)
#p<-oncoPrint(drawdata[,sample_order], get_type = function(x) strsplit(x, ";")[[1]],show_column_names=TRUE,column_order=sample_order,column_title = "Chinese Heatmap",axis_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize=pc),pct_gp = gpar(fontsize=pr),row_names_gp = gpar(fontsize=pr),alter_fun = alter_fun_my, col = col,heatmap_legend_param = list(title = "Mutation type"),remove_empty_columns=FALSE,show_row_barplot = FALSE,top_annotation = HeatmapAnnotation(TMB = anno_barplot(top_df_2$TMB,axis = TRUE),height = unit(2,"cm"),show_annotation_name  = TRUE),bottom_annotation=bot_ha)



formatp<-read.table(opt$k,header=T,row.names=1,sep="\t")
#matchp<-match(rownames(drawdata),rownames(formatp),nomatch=0)
matchp<-match(rownames(drawdata),rownames(formatp),nomatch=0)
newp<-formatp[matchp,]
#rownames(drawdata)
#rownames(newp)
finalp<-newp[,which(colSums(newp) > 0),drop=FALSE]
ha_cn = HeatmapAnnotation(cn = anno_text(colnames(finalp), rot = 90, just = "left", offset = unit(0, "npc") + unit(5, "mm"), gp = gpar(fontsize = 7)), annotation_height = unit(3, "cm"))
hk=Heatmap(finalp, col = c("0" = "white", "1" = "purple"), rect_gp = gpar(col = "grey"), show_row_names = FALSE, cluster_columns = TRUE,show_column_dend = FALSE, bottom_annotation = ha_cn, show_column_names = FALSE,show_heatmap_legend = FALSE, width = unit(4, "cm"), column_title = "KEGG pathway")

#pdf(out_pdf,width=12,height=8)
draw(p+hk,merge_legend = TRUE,heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
q()

#}else{
#out_pdf<-paste(opt$prefix_out,".complexheatmap.pdf",sep="")
#pdf(out_pdf,width=12,height=8)
#plot(1:5, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
#dev.off()
#q()
#}
#if (number < 50){
#draw(p,heatmap_legend_side = "bottom")
#}else{
#pdf(out_pdf,height=20,width = 50)
#draw(p)
#dev.off()
#png(out_png,type="cairo-png",width = 1800,height=1600)
#draw(p)
#dev.off()
#q()

