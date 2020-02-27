#Information ####
#This code aims to show and summerize the results obtained from deepmg
# Date: 08/02/2020
#step 1: ini####
#set path 
#updated: 17g, 07/02/2020
#path_r = "/data/projects/deepmg/analyses/analyze_res/"
if (Sys.info()['sysname'] == "Darwin"){
  path_machine = "/Users/dcongtinh/gene-abundance/experiment/results/excute_time"
  #path_machine = "/Users/dcongtinh/gene-abundance/experiment/results/fc_model/qtf_pc576_10_fillseqf_nb10_auy_gray"
}else{
  path_machine = "/data/projects/deepmg/"
} 
print(path_machine)
#path_r =paste0(path_machine,"sci_rep/read_results/")
path_r =path_machine
#path_r = "/Users/hainguyen/Documents/nthai/PhD/workspace/deepMG_tf/utils/read_results"

#Libraries
library(ggplot2)
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# Group to compare t-test
#finding the best then find the significant differene among others

ttest_group <- function(res,
                        v_groups = "dataset",
                        v_comparison =  "del0",
                        v_value1 = "train_acc",
                        v_value2 = "val_acc"
)
{
  
  
  mean_values <- aggregate(x = res, 
                           by = list(v_groups = res[,v_groups], v_comparison = res[,v_comparison]), 
                           FUN = "mean")
  
  
  #compute t-test pair each dataset
  
  
  list_group = unique(mean_values$v_groups)
  
  
  pvalue_v_value = list()
  j=1
  for( j in c(1:length(list_group))){
    
    ### for v_value 1  ###################
    #compute t-test in training
    this_group_values <- mean_values[mean_values$v_groups==list_group[j], c("v_comparison",v_value1)]
    #get the best representations then compare t-test to others
    #the name of the best color_img
    this_group_name_best = this_group_values[this_group_values[,v_value1] == max(this_group_values[,v_value1]),"v_comparison" ]
    
    
    #the vectors to store pvalue-ttest_score and significant results-signi_res
    ttest_score = vector()
    signi_res = vector()
    #get the names of others than the best color_img
    if(length(this_group_name_best)>1)
    {
      others <- this_group_values[!(this_group_values$v_comparison %in% this_group_name_best),"v_comparison"] #i=1
    }else{
      others <- this_group_values[this_group_values$v_comparison != this_group_name_best,"v_comparison"] #i=1  
    }
    
    i = 2
    signi_res = vector()
    #compare t-test between the best and the other (each by each)
    if (length(others)>0){
      
      for( i in c(1:length(others))){
        #print (others[i])
        #get train_acc of the best
        t_best = res[res[,v_comparison] == this_group_name_best & res[,v_groups]==list_group[j],v_value1]
        #get train_acc of the one other than the best
        t_other = res[res[,v_comparison] ==  others[i] & res[,v_groups]==list_group[j],v_value1]
        #compute t-test and get p-value
        ttest_score_temp = t.test(t_best,t_other)
        ttest_score[i]  = ttest_score_temp$p.value
        #print (names(ttest_score))
        #print(ttest_score$p.value)
        print(paste0(j,"_",i,list_group[j],"_",others[i],v_value1))
        if (ttest_score[i] < 0.05){   #if significant
          print (paste0(this_group_name_best," is significant improvement compared to ",others[i], " with p-value=",ttest_score[i]))
          signi_res =rbind(signi_res, "significant")
        }else{  #if not significant
          print (paste0(this_group_name_best," not is significant difference compared to ",others[i], " , p-value=",ttest_score[i]))
          signi_res =rbind(signi_res, "no")
        }
      }
      
      pvalue_v1 = data.frame(g1=others,ttest_s = ttest_score,signi_note=signi_res)
      pvalue_v1 = rbind(pvalue_v1, data.frame(g1=this_group_name_best,ttest_s=0,signi_note="best"))
      pvalue_v1 = data.frame(pvalue_v1, groups = list_group[j])
      pvalue_v_value =rbind(pvalue_v_value,cbind(pvalue_v1,type=v_value1))
      
      
    }
    
    
    
    ### for v_value 2  ###################
    #compute t-test in training
    this_group_values <- mean_values[mean_values$v_groups==list_group[j], c("v_comparison",v_value2)]
    #get the best representations then compare t-test to others
    #the name of the best color_img
    this_group_name_best = this_group_values[this_group_values[,v_value2] == max(this_group_values[,v_value2]),"v_comparison" ]
    
    print("begin v2")
    #the vectors to store pvalue-ttest_score and significant results-signi_res
    ttest_score = vector()
    signi_res = vector()
    #get the names of others than the best color_img
    #others <- this_group_values[this_group_values$v_comparison != this_group_name_best,"v_comparison"] #i=1
    if(length(this_group_name_best)>1)
    {
      others <- this_group_values[!(this_group_values$v_comparison %in% this_group_name_best),"v_comparison"] #i=1
    }else{
      others <- this_group_values[this_group_values$v_comparison != this_group_name_best,"v_comparison"] #i=1  
    }
    
    i = 1
    signi_res = vector()
    #compare t-test between the best and the other (each by each)
    if (length(others)>0){
      for( i in c(1:length(others))){
        #print (others[i])
        #get train_acc of the best
        t_best = res[res[,v_comparison] == this_group_name_best & res[,v_groups]==list_group[j],v_value2]
        #get train_acc of the one other than the best
        t_other = res[res[,v_comparison] ==  others[i] & res[,v_groups]==list_group[j],v_value2]
        #compute t-test and get p-value
        ttest_score_temp = t.test(t_best,t_other)
        ttest_score[i]  = ttest_score_temp$p.value
        #print (names(ttest_score))
        #print(ttest_score$p.value)
        print(paste0(j,"_",i,list_group[j],"_",others[i],v_value2))
        if (ttest_score[i] < 0.05){   #if significant
          print (paste0(this_group_name_best," is significant improvement compared to ",others[i], " with p-value=",ttest_score[i]))
          signi_res =rbind(signi_res, "significant")
        }else{  #if not significant
          print (paste0(this_group_name_best," not is significant difference compared to ",others[i], " , p-value=",ttest_score[i]))
          signi_res =rbind(signi_res, "no")
        }
      }
      
      pvalue_v2 = data.frame(g1=others,ttest_s = ttest_score,signi_note=signi_res)
      pvalue_v2 = rbind(pvalue_v2, data.frame(g1=this_group_name_best,ttest_s=0,signi_note="best"))
      pvalue_v2 = data.frame(pvalue_v2, groups = list_group[j])
      pvalue_v_value =rbind(pvalue_v_value,cbind(pvalue_v2,type=v_value2))
    }
    
    
  }
  return(pvalue_v_value)
}
#step 2a: read data:  cross-validation: cv####
setwd(path_r)
#read the results
#rSum = read.table("results_sum.txt", sep = "\t", header = T)
# r1 = read.table("results_ext.txt", sep = "\t", header = T) #collect internal validation results
r1 = list()
r2 = read.table("results_folds.txt", sep = "\t", header = T) #collect external validation results
#all_results = rbind(cbind(r2,external_vali="no"),cbind(r1,external_vali = "yes"))
all_results = rbind(r2,r1)

all_results = cbind(all_results, date_time = substr(all_results$filename, regexpr('_\\d{8}_\\d{6}', all_results$filename)+1, 
       regexpr('_\\d{8}_\\d{6}', all_results$filename)+15))
unique(all_results$filename)
dim(all_results)
colnames(all_results)
head(all_results)

#declare vectors to filter results with some conditions
representation = vector()
color_img = vector()
model = vector()
bin = vector()
dataset = vector()
taxa = vector()
whole = vector()
estop = vector()
padding = vector()
type_evaluate = vector()
resolution = vector()

representation = rep(NA,nrow(all_results))
representation[regexpr('rnd1_fill', all_results[,"filename"]) >= 1 & is.na(representation)] <- "fill_rnd"
representation[regexpr('_tsne_', all_results[,"filename"]) >= 1 & is.na(representation)] <- "tsne"
representation[regexpr('_nmf_', all_results[,"filename"]) >= 1 & is.na(representation)] <- "nmf"
representation[regexpr('_isomap_', all_results[,"filename"]) >= 1 & is.na(representation)] <- "isomap"
representation[regexpr('_mds_', all_results[,"filename"]) >= 1 & is.na(representation)] <- "mds"
representation[regexpr('_pca_', all_results[,"filename"]) >= 1 & is.na(representation)] <- "pca"
representation[regexpr('_fill', all_results[,"filename"]) >= 1 & is.na(representation)] <- "fill_phy"
representation[regexpr('_raw', all_results[,"filename"]) >= 1 & is.na(representation)] <- "raw-1D"
representation[regexpr('_lle_', all_results[,"filename"]) >= 1 & is.na(representation)] <- "lle"
unique(representation)

color_img =  rep(NA,nrow(all_results))
color_img[regexpr('custom', all_results[,"filename"]) >= 1 & is.na(color_img)] <- "color"
color_img[regexpr('gray', all_results[,"filename"]) >= 1 & is.na(color_img)] <- "gray"



model =  rep(NA,nrow(all_results))
model[regexpr('_vgg_', all_results[,"filename"]) >= 1 & is.na(model)] <- "vgg"
model[regexpr('l1f16dfc0.0p1dcnn', all_results[,"filename"]) >= 1 & is.na(model)] <- "cnn_l1f16"
model[regexpr('l1f64dfc0.0p1dcnn', all_results[,"filename"]) >= 1 & is.na(model)] <- "cnn_l1f64"
model[regexpr('l1f128dfc0.0p1dcnn', all_results[,"filename"]) >= 1 & is.na(model)] <- "cnn_l1f128"
model[regexpr('l2f16dfc0.0p1dcnn', all_results[,"filename"]) >= 1 & is.na(model)] <- "cnn_l2f16"
model[regexpr('l2f64dfc0.0p1dcnn', all_results[,"filename"]) >= 1 & is.na(model)] <- "cnn_l2f64"
model[regexpr('l2f128dfc0.0p1dcnn', all_results[,"filename"]) >= 1 & is.na(model)] <-"cnn_l2f128"
model[regexpr('l1f32dfc0.0p1dcnn', all_results[,"filename"]) >= 1 & is.na(model)] <-"cnn_l1f32"
model[regexpr('cn1d', all_results[,"filename"]) >= 1 & is.na(model)] <- "cn1d"
model[regexpr('fc', all_results[,"filename"]) >= 1 & is.na(model)] <- "fc"
model[regexpr('gbc', all_results[,"filename"]) >= 1 & is.na(model)] <- "gbc"
model[regexpr('svm', all_results[,"filename"]) >= 1 & is.na(model)] <- "svm"
model[regexpr('rf', all_results[,"filename"]) >= 1 & is.na(model)] <- "rf"

bin =  rep(NA,nrow(all_results))
bin[regexpr('spb', all_results[,"filename"]) >= 1 & is.na(bin)] <- "spb"
bin[regexpr('eqw', all_results[,"filename"]) >= 1 & is.na(bin)] <- "eqw"

scale_mode =  rep(NA,nrow(all_results))
scale_mode[regexpr('qtf', all_results[,"filename"]) >= 1 & is.na(scale_mode)] <- "qtf"
scale_mode[regexpr('mms', all_results[,"filename"]) >= 1 & is.na(scale_mode)] <- "Min_Max"


redim =  rep(NA,nrow(all_results))
redim[regexpr('ridge', all_results[,"filename"]) >= 1 & is.na(redim)] <- "ridge"
redim[regexpr('rigde', all_results[,"filename"]) >= 1 & is.na(redim)] <- "ridge"
redim[regexpr('lasso', all_results[,"filename"]) >= 1 & is.na(redim)] <- "lasso"
redim[regexpr('vf', all_results[,"filename"]) >= 1 & is.na(redim)] <- "lvf"
redim[regexpr('raw', all_results[,"filename"]) >= 1 & is.na(redim)] <- "Raw"
redim[regexpr('pc', all_results[,"filename"]) >= 1 & is.na(redim)] <- "Perceptron Weight Based Filter"
redim[is.na(redim)] <- "none"


dataset = rep(NA,nrow(all_results))
dataset[regexpr('wt2dgene', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "wt2"
dataset[regexpr('t2dgene', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "t2d"
dataset[regexpr('cirgene', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "cir"
dataset[regexpr('colgene', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "col"
dataset[regexpr('ibdgene', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "ibd"
dataset[regexpr('obegene', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "obe"
# dataset[regexpr('qinn.stage1', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "qinn.stage1"
# dataset[regexpr('karlsson', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "karlsson"
# dataset[regexpr('lechat', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "lechat"
# dataset[regexpr('qinj', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "qinj"
# dataset[regexpr('nielsen', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "nielsen"
# dataset[regexpr('zellerg', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "zellerg"
# dataset[regexpr('db_obesity', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "obesity"
# dataset[regexpr('db_t2dw', all_results[,"filename"]) >= 1 & is.na(dataset)] <- "t2dw"


taxa = rep(NA,nrow(all_results))
taxa[regexpr('freq.bug_all', all_results[,"filename"]) >= 1 & is.na(taxa)] <- "freq.bug_all"
taxa[regexpr('freq.bug_class', all_results[,"filename"]) >= 1 & is.na(taxa)] <- "freq.bug_class"
taxa[regexpr('freq.bug_family', all_results[,"filename"]) >= 1 & is.na(taxa)] <- "freq.bug_family"
taxa[regexpr('freq.bug_genus', all_results[,"filename"]) >= 1 & is.na(taxa)] <- "freq.bug_genus"
taxa[regexpr('freq.bug_order', all_results[,"filename"]) >= 1 & is.na(taxa)] <- "freq.bug_order"
taxa[regexpr('freq.bug_phylum', all_results[,"filename"]) >= 1 & is.na(taxa)] <- "freq.bug_phylum"
taxa[regexpr('freq.bug_species', all_results[,"filename"]) >= 1 & is.na(taxa)] <- "freq.bug_species"
taxa[regexpr('pathway', all_results[,"filename"]) >= 1 & is.na(taxa)] <- "pathway"
taxa[regexpr('marker', all_results[,"filename"]) >= 1 & is.na(taxa)] <- "marker"

whole = rep(NA,nrow(all_results))
whole[regexpr('whole', all_results[,"filename"]) >= 1 & is.na(whole)] <- "yes"
whole[regexpr('whole', all_results[,"filename"]) < 1 & is.na(whole)] <- "no"

estop = rep(NA,nrow(all_results))
estop[regexpr('estopc500', all_results[,"filename"]) >= 1 & is.na(estop)] <- "Without-early-Stopping"
estop[regexpr('estopc50', all_results[,"filename"]) >= 1 & is.na(estop)] <- "estop50"
estop[regexpr('estopc5', all_results[,"filename"]) >= 1 & is.na(estop)] <- "estop5"
estop[regexpr('estopc10', all_results[,"filename"]) >= 1 & is.na(estop)] <- "estop10"
estop[regexpr('estopc1_', all_results[,"filename"]) >= 1 & is.na(estop)] <- "estop1"
estop[regexpr('estopc-1', all_results[,"filename"]) >= 1 & is.na(estop)] <- "Without-early-Stopping"

padding = rep(NA,nrow(all_results))
padding[regexpr('pady', all_results[,"filename"]) >= 1 & is.na(padding)] <- "yes"
padding[is.na(padding)] <- "no"

type_evaluate = rep(NA,nrow(all_results))
type_evaluate[regexpr('_sum', all_results[,"filename"]) >= 1 & is.na(type_evaluate)] <- "internal"
type_evaluate[regexpr('_ext.', all_results[,"filename"]) >= 1 & is.na(type_evaluate)] <- "external"

resolution = rep(NA,nrow(all_results))
resolution[regexpr('_r48', all_results[,"filename"]) >= 1 & is.na(resolution)] <- "48x48"
resolution[regexpr('_r96', all_results[,"filename"]) >= 1 & is.na(resolution)] <- "96x96"
resolution[regexpr('_r0', all_results[,"filename"]) >= 1 & is.na(resolution)] <- "autofit"

autofit_data =  rep(NA,nrow(all_results))
autofit_data[regexpr('_auy_', all_results[,"filename"]) >= 1 & is.na(autofit_data)] <- "scale_range_min_max"
autofit_data[regexpr('_aun_', all_results[,"filename"]) >= 1 & is.na(autofit_data)] <- "scale_range_0_1"
autofit_data[is.na(autofit_data)] <- "NA"

del0 = rep(NA,nrow(all_results))
del0[regexpr('_del0_', all_results[,"filename"]) >= 1 & is.na(del0)] <- "yes"
del0[is.na(del0)] <- "no"

#check the length to sure vectors filled completely
length(representation)
length(color_img)
length(model)
length(bin)
length(dataset)
length(taxa)
length(resolution)

bin_scale = paste0(bin,"_",scale_mode)
#combine these vectors into the all_results_combined, then use all_results_combined for the analyse
all_results_combined = list()
all_results_combined <- data.frame(all_results,dataset,representation,color_img,model,bin,taxa,estop,padding,type_evaluate,
                                   whole,resolution, 
                                   autofit_data,scale_mode, bin_scale,del0,
                                   redim)
length( unique(all_results_combined$filename))
#View(all_results_combined)
#		tinh:wt2d####
#set condition to filter the results:
#set condition:
res_external1=list()
res_internal1=list()
v_dataset = "wt2"
v_bin = "eqf"
v_bin_scale = "eqf_qtf"
v_color_img = "gray"
v_model = "fc"
v_taxa = "freq.bug_species"
v_estop = "estop5"
v_padding = "no"
v_autofit_data = "scale_range_min_max"
#name for the chart
#title_chart <- paste("various representations: ",v_dataset,
#                     v_bin,
#                    v_color_img ,
#                     v_model,
#                     v_taxa ,
#                     v_estop )
#colnames(all_results_combined)
#filter results in internal validation
res_internal1 = all_results_combined[ all_results_combined$type_evaluate == "internal" & all_results_combined$whole == "no" 
  
  #all_results_combined$dataset== v_dataset
                                      #   all_results_combined$color_img== v_color_img 
                                         
                                     &  all_results_combined$model== v_model
                                         
                                    # &  all_results_combined$bin_scale==v_bin_scale 
                                       #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                        # all_results_combined$padding == v_padding & all_results_combined$autofit_data== v_autofit_data &
                                        
                                       ,]
table(res_internal1[,c("bin_scale","dataset")])
unique(res_internal1$filename)
unique(res_internal1$representation)
res_internal_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "no" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "internal"
                                        ,]
res_internal1 = rbind(res_internal1,res_internal_raw)
#filter results in external validation
res_external1 = all_results_combined[all_results_combined$type_evaluate == "external" & all_results_combined$whole == "yes" 
                                     
                                     #all_results_combined$dataset== v_dataset
                                     #   all_results_combined$color_img== v_color_img 
                                     
                                     &  all_results_combined$model== v_model
                                     
                                     # &  all_results_combined$bin_scale==v_bin_scale 
                                     #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                     # all_results_combined$padding == v_padding & all_results_combined$autofit_data== v_autofit_data &
                                     
                                     ,]
unique(res_external1$filename)
unique(res_internal1$filename)
res_external_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "yes" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "external"
                                        ,]
res_external1 = rbind(res_external1,res_external_raw)
v_groups = "dataset"
v_comparison = "model"
mean_values <- aggregate(x = res_internal1, 
                         by = list(v_groups = res_internal1[,v_groups], v_comparison = res_internal1[,v_comparison]), 
                         FUN = "mean")

ratio_majority_class = 0.5414365
#training
#position y of significant results to show in the charts
y_signifi_position = 1.01
#get min, max for setting the charts
min_value = min(min(res_internal1$train_acc),min(res_internal1$val_acc))
max_value = max(max(res_internal1$train_acc),max(res_internal1$val_acc),y_signifi_position)
#the chart of training performance
res_internal1 <- res_internal1[!is.na(res_internal1$dataset),]
names(mean_values)[1] <- 'dataset'
p1 = 
  ggplot(data=res_internal1, aes(x=dataset, y=train_acc, fill=dataset)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))+
  #ggtitle(paste("Training on ",title_chart))  +
  geom_text(data =  mean_values[,c("dataset","train_acc")], aes(label = round(train_acc,3), y = train_acc - 0.055), size=5) +
  
 # geom_text(data = means_train, aes(label = train_acc, y = train_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_train,shape =8, size =3,
   #          aes(x=representation, y=y_signifi_position)) +
  theme(legend.position="none") +
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
 # geom_text(aes(2,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  ylim(min_value,max_value)

#testing and external, significant results


names(mean_values)[1] <- 'dataset'
#mean_values_val <- mean_values[,c("dataset","val_acc")]
p2 = 
  ggplot(data=res_internal1, aes(x=dataset, y=val_acc, fill=dataset))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))+
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
 # ggtitle(paste("Testing on ",title_chart))+ 
  geom_boxplot() +
 # stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3)   +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) +
  #ggtitle(paste("Test on ",title_chart))  +
  geom_text(data =  mean_values[,c("dataset","val_acc")], aes(label = round(val_acc,3), y = val_acc - 0.053), size=5) +
  #external+mean
  #geom_point(data=val_int_ext, size = 3,fill="black",
   #          aes(x=representation, y=value, shape = type_evaluate ))+
  #scale_shape_manual(values=c(4, 24))+
  #geom_text(data = means_val, aes(label = val_acc, y = val_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_val,shape =8, size =3,
  #           aes(x=representation, y=y_signifi_position)) +
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
  #geom_text(aes(1,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  ylim(min_value,max_value)
multiplot(p1, p2,cols=2)
table(res_internal1[,c("model","dataset")])
dim(res_internal1)
mean(res_internal1$val_acc)
table(res_internal1[,"model"])
table(res_external1[,"representation"])

#		cn1d:wt2d####
#set condition to filter the results:
#set condition:
res_external1=list()
res_internal1=list()
v_dataset = "wt2"
v_bin = "eqw"
v_bin_scale = "eqw_qtf"
v_color_img = "gray"
v_model = "cn1d"
v_taxa = "freq.bug_species"
v_estop = "estop5"
v_padding = "no"
v_autofit_data = "scale_range_min_max"
#name for the chart
#title_chart <- paste("various representations: ",v_dataset,
#                     v_bin,
#                    v_color_img ,
#                     v_model,
#                     v_taxa ,
#                     v_estop )
#colnames(all_results_combined)
#filter results in internal validation
res_internal1 = all_results_combined[ all_results_combined$type_evaluate == "internal" & all_results_combined$whole == "no" 
                                      
                                      #all_results_combined$dataset== v_dataset
                                      #   all_results_combined$color_img== v_color_img 
                                      
                                      &  all_results_combined$model== v_model
                                      
                                      # &  all_results_combined$bin_scale==v_bin_scale 
                                      #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                      # all_results_combined$padding == v_padding & all_results_combined$autofit_data== v_autofit_data &
                                      
                                      ,]
table(res_internal1[,c("bin_scale","dataset")])
unique(res_internal1$filename)
unique(res_internal1$representation)
res_internal_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "no" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "internal"
                                        ,]
res_internal1 = rbind(res_internal1,res_internal_raw)
#filter results in external validation
res_external1 = all_results_combined[all_results_combined$type_evaluate == "external" & all_results_combined$whole == "yes" 
                                     
                                     #all_results_combined$dataset== v_dataset
                                     #   all_results_combined$color_img== v_color_img 
                                     
                                     &  all_results_combined$model== v_model
                                     
                                     # &  all_results_combined$bin_scale==v_bin_scale 
                                     #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                     # all_results_combined$padding == v_padding & all_results_combined$autofit_data== v_autofit_data &
                                     
                                     ,]
unique(res_external1$filename)
unique(res_internal1$filename)
res_external_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "yes" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "external"
                                        ,]
res_external1 = rbind(res_external1,res_external_raw)
v_groups = "dataset"
v_comparison = "model"
mean_values <- aggregate(x = res_internal1, 
                         by = list(v_groups = res_internal1[,v_groups], v_comparison = res_internal1[,v_comparison]), 
                         FUN = "mean")

ratio_majority_class = 0.5414365
#training
#position y of significant results to show in the charts
y_signifi_position = 1.01
#get min, max for setting the charts
min_value = min(min(res_internal1$train_acc),min(res_internal1$val_acc))
max_value = max(max(res_internal1$train_acc),max(res_internal1$val_acc),y_signifi_position)
#the chart of training performance
res_internal1 <- res_internal1[!is.na(res_internal1$dataset),]
names(mean_values)[1] <- 'dataset'
p1 = 
  ggplot(data=res_internal1, aes(x=dataset, y=train_acc, fill=dataset)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) +
  #ggtitle(paste("Training on ",title_chart))  +
  geom_text(data =  mean_values[,c("dataset","train_acc")], aes(label = train_acc, y = train_acc - 0.05)) +
  
  # geom_text(data = means_train, aes(label = train_acc, y = train_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_train,shape =8, size =3,
  #          aes(x=representation, y=y_signifi_position)) +
  theme(legend.position="none") +
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
  # geom_text(aes(2,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  ylim(min_value,max_value)

#testing and external, significant results


names(mean_values)[1] <- 'dataset'
#mean_values_val <- mean_values[,c("dataset","val_acc")]
p2 = 
  ggplot(data=res_internal1, aes(x=dataset, y=val_acc, fill=dataset))+
  # ggtitle(paste("Testing on ",title_chart))+ 
  geom_boxplot() +
  # stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3)   +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) +
  #ggtitle(paste("Test on ",title_chart))  +
  geom_text(data =  mean_values[,c("dataset","val_acc")], aes(label = val_acc, y = val_acc - 0.05)) +
  #external+mean
  #geom_point(data=val_int_ext, size = 3,fill="black",
  #          aes(x=representation, y=value, shape = type_evaluate ))+
  #scale_shape_manual(values=c(4, 24))+
  #geom_text(data = means_val, aes(label = val_acc, y = val_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_val,shape =8, size =3,
  #           aes(x=representation, y=y_signifi_position)) +
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
  #geom_text(aes(1,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  ylim(min_value,max_value)
multiplot(p1, p2,cols=2)
table(res_internal1[,c("model","dataset")])
dim(res_internal1)
mean(res_internal1$val_acc)
table(res_internal1[,"model"])
table(res_external1[,"representation"])

#		ridge:wt2d####
#set condition to filter the results:
#set condition:
res_external1=list()
res_internal1=list()
v_dataset = "wt2"
v_bin = "eqf"
v_bin_scale = "eqf_qtf"
v_color_img = "gray"
v_model = "fc"
#v_model = "fc"
v_taxa = "freq.bug_species"
v_estop = "estop5"
v_padding = "no"
v_autofit_data = "scale_range_min_max"
v_redim = "pc"
#name for the chart
#title_chart <- paste("various representations: ",v_dataset,
#                     v_bin,
#                    v_color_img ,
#                     v_model,
#                     v_taxa ,
#                     v_estop )
#colnames(all_results_combined)
#filter results in internal validation
res_internal1 = all_results_combined[ all_results_combined$type_evaluate == "internal" & all_results_combined$whole == "no" 
                                      
                                      #all_results_combined$dataset== v_dataset
                                      #   all_results_combined$color_img== v_color_img 
                                      
                                    #  &  all_results_combined$model== v_model
                                      & all_results_combined$redim == v_redim
                                       &  all_results_combined$bin_scale==v_bin_scale 
                                      #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                      # all_results_combined$padding == v_padding & all_results_combined$autofit_data== v_autofit_data &
                                      
                                      ,]
table(res_internal1[,c("redim","dataset")])
unique(res_internal1$filename)
unique(res_internal1$representation)
res_internal_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "no" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "internal"
                                        ,]
res_internal1 = rbind(res_internal1,res_internal_raw)
#filter results in external validation
res_external1 = all_results_combined[all_results_combined$type_evaluate == "external" & all_results_combined$whole == "yes" 
                                     
                                     #all_results_combined$dataset== v_dataset
                                     #   all_results_combined$color_img== v_color_img 
                                     
                                     &  all_results_combined$model== v_model
                                     
                                     # &  all_results_combined$bin_scale==v_bin_scale 
                                     #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                     # all_results_combined$padding == v_padding & all_results_combined$autofit_data== v_autofit_data &
                                     
                                     ,]
unique(res_external1$filename)
unique(res_internal1$filename)
res_external_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "yes" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "external"
                                        ,]
res_external1 = rbind(res_external1,res_external_raw)
v_groups = "dataset"
v_comparison = "model"
mean_values <- aggregate(x = res_internal1, 
                         by = list(v_groups = res_internal1[,v_groups], v_comparison = res_internal1[,v_comparison]), 
                         FUN = "mean")

ratio_majority_class = 0.5414365
#training
#position y of significant results to show in the charts
y_signifi_position = 1.01
#get min, max for setting the charts
min_value = min(min(res_internal1$train_acc),min(res_internal1$val_acc))
max_value = max(max(res_internal1$train_acc),max(res_internal1$val_acc),y_signifi_position)
#the chart of training performance
res_internal1 <- res_internal1[!is.na(res_internal1$dataset),]
names(mean_values)[1] <- 'dataset'
p1 = 
  ggplot(data=res_internal1, aes(x=dataset, y=train_acc, fill=model)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) +
  #ggtitle(paste("Training on ",title_chart))  +
  #geom_text(data =  mean_values[,c("dataset","train_acc")], aes(label = train_acc, y = train_acc - 0.05)) +
  
  # geom_text(data = means_train, aes(label = train_acc, y = train_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_train,shape =8, size =3,
  #          aes(x=representation, y=y_signifi_position)) +
  theme(legend.position="none") +
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
  # geom_text(aes(2,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  ylim(min_value,max_value)

#testing and external, significant results


names(mean_values)[1] <- 'dataset'
#mean_values_val <- mean_values[,c("dataset","val_acc")]
p2 = 
  ggplot(data=res_internal1, aes(x=dataset, y=val_acc, fill=model))+
  # ggtitle(paste("Testing on ",title_chart))+ 
  geom_boxplot() +
  # stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3)   +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) +
  #ggtitle(paste("Test on ",title_chart))  +
  #geom_text(data =  mean_values[,c("dataset","val_acc")], aes(label = val_acc, y = val_acc - 0.05)) +
  #external+mean
  #geom_point(data=val_int_ext, size = 3,fill="black",
  #          aes(x=representation, y=value, shape = type_evaluate ))+
  #scale_shape_manual(values=c(4, 24))+
  #geom_text(data = means_val, aes(label = val_acc, y = val_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_val,shape =8, size =3,
  #           aes(x=representation, y=y_signifi_position)) +
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
  #geom_text(aes(1,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  ylim(min_value,max_value)
multiplot(p1, p2,cols=2)
table(res_internal1[,c("redim","dataset")])
dim(res_internal1)
mean(res_internal1$val_acc)
table(res_internal1[,"model"])
table(res_external1[,"representation"])


ggplot(data=res_internal1, aes(x=dataset, y=time, fill=model)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3)# +
  #ggtitle(paste("Training on ",title_chart))  +
  #geom_text(data =  mean_values[,c("dataset","train_acc")], aes(label = train_acc, y = train_acc - 0.05)) +
  
  # geom_text(data = means_train, aes(label = train_acc, y = train_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_train,shape =8, size =3,
  #          aes(x=representation, y=y_signifi_position)) +
  #theme(legend.position="none") 
#+
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
  # geom_text(aes(2,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  

#		redim-fc:wt2d####
#set condition to filter the results:
#set condition:
res_external1=list()
res_internal1=list()
v_dataset = "wt2"
v_bin = "eqf"
v_bin_scale = "NA_NA"
v_color_img = "gray"
#v_model = "cnn_l1f32"
v_model = "fc"
v_taxa = "freq.bug_species"
v_estop = "estop5"
v_padding = "no"
v_autofit_data = "scale_range_min_max"
v_redim = "pc"

#name for the chart
#title_chart <- paste("various representations: ",v_dataset,
#                     v_bin,
#                    v_color_img ,
#                     v_model,
#                     v_taxa ,
#                     v_estop )
#colnames(all_results_combined)
#filter results in internal validation
res_internal1 = all_results_combined[ all_results_combined$type_evaluate == "internal" & all_results_combined$whole == "no" 
                                      
                                      #all_results_combined$dataset== v_dataset
                                      #   all_results_combined$color_img== v_color_img 
                                      
                                        &  all_results_combined$model== v_model
                                    #  & all_results_combined$redim == v_redim
                                      &  all_results_combined$bin_scale==v_bin_scale 
                                      #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                      # all_results_combined$padding == v_padding
                                      & all_results_combined$autofit_data== v_autofit_data 
                                      
                                      ,]
table(res_internal1[,c("autofit_data","dataset")])
unique(res_internal1$filename)
unique(res_internal1$representation)
res_internal_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "no" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "internal"
                                        ,]
res_internal1 = rbind(res_internal1,res_internal_raw)
#filter results in external validation
res_external1 = all_results_combined[all_results_combined$type_evaluate == "external" & all_results_combined$whole == "yes" 
                                     
                                     #all_results_combined$dataset== v_dataset
                                     #   all_results_combined$color_img== v_color_img 
                                     
                                     &  all_results_combined$model== v_model
                                     
                                     # &  all_results_combined$bin_scale==v_bin_scale 
                                     #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                     # all_results_combined$padding == v_padding & all_results_combined$autofit_data== v_autofit_data &
                                     
                                     ,]
unique(res_external1$filename)
unique(res_internal1$filename)
res_external_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "yes" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "external"
                                        ,]
res_external1 = rbind(res_external1,res_external_raw)
v_groups = "dataset"
v_comparison = "model"
mean_values <- aggregate(x = res_internal1, 
                         by = list(v_groups = res_internal1[,v_groups], v_comparison = res_internal1[,v_comparison]), 
                         FUN = "mean")

ratio_majority_class = 0.5414365
#training
#position y of significant results to show in the charts
y_signifi_position = 1.01
#get min, max for setting the charts
min_value = min(min(res_internal1$train_acc),min(res_internal1$val_acc))
max_value = max(max(res_internal1$train_acc),max(res_internal1$val_acc),y_signifi_position)
#the chart of training performance
res_internal1 <- res_internal1[!is.na(res_internal1$dataset),]
names(mean_values)[1] <- 'dataset'
p3 = 
  ggplot(data=res_internal1, aes(x=dataset, y=train_acc, fill=redim)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) + #+geom_point() +
#  stat_summary(fun.y=mean, geom="line", size=1.5,
 #              linetype="dotted")
  ggtitle("Training on FC model") + 
  #ggtitle(paste("Training on ",title_chart))  +
  #geom_text(data =  mean_values[,c("dataset","train_acc")], aes(label = train_acc, y = train_acc - 0.05)) +
  
  # geom_text(data = means_train, aes(label = train_acc, y = train_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_train,shape =8, size =3,
  #          aes(x=representation, y=y_signifi_position)) +
  theme(legend.position="none") +
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
  # geom_text(aes(2,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  ylim(min_value,max_value)

#testing and external, significant results


names(mean_values)[1] <- 'dataset'
#mean_values_val <- mean_values[,c("dataset","val_acc")]
p4 = 
  ggplot(data=res_internal1, aes(x=dataset, y=val_acc, fill=redim))+
  ggtitle("Training on FC model") + 
  # ggtitle(paste("Testing on ",title_chart))+ 
  geom_boxplot() +
  # stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3)   +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) +
  theme(legend.position="right") + labs(fill = "FS methods")+
  ylim(min_value,max_value)


multiplot(p3, p4,cols=2)
#table(res_internal1[,c("redim","dataset")])
#dim(res_internal1)
#mean(res_internal1$val_acc)
#table(res_internal1[,"model"])
#table(res_external1[,"representation"])


ggplot(data=res_internal1, aes(x=dataset, y=time, fill=model)) + geom_boxplot() +
            stat_summary(fun.y=mean, colour="black", geom="point", 
            shape=4, size=3)# +
#multiplot(p1,p2,p3, p4,cols=4)
#		redim-cnn:wt2d####
#set condition to filter the results:
#set condition:
res_external1=list()
res_internal1=list()
v_dataset = "wt2"
v_bin = "eqw"
v_bin_scale = "eqw_NA"
v_color_img = "gray"
v_model = "cnn_l1f32"
#v_model = "fc"
v_taxa = "freq.bug_species"
v_estop = "estop5"
v_padding = "no"
v_autofit_data = "scale_range_min_max"
v_redim = "ridge"

#name for the chart
#title_chart <- paste("various representations: ",v_dataset,
#                     v_bin,
#                    v_color_img ,
#                     v_model,
#                     v_taxa ,
#                     v_estop )
#colnames(all_results_combined)
#filter results in internal validation
res_internal1 = all_results_combined[ all_results_combined$type_evaluate == "internal" & all_results_combined$whole == "no" 
                                      
                                      #all_results_combined$dataset== v_dataset
                                      #   all_results_combined$color_img== v_color_img 
                                      
                                      &  all_results_combined$model== v_model
                                      #  & all_results_combined$redim == v_redim
                                      &  all_results_combined$bin_scale==v_bin_scale 
                                      #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                      # all_results_combined$padding == v_padding
                                      & all_results_combined$autofit_data== v_autofit_data 
                                      
                                      ,]
table(res_internal1[,c("autofit_data","dataset")])
unique(res_internal1$filename)
unique(res_internal1$representation)
res_internal_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "no" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "internal"
                                        ,]
res_internal1 = rbind(res_internal1,res_internal_raw)
#filter results in external validation
res_external1 = all_results_combined[all_results_combined$type_evaluate == "external" & all_results_combined$whole == "yes" 
                                     
                                     #all_results_combined$dataset== v_dataset
                                     #   all_results_combined$color_img== v_color_img 
                                     
                                     &  all_results_combined$model== v_model
                                     
                                     # &  all_results_combined$bin_scale==v_bin_scale 
                                     #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                     # all_results_combined$padding == v_padding & all_results_combined$autofit_data== v_autofit_data &
                                     
                                     ,]
unique(res_external1$filename)
unique(res_internal1$filename)
res_external_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "yes" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "external"
                                        ,]
res_external1 = rbind(res_external1,res_external_raw)
v_groups = "dataset"
v_comparison = "model"
mean_values <- aggregate(x = res_internal1, 
                         by = list(v_groups = res_internal1[,v_groups], v_comparison = res_internal1[,v_comparison]), 
                         FUN = "mean")

ratio_majority_class = 0.5414365
#training
#position y of significant results to show in the charts
y_signifi_position = 1.01
#get min, max for setting the charts
min_value = min(min(res_internal1$train_acc),min(res_internal1$val_acc))
max_value = max(max(res_internal1$train_acc),max(res_internal1$val_acc),y_signifi_position)
#the chart of training performance
res_internal1 <- res_internal1[!is.na(res_internal1$dataset),]
names(mean_values)[1] <- 'dataset'
p3 = 
  ggplot(data=res_internal1, aes(x=dataset, y=train_acc, fill=redim)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) + #+geom_point() +
  #  stat_summary(fun.y=mean, geom="line", size=1.5,
  #              linetype="dotted")
  ggtitle("Training on CNN1D model") + 
  #ggtitle(paste("Training on ",title_chart))  +
  #geom_text(data =  mean_values[,c("dataset","train_acc")], aes(label = train_acc, y = train_acc - 0.05)) +
  
  # geom_text(data = means_train, aes(label = train_acc, y = train_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_train,shape =8, size =3,
  #          aes(x=representation, y=y_signifi_position)) +
  theme(legend.position="none") +
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
  # geom_text(aes(2,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  ylim(min_value,max_value)

#testing and external, significant results


names(mean_values)[1] <- 'dataset'
#mean_values_val <- mean_values[,c("dataset","val_acc")]
p4 = 
  ggplot(data=res_internal1, aes(x=dataset, y=val_acc, fill=redim))+
  ggtitle("Training on CNN1D model") +
  # ggtitle(paste("Testing on ",title_chart))+ 
  geom_boxplot() +
  # stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3)   +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) +
  theme(legend.position="right") + labs(fill = "FS methods")+
  ylim(min_value,max_value)

multiplot(p3, p4,cols=2)
#table(res_internal1[,c("redim","dataset")])
#dim(res_internal1)
#mean(res_internal1$val_acc)
#table(res_internal1[,"model"])
#table(res_external1[,"representation"])


#ggplot(data=res_internal1, aes(x=dataset, y=time, fill=model)) + geom_boxplot() +
# stat_summary(fun.y=mean, colour="black", geom="point", 
#             shape=4, size=3)# +
#multiplot(p1,p2,p3, p4,cols=4)
#		time-fc:wt2d####
#set condition to filter the results:
#set condition:
res_external1=list()
res_internal1=list()
v_dataset = "wt2"
v_bin = "eqf"
v_bin_scale = "eqf_NA"
v_color_img = "gray"
#v_model = "cnn_l1f32"
v_model = "fc"
v_taxa = "freq.bug_species"
v_estop = "estop5"
v_padding = "no"
v_autofit_data = "scale_range_min_max"
v_redim = "pc"

#name for the chart
#title_chart <- paste("various representations: ",v_dataset,
#                     v_bin,
#                    v_color_img ,
#                     v_model,
#                     v_taxa ,
#                     v_estop )
#colnames(all_results_combined)
#filter results in internal validation
res_internal1 = all_results_combined[ all_results_combined$type_evaluate == "internal" & all_results_combined$whole == "no" 
                                      
                                      #all_results_combined$dataset== v_dataset
                                      #   all_results_combined$color_img== v_color_img 
                                      
                                      #  &  all_results_combined$model== v_model
                                      #  & all_results_combined$redim == v_redim
                                      # &  all_results_combined$bin_scale==v_bin_scale 
                                      #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                      # all_results_combined$padding == v_padding
                                      #  & all_results_combined$autofit_data== v_autofit_data 
                                      
                                      ,]
table(res_internal1[,c("redim","dataset")])
unique(res_internal1$filename)
unique(res_internal1$representation)
res_internal_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "no" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "internal"
                                        ,]
res_internal1 = rbind(res_internal1,res_internal_raw)
#filter results in external validation
res_external1 = all_results_combined[all_results_combined$type_evaluate == "external" & all_results_combined$whole == "yes" 
                                     
                                     #all_results_combined$dataset== v_dataset
                                     #   all_results_combined$color_img== v_color_img 
                                     
                                     &  all_results_combined$model== v_model
                                     
                                     # &  all_results_combined$bin_scale==v_bin_scale 
                                     #  all_results_combined$taxa == v_taxa & all_results_combined$estop == v_estop & 
                                     # all_results_combined$padding == v_padding & all_results_combined$autofit_data== v_autofit_data &
                                     
                                     ,]
unique(res_external1$filename)
unique(res_internal1$filename)
res_external_raw = all_results_combined[all_results_combined$dataset== v_dataset 
                                        &  all_results_combined$representation == "raw-1D" 
                                        &  all_results_combined$whole == "yes" 
                                        &  all_results_combined$taxa == v_taxa 
                                        & all_results_combined$type_evaluate == "external"
                                        ,]
res_external1 = rbind(res_external1,res_external_raw)
v_groups = "dataset"
v_comparison = "model"
mean_values <- aggregate(x = res_internal1, 
                         by = list(v_groups = res_internal1[,v_groups], v_comparison = res_internal1[,v_comparison]), 
                         FUN = "mean")

ratio_majority_class = 0.5414365
#training
#position y of significant results to show in the charts
y_signifi_position = 1.01
#get min, max for setting the charts
min_value = min(min(res_internal1$train_acc),min(res_internal1$val_acc))
max_value = max(max(res_internal1$train_acc),max(res_internal1$val_acc),y_signifi_position)
#the chart of training performance
res_internal1 <- res_internal1[!is.na(res_internal1$dataset),]
names(mean_values)[1] <- 'dataset'
p3 = 
  ggplot(data=res_internal1, aes(x=dataset, y=train_acc, fill=time)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) + #+geom_point() +
  #  stat_summary(fun.y=mean, geom="line", size=1.5,
  #              linetype="dotted")
  #ggtitle("Training on CNN1D model") + 
  #ggtitle(paste("Training on ",title_chart))  +
  #geom_text(data =  mean_values[,c("dataset","train_acc")], aes(label = train_acc, y = train_acc - 0.05)) +
  
  # geom_text(data = means_train, aes(label = train_acc, y = train_acc - 0.05)) +
  #sinificant
  #geom_point(data=signifi_comp_train,shape =8, size =3,
  #          aes(x=representation, y=y_signifi_position)) +
  theme(legend.position="none") +
  #geom_hline(aes(yintercept=ratio_majority_class),linetype = "dashed",show.legend=TRUE)+
  # geom_text(aes(2,ratio_majority_class,label = "majority_class", vjust = -1), size = 3)+
  ylim(min_value,max_value)

#testing and external, significant results


names(mean_values)[1] <- 'dataset'
#mean_values_val <- mean_values[,c("dataset","val_acc")]
p4 = 
  ggplot(data=res_internal1, aes(x=dataset, y=val_acc, fill=time))+
  #ggtitle("Training on CNN1D model") +
  # ggtitle(paste("Testing on ",title_chart))+ 
  geom_boxplot() +
  # stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3)   +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) +
  theme(legend.position="right") + labs(fill = "FS methods")+
  ylim(min_value,max_value)


multiplot(p3, p4,cols=2)
#table(res_internal1[,c("redim","dataset")])
#dim(res_internal1)
#mean(res_internal1$val_acc)
#table(res_internal1[,"model"])
#table(res_external1[,"representation"])

tmp_levels <- levels(res_internal1$redim)
tmp <- tmp_levels[3]
tmp_levels[3] <- tmp_levels[4]
tmp_levels[4] <- tmp
#levels(res_internal1$redim)
res_internal1$redim <- factor(res_internal1$redim, levels = tmp_levels)
ggplot(data=res_internal1, aes(x=dataset, y=time, fill=redim)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=3) + labs(fill = "Methods")+
    theme(axis.text=element_text(size=16),
            axis.title=element_text(size=16,face="bold"),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
table (res_internal1[,c("dataset","redim")])
#multiplot(p1,p2,p3, p4,cols=4)
#   tinh:p-value####
#fc_pp <- res_internal1
#fc_raw <- res_internal1
#fc_pp = res_internal1

da_v = "cir" # 2.2e-16
t.test(fc_raw[fc_raw$dataset==da_v,"val_acc"],fc_pp[fc_pp$dataset==da_v,"val_acc"])

da_v = "col" # 2.2e-16
t.test(fc_raw[fc_raw$dataset==da_v,"val_acc"],fc_pp[fc_pp$dataset==da_v,"val_acc"])

da_v = "obe" # 0.061
t.test(fc_raw[fc_raw$dataset==da_v,"val_acc"],fc_pp[fc_pp$dataset==da_v,"val_acc"])

da_v = "ibd" # 8.585e-15
t.test(fc_raw[fc_raw$dataset==da_v,"val_acc"],fc_pp[fc_pp$dataset==da_v,"val_acc"])

da_v = "t2d" # 7.934e-05
t.test(fc_raw[fc_raw$dataset==da_v,"val_acc"],fc_pp[fc_pp$dataset==da_v,"val_acc"])

da_v = "wt2" # 0.00615
t.test(fc_raw[fc_raw$dataset==da_v,"val_acc"],fc_pp[fc_pp$dataset==da_v,"val_acc"])
