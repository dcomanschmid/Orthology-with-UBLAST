###########################################################################
# Reciprocal Best Hits from "USEARCH/UBLAST" output                       #
# proteome wide orthologs                                                 #       
#	1. select only orthologs with alignment length of at least 50 AA  #
# 	2. sort (on E-values) the query <> target ublast output           #
# 		! BHQ: best hit/s of the query = target protein/s with    # 
#                 the smallest E-value                                    #
#               ! more than one BHQ = co-orthologs                        #
#       3. for each query protein make a list with BHQ target proteins    #
#       4. for each query protein retrieve the orthologs of the BHQ       #
#          target proteins                                                #
#       	4.1. sort (on E-values) the target <> query ublast output #
#               	! BHT: best hit/s of the target                   #
# 	5. RBH: reciprocal best hit/s: "BHQ = BHT"                        #
#                                                                         #
# Run R from terminal :)                                                  #
#                                                                         #
# Diana Coman Schmid                                                      #
# Eawag 2015                                                              #
# diana.comanschmid@eawag.ch                                              #
###########################################################################


rm(list=ls())

workDir <- "/media/dianacs/data1/ortho_usearch/"
setwd(workDir)

# read in the usearch>ublast output file for the zebrafish <> trout comparison

zt.ublast <- read.table(file.path(workDir,"ZT.u8"),header=FALSE)
dim(zt.ublast)

# add column names (standard usearch>ublast file output: http://www.drive5.com/usearch/manual/blast6out.html

colnames(zt.ublast) <- c("Query_id","Target_id","Perc_ident","Aln_len","No_mism","No_gap","Start_pos_query","End_pos_query","Start_pos_target","End_pos_target","E_value","Bit_score")
head(zt.ublast)

# select alignments with length of at least 50 AA

zt.s <- zt.ublast[which(zt.ublast$Aln_len >= 50),c(1,2,3,4,11)]


# for each query protein, extract the usearch>ublast output for their corresponding target/s 
# sort ascending based on E-values
# store matrices in a list

queries.z <- unique(zt.s$Query_id)
str(queries.z)

atBeginning <- Sys.time()
then <- Sys.time()

queries.z.l <- list()
for (q in 1:length(queries.z)){
  q.m <- zt.s[which(zt.s$Query_id == queries.z[q]),]
  queries.z.l[[q]] <- q.m[order(q.m$E_value),]
  names(queries.z.l)[q] <- paste(queries.z[q])
}

print(Sys.time()-then)

# Time difference of 2.024467 hours

# for each query protein, extract the BHQ i.e. target protein/s with the smallest E-value (store as list)

atBeginning <- Sys.time()
then <- Sys.time()

targets.z.l <- list()
for (tg in names(queries.z.l)){
  targets.z.l[[tg]] <- queries.z.l[[tg]][which(queries.z.l[[tg]]$E_value == min(queries.z.l[[tg]]$E_value)),2]
}
print(Sys.time()-then)

# Time difference of 46.98568 secs


# read in the usearch>ublast output file for the trout <> zebrafish comparison

tz.ublast <- read.table(file.path(workDir,"TZ.u8"),header=FALSE)
dim(tz.ublast)

# add column names (standard usearch>ublast file output: http://www.drive5.com/usearch/manual/blast6out.html

colnames(tz.ublast) <- c("Query_id","Target_id","Perc_ident","Aln_len","No_mism","No_gap","Start_pos_query","End_pos_query","Start_pos_target","End_pos_target","E_value","Bit_score")
head(tz.ublast)

# select alignments with length of at least 50 AA

tz.s <- tz.ublast[which(tz.ublast$Aln_len >= 50),c(1,2,3,4,11)]
tz.s <- tz.s[order(tz.s$E_value),]

# for each query protein, extract the usearch>ublast output for their corresponding BHQ
# sort ascending based on E-values
# store matrices in a list

atBeginning <- Sys.time()
then <- Sys.time()

queries.t.l <- list()
for (q in names(targets.z.l)){
  q.m <- tz.s[which(tz.s$Query_id %in% targets.z.l[[q]]),]
  queries.t.l[[q]] <- q.m[order(q.m$E_value),]
}

print(Sys.time()-then)

# Time difference of 34.72711 mins

# for the BHQ of each query protein, if the BHT (protein/s with the smallest E-value) match the query protein => assign orthology as RBH
# store matrices in a list

atBeginning <- Sys.time()
then <- Sys.time()

rbh.l <- list()
for (r in names(queries.t.l)){
    mini <- queries.t.l[[r]][which(queries.t.l[[r]]$E_value == min(queries.t.l[[r]]$E_value)),]
    rbh.l[[r]] <- mini[which(mini$Target_id %in% unique(queries.z.l[[paste(r)]]$Query_id)),]
}

print(Sys.time()-then)

# Time difference of 8.681766 mins

summary(rbh.l)

# store the the RBH as matrix 
atBeginning <- Sys.time()
then <- Sys.time()

rbh.m <- do.call("rbind",rbh.l)
rownames(rbh.m) <- NULL

print(Sys.time()-then)

# Time difference of 3.460173 mins

head(rbh.m)

# save the RBH matrix

write.table(rbh.m,file.path(workDir,"zebrafish_trout_RBHublast.txt"),quote=FALSE)


