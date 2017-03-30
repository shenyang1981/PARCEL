a = read.table("mapping_combined.out",header=F,sep="\t")
colnames(a) = c("sample","Total","After Trimming","Mapped","condition")

a[,"id"] = gsub("(.*)_(.*)","\\1",a[,1],perl=T)
a[1:10,]
total = aggregate(a[,"V4"],list(a[,"id"]),sum)
total[1:10,]
ratio = aggregate(a[,"V4"],list(a[,"id"]),function(x){sd(x)/mean(x)});
mycv = aggregate(a[,"V4"],list(a[,"id"]),function(x){sd(x)/mean(x)});
mycv[1:10,]
dim(mycv)
unique(mycv[,1])
sum(mycv[,2]>0.1)
sum(mycv[,2]>0.2)
result = merge(final,mycv,by="Group.1")
result = merge(total,mycv,by="Group.1")
result[1:10,]
result[order(result[,2],decreasing=T),][1:10,]
result[order(result[,2],decreasing=T),][1:30,]
result[order(result[,2],decreasing=T),]
lowcv = result[result[,3]<0.2,]
lowcv[1:10,]
lowcv[order(lowcv[,2],decreasing=T),][1:10,]
candidates = c("tpp","sam","atp",lowcv[order(lowcv[,2],decreasing=T),][1:10,2])
candidates
candidates = c("tpp","sam","atp",lowcv[order(lowcv[,2],decreasing=T),][1:10,1])

