myel_CGA01<- cbind(table(Idents(myel_s$CGA01_S)), table(Idents(myel_s$CGA01_T)))
head(myel_CGA01)
colnames(myel_CGA01)<- c("D", "P")
write.table(myel_CGA01, "myel_CGA01.txt")

write.table(prop.table(myel_CGA01,  margin = 2), "prop_CGA01_myel.txt")



a<-read.csv2("cell_CGA01.csv", row.names = 1)
colnames(a)<-c("P", "D")
a<- as.matrix(a)

chi1<- chisq.test(a)

chi1$p.value

write.csv2(chi1$expected, "chitest_CGA01.csv")

b<-read.csv2("cell_CGA02.csv", row.names = 1)
colnames(b)<-c("P", "D")
b<- as.matrix(b)

chi2<- chisq.test(b)

chi2$p.value

write.csv2(chi2$expected, "chitest_CGA02.csv")

c<-read.csv2("cell_CGA03.csv", row.names = 1)
colnames(c)<-c("P", "D")
c<- as.matrix(c)

chi3<- chisq.test(c)

chi3$p.value

write.csv2(chi3$expected, "chitest_CGA03.csv")
### willcox test
#creo oggetto con i tre soggetti divisi per P e D 
#data la variabilita delle conte cellulari uso le frequenze

a<-read.csv2("freq_CGA01.csv", row.names = 1)
b<-read.csv2("freq_CGA02.csv", row.names = 1)
c<-read.csv2("freq_CGA03.csv", row.names = 1)


W<- cbind(a$V2, b$V2, c$V2, a$V1, b$V1, c$V1)
colnames(W)<- c("01D", "02D", "03D","01P", "02P", "03P")
rownames(W)<- as.numeric(0:19)

 wilcox.test(W[21, c(1:3)], W[21, c(4:6)], paired = TRUE)


 wilcox.test(W[3, c(1:3)], W[3, c(4:6)], paired = TRUE)


#W.list <- sapply(X = W, FUN = function(x) {
 #  x<- wilcox.test(W[1, c(1:3)], W[1, c(4:6)], paired = TRUE)
#})
#for (n in c(rownames(W))) w<-(wilcox.test(W[n, c(1:3)], W[n, c(4:6)], paired = TRUE))
