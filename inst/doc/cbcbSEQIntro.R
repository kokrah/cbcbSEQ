### R code from vignette source 'cbcbSEQIntro.Rnw'

###################################################
### code chunk number 1: options (eval = FALSE)
###################################################
## options(width=70)


###################################################
### code chunk number 2: cbcbSEQIntro.Rnw:65-73
###################################################
require(pasilla)
# locate the path of the dataset and read in the dataset
datafile = system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
counts = read.table(datafile, header=TRUE, row.names=1)
head(counts)
dim(counts)
counts = counts[rowSums(counts) > ncol(counts),]
dim(counts)


###################################################
### code chunk number 3: cbcbSEQIntro.Rnw:78-84
###################################################
design = data.frame(row.names=colnames(counts),
                    condition=c("untreated","untreated","untreated",
                                "untreated","treated","treated","treated"),
                    libType=c("single-end","single-end","paired-end",
                              "paired-end","single-end","paired-end","paired-end"))
design


###################################################
### code chunk number 4: cbcbSEQIntro.Rnw:89-103
###################################################
# load batch package
require(cbcbSEQ)
#
# quantile normalize: adjust counts for library size.
qcounts = qNorm(counts)
# convert counts to log2 counts per milliom. (voom scale)
cpm = log2CPM(qcounts)
names(cpm)
libsize = cpm$lib.size
cpm = cpm$y
#
# PCA analysis
# returns a list with two components v and d.
res = makeSVD(cpm)


###################################################
### code chunk number 5: cbcbSEQIntro.Rnw:108-109
###################################################
pcRes(res$v,res$d, design$condition, design$libType)


###################################################
### code chunk number 6: cbcbSEQIntro.Rnw:112-118
###################################################
plotPC(res$v,res$d,
       col=design$condition, # color by batch
       pch=19, main="PCA plot",
       xlim=c(min(res$v[,1])-.08,max(res$v[,1])+.08),
	     ylim=c(min(res$v[,2])-.08,max(res$v[,2])+.08))
text(res$v[,1], res$v[,2], design$libType, pos=1, cex=0.6)


###################################################
### code chunk number 7: cbcbSEQIntro.Rnw:129-137
###################################################
# combatMod function
# noScale=TRUE option not to scale adjust
tmp = combatMod(cpm, batch=design$libType, mod=design$condition, noScale=TRUE)

# look at PCA results again
res = makeSVD(tmp)
# batch effect is reduced
pcRes(res$v,res$d, design$condition, design$libType)


###################################################
### code chunk number 8: cbcbSEQIntro.Rnw:139-145
###################################################
plotPC(res$v,res$d,
       col=design$condition, # color by batch
       pch=19, main="PCA plot",
       xlim=c(min(res$v[,1])-.08,max(res$v[,1])+.08),
       ylim=c(min(res$v[,2])-.08,max(res$v[,2])+.08))
text(res$v[,1], res$v[,2], design$libType, pos=1, cex=0.6)


###################################################
### code chunk number 9: cbcbSEQIntro.Rnw:149-151
###################################################
v = voomMod(tmp, model.matrix(~design$condition), lib.size=libsize)
v$plot


###################################################
### code chunk number 10: cbcbSEQIntro.Rnw:154-158
###################################################
summary(v)
fit = lmFit(v)
eb = eBayes(fit)
top = topTable(eb, coef=2, n=nrow(v$E))


###################################################
### code chunk number 11: cbcbSEQIntro.Rnw:161-165
###################################################
sel = top$adj.P.Val < 0.05
plot(top$logFC, -log10(top$adj.P.Val), pch=16, cex=0.3,
     main=paste(sum(sel), "/", length(sel)),col=ifelse(sel,"red","black"))
abline(v=c(-1,1), h=-log10(0.05), col="blue")


###################################################
### code chunk number 12: cbcbSEQIntro.Rnw:169-177
###################################################
cond=design$condition
batch=design$libType
mod = model.matrix(~cond+batch ,
                   contrasts.arg=list(cond="contr.treatment", batch="contr.sum"))
v1 = voom(counts, mod)
fit1 = lmFit(v1)
eb1 = eBayes(fit1)
top1 = topTable(eb1, coef=2, n=nrow(v1$E))


###################################################
### code chunk number 13: cbcbSEQIntro.Rnw:180-184
###################################################
top$ID = rownames(top)
top1$ID = rownames(top1)
tab = merge(top[,c("ID", "adj.P.Val")], top1[,c("ID", "adj.P.Val")], by="ID")
as.data.frame(table(combat = tab[,2] < 0.05, model = tab[,3] < 0.05))


###################################################
### code chunk number 14: sessionInfo
###################################################
toLatex(sessionInfo())


