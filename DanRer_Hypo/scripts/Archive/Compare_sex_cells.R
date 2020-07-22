LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                            label = i, size = text.size)
    plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                              adj.y.s, yend = y1, size = segment.size)
  }
  return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
                    adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                    adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}



cluster.36 <- SubsetData(hypo, ident.use = 36, subset.raw = T)
cluster.36 <- SetAllIdent(cluster.36, id = "sex")
avg.cluster.36 <- log1p(AverageExpression(cluster.36, show.progress = FALSE))
avg.cluster.36$gene <- rownames(avg.cluster.36)

cluster.9 <- SubsetData(hypo, ident.use = 9, subset.raw = T)
cluster.9 <- SetAllIdent(cluster.9, id = "sex")
avg.cluster.9 <- log1p(AverageExpression(cluster.9, show.progress = FALSE))
avg.cluster.9$gene <- rownames(avg.cluster.9)

cluster.0 <- SubsetData(hypo, ident.use = 0, subset.raw = T)
cluster.0 <- SetAllIdent(cluster.0, id = "sex")
avg.cluster.0 <- log1p(AverageExpression(cluster.0, show.progress = FALSE))
avg.cluster.0$gene <- rownames(avg.cluster.0)

cluster.4 <- SubsetData(hypo, ident.use = 4, subset.raw = T)
cluster.4 <- SetAllIdent(cluster.4, id = "sex")
avg.cluster.4 <- log1p(AverageExpression(cluster.4, show.progress = FALSE))
avg.cluster.4$gene <- rownames(avg.cluster.4)

genes.to.label1 = c("rps29", "rpl29", "histh1l", "fkb5")
genes.to.label2 = c("pnrc2", "hsp70l", "mvp", "cbx7a")
genes.to.label3 = c("ncf1", "dnajb1b", "thy1", "mtbl")
p1 <- ggplot(avg.cluster.36, aes(male, female)) + geom_point() + ggtitle("Cluster.36")
p1 <- LabelUR(p1, genes = c(genes.to.label1, genes.to.label2), avg.cluster.36, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p1 <- LabelUL(p1, genes = genes.to.label3, avg.cluster.36, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)

p2 <- ggplot(avg.cluster.9, aes(male, female)) + geom_point() + ggtitle("Cluster.9")
p2 <- LabelUR(p2, genes = c(genes.to.label1, genes.to.label3), avg.cluster.9, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p2 <- LabelUL(p2, genes = genes.to.label2, avg.cluster.9, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)

p3 <- ggplot(avg.cluster.0, aes(male, female)) + geom_point() + ggtitle("Cluster.0")
p3 <- LabelUR(p3, genes = c(genes.to.label1, genes.to.label3), avg.cluster.0, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p3 <- LabelUL(p3, genes = genes.to.label2, avg.cluster.0, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)

p4 <- ggplot(avg.cluster.4, aes(male, female)) + geom_point() + ggtitle("Cluster.4")
p4 <- LabelUR(p4, genes = c(genes.to.label1, genes.to.label3), avg.cluster.4, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p4 <- LabelUL(p4, genes = genes.to.label2, avg.cluster.4, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)

pdf("Figures/test.pdf", height = 14, width = 14)
plot_grid(p1, p2, p3, p4)
dev.off()




cluster.36 <- SetAllIdent(cluster.36, "sex")
cluster.36.sex.markers <- FindAllMarkers(cluster.36, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

cluster.9 <- SetAllIdent(cluster.9, "sex")
cluster.9.sex.markers <- FindAllMarkers(cluster.9, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)



