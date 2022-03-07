library(survminer)
library(survival)
library(preprocessCore)

OS.analysis <- function(sub.clinical.group, colors, sig.reg = F, cutoff = T, prefix = 'OS', title = NULL, out.figs.dir = './', risk.table = T) {
    if (cutoff) {
        cutoff <- surv_cutpoint(
            sub.clinical.group,
            time = "Time",
            event = "Sur",
            variables = c("Groups")
        )$cutpoint[[1]]
        sub.clinical.group$Groups <- ifelse(sub.clinical.group$Groups <= cutoff, 'Low risk', 'High risk')
    }
    if (sig.reg) {
        groups.diff <- survdiff(Surv(Time, Sur) ~ Groups, data = as.data.frame(sub.clinical.group))
        p.val <- pchisq(groups.diff$chisq, length(groups.diff$n)-1, lower.tail = FALSE)
        return(list(p.val = p.val, cut.off = cutoff, group = sub.clinical.group$Groups))
    }

    fit <- survfit(Surv(Time, Sur) ~ Groups, data = as.data.frame(sub.clinical.group))
    out.file <- file.path(out.figs.dir, prefix)
	pdf(out.file, width = 8, height = 8)
    ggsurv <- ggsurvplot(
        fit,
        data = sub.clinical.group,
        pval = T,
        risk.table.col="strata",
        risk.table=risk.table,
        palette = colors,
        title = title,
        ggtheme = theme_classic2(base_size=18) + theme(plot.title=element_text(hjust=0.5)))
	print(ggsurv)
    dev.off()
	print(ggsurv)
    return(list(gplot = ggsurv, cutoff = cutoff, group = sub.clinical.group$Groups))
}

quantileNorm <- function(x) {
    Y <- as.matrix(x)
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
    return(Y)
}

pairwiseComp <- function(lab.vec) {
    combn(lab.vec, 2) %>% { lapply(1 : dim(.)[2], function(x) {.[, x]}) }
}




