args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6L || !identical(args[c(1L,3L,5L)],
    c("--h1", "--mechanisms", "--out")))
  stop("Usage: Rscript scripts/render_figure3.R --h1 H1_DIR --mechanisms MECHANISMS_DIR --out NEW_DIR")
paths <- vapply(args[c(2L,4L)], normalizePath, character(1), mustWork=TRUE)
out <- args[6L]
if (file.exists(out)) stop("Output path already exists: ", out)
dir.create(out, recursive=TRUE)
out <- normalizePath(out)
parents <- read.csv(file.path(paths[1L], "h1_parents.csv"))
h2 <- read.csv(file.path(paths[2L], "h2_seasons_raw.csv"))
h3_raw <- read.csv(file.path(paths[2L], "h3_seasons_raw.csv"))
h3_blend <- read.csv(file.path(paths[2L], "h3_seasons_blend.csv"))
years <- 2013:2022
fraction <- vapply(years, function(year) {
  q <- parents[parents$policy == "raw" & parents$season == year, ]
  n <- sum(q$timing_event)
  if (!n) return(NA_real_)
  sum(q$timing_event & q$displacement & !q$mixed & !q$unresolved)/n
}, numeric(1))
labels <- vapply(years, function(year) {
  q <- parents[parents$policy == "raw" & parents$season == year, ]
  paste0(sum(q$timing_event & q$displacement & !q$mixed & !q$unresolved),
         "/", sum(q$timing_event))
}, character(1))
grDevices::pdf(file.path(out, "figure_3_1.pdf"), width=11, height=10)
graphics::par(mfrow=c(3,1), mar=c(5,5,3,2))
x <- seq_along(years)
graphics::barplot(fraction, names.arg=paste0(years,"/",substr(years+1L,3L,4L)),
                  ylim=c(0,1.15), las=2, col="steelblue",
                  ylab="Pure displacement fraction", main="(a) H1 parent events")
graphics::abline(h=.5, lty=2)
graphics::text(x*1.2-.5, fraction+.06, labels=labels, cex=.7)
graphics::plot(x, h2$E_S, type="b", pch=16, col="steelblue",
               ylim=range(c(h2$E_S,h2$E_beta), na.rm=TRUE),
               xaxt="n", xlab="Season", ylab="Weeks per log change",
               main="(b) H2 local timing sensitivity")
graphics::axis(1, at=x, labels=years, las=2)
graphics::lines(x, h2$E_beta, type="b", pch=17, col="firebrick")
graphics::abline(h=0, lty=2)
graphics::legend("bottomright", legend=c("Initial S", "Transmission beta"),
                 col=c("steelblue", "firebrick"), lty=1, pch=c(16,17))
graphics::plot(x, h3_raw$E_I, type="b", pch=16, col="steelblue",
               ylim=range(c(h3_raw$E_I,h3_blend$E_I), na.rm=TRUE),
               xaxt="n", xlab="Season", ylab="Weeks per log change",
               main="(c) H3 initial I sensitivity")
graphics::axis(1, at=x, labels=years, las=2)
graphics::lines(x, h3_blend$E_I, type="b", pch=17, col="firebrick")
graphics::abline(h=0, lty=2)
graphics::legend("bottomright", legend=c("Primary raw", "Boundary blend"),
                 col=c("steelblue", "firebrick"), lty=1, pch=c(16,17))
grDevices::dev.off()
write.csv(data.frame(year=years, pure_displacement_fraction=fraction,
                     parent_event_label=labels, E_S=h2$E_S, E_beta=h2$E_beta,
                     E_I_raw=h3_raw$E_I, E_I_blend=h3_blend$E_I),
          file.path(out, "figure_3_1_data.csv"), row.names=FALSE)
