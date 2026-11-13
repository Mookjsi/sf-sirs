# Direct simulation from the selected Phase 26 M2 evaluation. The evaluation RDS is
# supplied by an authorized reviewer; this module does not distribute fitted data.
frozen_m2_load <- function(path, expected_sha256=NULL) {
  if (!file.exists(path)) stop("Missing selected M2 evaluation: ", path)
  if (!is.null(expected_sha256)) {
    actual <- unname(tools::sha256sum(path))
    if (!identical(actual, expected_sha256)) stop("M2 evaluation SHA-256 mismatch")
  }
  x <- readRDS(path)
  if (!identical(x$metadata$phase, 26L) ||
      !identical(x$metadata$stage, "m2-evaluate") ||
      !identical(x$result$status, "EVALUATED") ||
      !identical(x$result$model, "M2")) stop("Wrong frozen M2 evaluation")
  score <- x$result$score
  w <- score$weekly
  required <- c("year", "week", "seasonal_week", "seasonal_beta",
                "effective_beta", "S", "I", "R", "incidence", "expected")
  if (!all(required %in% names(w)) || nrow(w) != 573L ||
      length(score$mask) != 573L || sum(score$mask) != 572L)
    stop("Incomplete frozen M2 weekly trajectory")
  list(score=score, weekly=w, source=path,
       source_sha256=unname(tools::sha256sum(path)))
}

frozen_m2_fixed <- function(model) {
  p <- model$score$parameters
  list(pop=sum(unlist(model$score$states[1L, c("S", "I", "R")],
                      use.names=FALSE)), mu=p$mu, gamma=p$gamma,
       omega=p$omega, delta=p$delta, rho=p$rho)
}

frozen_m2_check <- function(model, tolerance=1e-6) {
  w <- model$weekly
  fixed <- frozen_m2_fixed(model)
  beta <- vapply(1:52, function(i) {
    found <- which(w$seasonal_week == i)
    if (!length(found)) stop("Missing seasonal week ", i)
    w$seasonal_beta[found[1L]]
  }, numeric(1))
  first <- unlist(w[1L, c("S", "I", "R")], use.names=TRUE)
  sim <- simulate_weekly_sirs(beta, first, fixed, w$seasonal_week,
                              rep(1, nrow(w)), initial_mass=sum(first),
                              beta_by_transition=w$effective_beta)
  fields <- c(S="S", I="I", R="R", incidence="incidence",
              expected_cases="expected")
  errors <- vapply(names(fields), function(k)
    max(abs(sim[[k]] - w[[fields[[k]]]])), numeric(1))
  if (any(!is.finite(errors)) || any(errors > tolerance))
    stop("Frozen M2 replay differs: ", paste(names(errors), signif(errors, 4), collapse=", "))
  errors
}

frozen_m2_blend <- function(model, historical_rows=364L,
                            boundary_blend_weeks=2L) {
  w <- model$weekly
  if (nrow(w) != 573L || historical_rows != 364L ||
      boundary_blend_weeks != 2L) stop("Unexpected blend calendar")
  correction <- w$log_correction
  for (i in which(w$week == 27L & seq_len(nrow(w)) <= historical_rows)) {
    if (i > 1L) for (j in seq_len(boundary_blend_weeks)) {
      at <- i+j-1L
      if (at <= historical_rows)
        correction[at] <- (j/(boundary_blend_weeks+1L))*correction[at] +
          (1-j/(boundary_blend_weeks+1L))*w$log_correction[i-1L]
    }
  }
  w$effective_beta <- w$seasonal_beta*exp(correction)*w$future_q
  beta <- vapply(1:52, function(i)
    w$seasonal_beta[which(w$seasonal_week == i)[1L]], numeric(1))
  first <- unlist(w[1L, c("S", "I", "R")], use.names=TRUE)
  sim <- simulate_weekly_sirs(beta, first, frozen_m2_fixed(model),
                              w$seasonal_week, rep(1, nrow(w)),
                              initial_mass=sum(first),
                              beta_by_transition=w$effective_beta)
  for (field in c("S", "I", "R", "incidence"))
    w[[field]] <- sim[[field]]
  w$expected <- sim$expected_cases
  model$weekly <- w
  model
}

frozen_m2_clock <- function(model) {
  w <- model$weekly
  data.frame(row=seq_len(nrow(w)), year=w$year, week=w$week,
             seasonal_week=w$seasonal_week)
}

frozen_m2_profile <- function(model, year, u=0, z=0, v=0) {
  w <- model$weekly
  clock <- frozen_m2_clock(model)
  ix <- p27_window(clock, year, season=TRUE)
  at <- ix[1L] - 1L
  if (at < 1L) stop("No saved state before season")
  saved <- unlist(w[at, c("S", "I", "R")], use.names=TRUE)
  state <- p27_restart_state(saved, u=u, z=z)
  effective <- w$effective_beta[at:nrow(w)]
  effective[seq.int(2L, length(ix)+1L)] <-
    effective[seq.int(2L, length(ix)+1L)] * exp(v)
  beta <- vapply(1:52, function(i)
    w$seasonal_beta[which(w$seasonal_week == i)[1L]], numeric(1))
  sim <- simulate_weekly_sirs(beta, state, frozen_m2_fixed(model),
                              w$seasonal_week[at:nrow(w)],
                              rep(1, length(effective)),
                              initial_mass=sum(saved),
                              beta_by_transition=effective)
  list(expected=c(w$expected[seq_len(at)], sim$expected_cases[-1L]),
       states=sim, restart=at, initial_state=state, window=ix,
       effective_beta=effective)
}

frozen_m2_peak <- function(model, year, expected, id) {
  p27_peak_for_sensitivity(expected, frozen_m2_clock(model), year, id,
                           season=TRUE)
}

frozen_m2_eligible <- function(baseline, variant) {
  if (baseline$status != "ELIGIBLE" || variant$status != "ELIGIBLE")
    return(FALSE)
  f <- p10_assignment(baseline$peaks, variant$peaks,
                      p10_defaults()$identity)
  reverse <- p10_assignment(baseline$peaks, variant$peaks,
                            p10_defaults()$identity, TRUE)
  matched <- f[f$status == "MATCHED" & f$left_peak == baseline$peak, ]
  nrow(matched) == 1L && matched$right_peak == variant$peak &&
    p10_identity_consistent(f, reverse)
}

frozen_m2_h2 <- function(model, years=2013:2022) {
  steps <- log(1.2)/2^(2:6)
  rows <- list()
  for (year in years) {
    base <- frozen_m2_peak(model, year, model$weekly$expected,
                           paste(year, "base"))
    for (h in steps) {
      variants <- list()
      for (axis in c("S", "beta")) for (sgn in c(-1, 1)) {
        key <- paste(axis, sgn, sep="/")
        variants[[key]] <- tryCatch({
          run <- frozen_m2_profile(model, year,
                                   u=if (axis == "S") sgn*h else 0,
                                   v=if (axis == "beta") sgn*h else 0)
          peak <- frozen_m2_peak(model, year, run$expected, key)
          if (!frozen_m2_eligible(base, peak)) peak$status <- "BRANCH_MISMATCH"
          peak
        }, error=function(e) list(status=paste0("INFEASIBLE: ", conditionMessage(e)),
                                  time=NA_real_))
      }
      statuses <- vapply(variants, `[[`, character(1), "status")
      ok <- base$status == "ELIGIBLE" && all(statuses == "ELIGIBLE")
      slope <- function(axis) if (ok)
        (variants[[paste0(axis, "/1")]]$time -
         variants[[paste0(axis, "/-1")]]$time)/(2*h) else NA_real_
      rows[[length(rows)+1L]] <- data.frame(year=year, h=h,
        eligible=ok, reason=if (ok) "" else paste(unique(c(base$status,
          statuses[statuses != "ELIGIBLE"])), collapse=";"),
        E_S=slope("S"), E_beta=slope("beta"))
    }
  }
  steps_table <- do.call(rbind, rows)
  seasons <- do.call(rbind, lapply(years, function(year) {
    q <- steps_table[steps_table$year == year & steps_table$eligible, ]
    if (nrow(q)) {
      q <- q[order(q$h), ]
      chosen <- q[(nrow(q)+1L)%/%2L, ]
      data.frame(year=year, eligible=TRUE, selected_h=chosen$h,
                 E_S=chosen$E_S, E_beta=chosen$E_beta)
    } else data.frame(year=year, eligible=FALSE, selected_h=NA_real_,
                      E_S=NA_real_, E_beta=NA_real_)
  }))
  list(steps=steps_table, seasons=seasons)
}

frozen_m2_h3 <- function(model, years=2013:2022) {
  steps <- log(1.2)/2^(2:6)
  rows <- list(); contrasts <- list()
  for (year in years) {
    base <- frozen_m2_peak(model, year, model$weekly$expected,
                           paste(year, "base"))
    for (h in c(steps, log(1.2))) {
      sides <- lapply(c(-1, 1), function(sgn) tryCatch({
        run <- frozen_m2_profile(model, year, z=sgn*h)
        peak <- frozen_m2_peak(model, year, run$expected,
                               paste(year, h, sgn))
        if (!frozen_m2_eligible(base, peak)) peak$status <- "BRANCH_MISMATCH"
        peak
      }, error=function(e) list(status=paste0("INFEASIBLE: ", conditionMessage(e)),
                                time=NA_real_)))
      ok <- base$status == "ELIGIBLE" &&
        all(vapply(sides, `[[`, character(1), "status") == "ELIGIBLE")
      value <- if (ok) sides[[2L]]$time-sides[[1L]]$time else NA_real_
      row <- data.frame(year=year, h=h, eligible=ok,
                        minus_status=sides[[1L]]$status,
                        plus_status=sides[[2L]]$status,
                        contrast=value)
      if (h %in% steps) rows[[length(rows)+1L]] <- row
      else contrasts[[length(contrasts)+1L]] <- row
    }
  }
  steps_table <- do.call(rbind, rows)
  seasons <- do.call(rbind, lapply(years, function(year) {
    q <- steps_table[steps_table$year == year & steps_table$eligible, ]
    if (nrow(q)) {
      q <- q[order(q$h), ]
      chosen <- q[(nrow(q)+1L)%/%2L, ]
      data.frame(year=year, eligible=TRUE, selected_h=chosen$h,
                 E_I=chosen$contrast/(2*chosen$h))
    } else data.frame(year=year, eligible=FALSE, selected_h=NA_real_,
                      E_I=NA_real_)
  }))
  list(steps=steps_table, seasons=seasons,
       finite=do.call(rbind, contrasts))
}

frozen_m2_h1 <- function(model, declared_grid, parents) {
  required_grid <- c("scenario", "target", "policy", "path", "v")
  required_parent <- c("policy", "season", "path", "parent",
                       "parent_left", "parent_right")
  if (!all(required_grid %in% names(declared_grid)) ||
      !all(required_parent %in% names(parents))) stop("Incomplete H1 design inputs")
  declared_grid <- declared_grid[declared_grid$path == "MATCHED_BETA_ONLY", ]
  parents <- parents[parents$path == "MATCHED_BETA_ONLY", ]
  if (anyDuplicated(declared_grid$scenario) || any(!declared_grid$policy %in% c("raw", "blend")) ||
      any(!as.integer(declared_grid$target) %in% 2013:2022) ||
      any(abs(declared_grid$v) > log(1.2)+1e-12)) stop("Invalid H1 design grid")
  all_edges <- list(); all_parents <- list()
  clock <- frozen_m2_clock(model)
  for (policy in c("raw", "blend")) {
    active <- if (policy == "raw") model else frozen_m2_blend(model)
    for (year in 2013:2022) {
      design <- declared_grid[declared_grid$policy == policy &
                              as.integer(declared_grid$target) == year, ]
      design <- design[order(design$v), ]
      parent <- parents[parents$policy == policy & parents$season == year, ]
      parent <- parent[order(parent$parent), ]
      if (nrow(parent) != 86L || nrow(design) < 87L ||
          anyDuplicated(design$v) ||
          abs(design$v[1L]+log(1.2)) > 1e-10 ||
          abs(tail(design$v, 1L)-log(1.2)) > 1e-10)
        stop("Incomplete H1 design for ", policy, " ", year)
      runs <- lapply(seq_len(nrow(design)), function(i) {
        d <- design[i, ]
        profile <- frozen_m2_profile(active, year, v=d$v)
        r <- p27_peak_record(profile$expected, clock, year, d$scenario,
                             padding=13L, season=TRUE)
        r$grid$v <- d$v; r$grid$path <- "MATCHED_BETA_ONLY"
        if (nrow(r$peaks)) {
          r$peaks$v <- d$v; r$peaks$path <- "MATCHED_BETA_ONLY"
        }
        r$curve$v <- d$v; r$curve$path <- "MATCHED_BETA_ONLY"
        r
      })
      names(runs) <- design$scenario
      tracked <- p25_track(runs, phase25_config())
      gaps <- p27_height_gaps(tracked$grid, tracked$peaks, phase25_config())
      edges <- p27_classify(tracked$grid, tracked$peaks, gaps)
      curves <- p27_bind(lapply(runs, `[[`, "curve"))
      edges <- m2_reclassify(edges, tracked$grid, tracked$peaks, curves)
      edges$policy <- policy
      all_edges[[paste(policy, year)]] <- edges
      boundaries <- c(parent$parent_left, tail(parent$parent_right, 1L))
      index <- findInterval(edges$v_left, boundaries, rightmost.closed=TRUE)
      rows <- lapply(seq_len(nrow(parent)), function(j) {
        q <- edges[index == j, ]
        changing <- q[q$timing_change, , drop=FALSE]
        known <- changing$classification %in%
          c("CONTINUOUS_DISPLACEMENT", "DOMINANCE_EXCHANGE")
        data.frame(policy=policy, season=year, path="MATCHED_BETA_ONLY",
          parent=parent$parent[j], parent_left=parent$parent_left[j],
          parent_right=parent$parent_right[j],
          timing_event=nrow(changing)>0L,
          displacement=any(changing$classification=="CONTINUOUS_DISPLACEMENT"),
          exchange=any(changing$classification=="DOMINANCE_EXCHANGE"),
          unresolved=any(!known),
          mixed=length(unique(changing$classification))>1L,
          absolute_movement=sum(changing$absolute_movement))
      })
      all_parents[[paste(policy, year)]] <- do.call(rbind, rows)
    }
  }
  list(edges=do.call(rbind, all_edges),
       parents=do.call(rbind, all_parents))
}
