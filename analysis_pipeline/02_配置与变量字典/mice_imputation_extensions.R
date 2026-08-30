# ==============================================================================
# Analysis: MICE二分类插补数值收敛扩展
# Date: 2026-07-31
#
# Method:
# Purpose:
#   保留mice::mice.impute.logreg的统计方法与数据增广步骤；针对两个锁定目标，
#   内部stats::glm.fit仅在精确的不收敛警告下提高迭代上限并重试一次。
# Source type: GitHub / package_source / method_paper
# Source name or URL:
#   https://github.com/amices/mice/blob/61083e667fa41cc4bc655492591b100bf8a7e443/R/mice.impute.logreg.R
# Citation or commit/version:
#   amices/mice commit 61083e667fa41cc4bc655492591b100bf8a7e443;
#   White IR, Daniel R, Royston P. Stat Med. 2010;29:2920-2931.
# Access date: 2026-07-31
# Adaptation notes:
#   原函数的数据增广、quasibinomial-logit拟合、系数后验随机抽取和
#   Bernoulli随机插补均保持不变。nephrotoxic_drug_exposure先使用100次；
#   source_control_status保留官方默认25次。仅在捕获精确的
#   “glm.fit: algorithm did not converge”且fit$converged=FALSE时，
#   以相同数据、公式和起始规则使用maxit=2000重试一次，并记录审计信息。
#   两个扩展仅绑定抗PSA的对应锁定目标；不改变m=20、MICE链maxit=50、
#   预测矩阵、访问顺序、donors、随机种子或其他变量的插补方法。
# Verification command:
#   Rscript --vanilla 05_Bootstrap/05_run_bootstrap.R --trial=anti_psa --b_start=16 --b_end=16
# ==============================================================================

# 定义单次Bootstrap重复内的logreg条件式重试审计环境。
.logreg100_retry_audit_environment <- new.env(parent = emptyenv())
.logreg100_retry_audit_environment$call_n <- 0L
.logreg100_retry_audit_environment$rows <- list()

# 重置当前Bootstrap重复的logreg调用和重试记录。
reset_logreg100_retry_audit <- function() {
  .logreg100_retry_audit_environment$call_n <- 0L
  .logreg100_retry_audit_environment$rows <- list()
  invisible(NULL)
}

# 读取当前Bootstrap重复的logreg条件式重试记录。
get_logreg100_retry_audit <- function() {
  if (length(.logreg100_retry_audit_environment$rows) == 0L) {
    return(data.frame(
      target_field = character(0),
      locked_method = character(0),
      call_index = integer(0),
      initial_maxit = integer(0),
      initial_iterations = integer(0),
      retry_maxit = integer(0),
      retry_iterations = integer(0),
      observed_n = integer(0),
      missing_n = integer(0),
      observed_event_n = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, .logreg100_retry_audit_environment$rows)
}

mice.impute.logreg100 <- function(y, ry, x, wy = NULL, ...) {
  # 记录本次MICE链中的logreg内部调用顺序。
  .logreg100_retry_audit_environment$call_n <-
    .logreg100_retry_audit_environment$call_n + 1L
  current_call_index <- .logreg100_retry_audit_environment$call_n

  # 未显式指定待插补位置时，使用非观测位置。
  if (is.null(wy)) {
    wy <- !ry
  }

  # 复用mice官方数据增广，降低完全预测造成无限估计的风险。
  augmented <- mice:::augment(y, ry, x, wy)
  x <- augmented$x
  y <- augmented$y
  ry <- augmented$ry
  wy <- augmented$wy
  weight <- augmented$w

  # 为二项logistic模型加入截距列。
  x <- cbind(1, as.matrix(x))
  # 首次拟合固定使用已审计的100次内部迭代上限。
  initial_warnings <- character(0)
  fit <- withCallingHandlers(
    stats::glm.fit(
      x = x[ry, , drop = FALSE],
      y = y[ry],
      family = stats::quasibinomial(link = "logit"),
      weights = weight[ry],
      control = stats::glm.control(maxit = 100)
    ),
    warning = function(condition) {
      initial_warnings <<- c(
        initial_warnings,
        conditionMessage(condition)
      )
      invokeRestart("muffleWarning")
    }
  )

  # 只有精确的不收敛警告允许以2000次上限重试一次。
  exact_nonconvergence <- (
    length(initial_warnings) == 1L &&
      identical(
        initial_warnings,
        "glm.fit: algorithm did not converge"
      ) &&
      !isTRUE(fit$converged)
  )
  if (exact_nonconvergence) {
    initial_iterations <- as.integer(fit$iter)
    retry_warnings <- character(0)
    retry_fit <- withCallingHandlers(
      stats::glm.fit(
        x = x[ry, , drop = FALSE],
        y = y[ry],
        family = stats::quasibinomial(link = "logit"),
        weights = weight[ry],
        control = stats::glm.control(maxit = 2000)
      ),
      warning = function(condition) {
        retry_warnings <<- c(
          retry_warnings,
          conditionMessage(condition)
        )
        invokeRestart("muffleWarning")
      }
    )
    if (
      length(retry_warnings) > 0L ||
        !isTRUE(retry_fit$converged)
    ) {
      stop(
        "logreg100内部glm.fit条件式重试失败；警告=",
        paste(unique(retry_warnings), collapse = " | "),
        "；converged=",
        isTRUE(retry_fit$converged)
      )
    }

    # 保存方法不变、仅提高优化迭代上限的重试审计信息。
    .logreg100_retry_audit_environment$rows[[
      length(.logreg100_retry_audit_environment$rows) + 1L
    ]] <- data.frame(
      target_field = "nephrotoxic_drug_exposure",
      locked_method = "logreg",
      call_index = current_call_index,
      initial_maxit = 100L,
      initial_iterations = initial_iterations,
      retry_maxit = 2000L,
      retry_iterations = as.integer(retry_fit$iter),
      observed_n = sum(ry),
      missing_n = sum(wy),
      observed_event_n = sum(as.numeric(y[ry]) == 2L),
      stringsAsFactors = FALSE
    )
    fit <- retry_fit
  } else if (
    length(initial_warnings) > 0L ||
      !isTRUE(fit$converged)
  ) {
    # 其他警告或无精确警告的不收敛仍然硬停止。
    stop(
      "logreg100内部glm.fit出现未批准异常；警告=",
      paste(unique(initial_warnings), collapse = " | "),
      "；converged=",
      isTRUE(fit$converged)
    )
  }

  fit_summary <- summary.glm(fit)
  beta <- stats::coef(fit)
  random_variation <- t(chol(mice:::sym(fit_summary$cov.unscaled)))
  beta_star <- beta + random_variation %*% stats::rnorm(ncol(random_variation))
  probability <- 1 / (
    1 + exp(-(x[wy, , drop = FALSE] %*% beta_star))
  )
  imputed <- stats::runif(nrow(probability)) <= probability
  imputed[imputed] <- 1

  # 恢复原二分类因子水平。
  if (is.factor(y)) {
    imputed <- factor(imputed, c(0, 1), levels(y))
  }

  imputed
}

# 定义仅用于anti-PSA source_control_status的标准logreg条件式重试实现。
mice.impute.sourcecontrollogreg <- function(y, ry, x, wy = NULL, ...) {
  # 记录本次MICE链中的目标内部调用顺序。
  .logreg100_retry_audit_environment$call_n <-
    .logreg100_retry_audit_environment$call_n + 1L
  current_call_index <- .logreg100_retry_audit_environment$call_n

  # 未显式指定待插补位置时，使用非观测位置。
  if (is.null(wy)) {
    wy <- !ry
  }

  # 完全复用mice官方logreg的数据增广步骤。
  augmented <- mice:::augment(y, ry, x, wy)
  x <- augmented$x
  y <- augmented$y
  ry <- augmented$ry
  wy <- augmented$wy
  weight <- augmented$w

  # 为二项logistic模型加入截距列。
  x <- cbind(1, as.matrix(x))
  # 首次拟合显式保留mice官方glm.fit默认的25次迭代上限。
  initial_warnings <- character(0)
  fit <- withCallingHandlers(
    stats::glm.fit(
      x = x[ry, , drop = FALSE],
      y = y[ry],
      family = stats::quasibinomial(link = "logit"),
      weights = weight[ry],
      control = stats::glm.control(maxit = 25)
    ),
    warning = function(condition) {
      initial_warnings <<- c(
        initial_warnings,
        conditionMessage(condition)
      )
      invokeRestart("muffleWarning")
    }
  )

  # 仅对精确且唯一的不收敛警告提高上限至2000并重试一次。
  exact_nonconvergence <- (
    length(initial_warnings) == 1L &&
      identical(
        initial_warnings,
        "glm.fit: algorithm did not converge"
      ) &&
      !isTRUE(fit$converged)
  )
  if (exact_nonconvergence) {
    initial_iterations <- as.integer(fit$iter)
    retry_warnings <- character(0)
    retry_fit <- withCallingHandlers(
      stats::glm.fit(
        x = x[ry, , drop = FALSE],
        y = y[ry],
        family = stats::quasibinomial(link = "logit"),
        weights = weight[ry],
        control = stats::glm.control(maxit = 2000)
      ),
      warning = function(condition) {
        retry_warnings <<- c(
          retry_warnings,
          conditionMessage(condition)
        )
        invokeRestart("muffleWarning")
      }
    )
    if (
      length(retry_warnings) > 0L ||
        !isTRUE(retry_fit$converged)
    ) {
      stop(
        "source_control_status内部glm.fit条件式重试失败；警告=",
        paste(unique(retry_warnings), collapse = " | "),
        "；converged=",
        isTRUE(retry_fit$converged)
      )
    }

    # 保存方法不变、仅提高优化迭代上限的重试审计信息。
    .logreg100_retry_audit_environment$rows[[
      length(.logreg100_retry_audit_environment$rows) + 1L
    ]] <- data.frame(
      target_field = "source_control_status",
      locked_method = "logreg",
      call_index = current_call_index,
      initial_maxit = 25L,
      initial_iterations = initial_iterations,
      retry_maxit = 2000L,
      retry_iterations = as.integer(retry_fit$iter),
      observed_n = sum(ry),
      missing_n = sum(wy),
      observed_event_n = sum(as.numeric(y[ry]) == 2L),
      stringsAsFactors = FALSE
    )
    fit <- retry_fit
  } else if (
    length(initial_warnings) > 0L ||
      !isTRUE(fit$converged)
  ) {
    # 其他警告或无精确警告的不收敛仍然硬停止。
    stop(
      "source_control_status内部glm.fit出现未批准异常；警告=",
      paste(unique(initial_warnings), collapse = " | "),
      "；converged=",
      isTRUE(fit$converged)
    )
  }

  # 以下后验抽样和Bernoulli插补步骤与mice官方logreg一致。
  fit_summary <- summary.glm(fit)
  beta <- stats::coef(fit)
  random_variation <- t(chol(mice:::sym(fit_summary$cov.unscaled)))
  beta_star <- beta + random_variation %*% stats::rnorm(ncol(random_variation))
  probability <- 1 / (
    1 + exp(-(x[wy, , drop = FALSE] %*% beta_star))
  )
  imputed <- stats::runif(nrow(probability)) <= probability
  imputed[imputed] <- 1

  # 恢复原二分类因子水平。
  if (is.factor(y)) {
    imputed <- factor(imputed, c(0, 1), levels(y))
  }

  imputed
}
