# ==============================================================================
# Analysis: 原样本20份插补数据的CBPS-ATE与br.logit ATO权重
# Date: 2026-07-29
# Random seed: 不新增随机过程
# R: 4.4.2
#
# Method:
# Purpose: 在每份锁定插补数据中拟合互补CBPS-ATE和当前投稿主ATO权重。
# Source type: GitHub / local_knowledge_base / method_paper
# Source name or URL:
#   https://github.com/ngreifer/WeightIt
#   https://github.com/ngreifer/cobalt
#   D:/source/obsidian/method_wiki/causal_inference/g_methods_ipw_standardization.md
# Citation or commit/version:
#   WeightIt commit 595d2812c2f4d56ff7c14aaf6fcde76c59aa8c8f
#   cobalt commit d955e34778a1ba44a6e437535eaa952ab855e0ee
#   Imai K, Ratkovic M. JRSSB. 2014;76:243-263.
#   Li F, Morgan KL, Zaslavsky AM. JASA. 2018;113:390-400.
# Access date: 2026-07-29
# Adaptation notes:
#   CBPS原始权重按试验内边际策略概率显式稳定化；ATO使用br.logit；
#   本脚本保存权重和设计诊断，不汇总治疗效应。
# Verification command:
#   Rscript --vanilla 04_倾向评分与权重/03_fit_cbps_ate_weights.R
# ==============================================================================

# 读取完整命令行以定位当前脚本。
# 使用用户锁定的绝对路径定位单次Boot-MI权重函数文件。
engine_path <- normalizePath(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "05_Bootstrap", "04_bootstrap_single_repeat.R"),
  winslash = "/",
  mustWork = TRUE
)
# 加载统一配置、锁定核心和br.logit权重函数。
source(engine_path, local = globalenv(), encoding = "UTF-8")
# 将工作目录固定到运行包根目录，保证交互式与命令行执行一致。
setwd(package_root)
# 核对当前脚本直接使用的包。
require_packages(c("mice", "digest"))

# 初始化逐插补权重诊断列表。
diagnostic_rows <- list()
# 初始化权重对象文件清单列表。
manifest_rows <- list()

# 逐个平行试验拟合原样本权重。
for (trial_index in seq_len(nrow(trial_registry))) {
  # 读取当前试验注册信息。
  trial_row <- trial_registry[trial_index, , drop = FALSE]
  # 保存当前试验标识。
  trial <- trial_row$trial
  # 构造当前试验MICE对象路径。
  mice_path <- file.path(
    output_root,
    "02_MICE",
    paste0("mids_", trial, "_m20_maxit50.rds")
  )
  # MICE对象不存在时停止并提示先运行插补。
  if (!file.exists(mice_path)) {
    # 报告缺失对象路径。
    stop("请先运行02_run_multiple_imputation.R；缺少：", mice_path)
  }
  # 读取当前试验MICE保存对象。
  saved_mice <- readRDS(mice_path)
  # 核对试验身份和插补数量。
  if (
    !identical(saved_mice$trial, trial) ||
      saved_mice$imputation$m != mice_m
  ) {
    # 身份或数量不一致时停止。
    stop(trial, "MICE对象身份或m值异常。")
  }
  # 逐份插补数据拟合两种预设权重。
  for (imputation_index in seq_len(mice_m)) {
    # 从mids对象提取当前完成数据。
    completed <- mice::complete(
      saved_mice$imputation,
      imputation_index
    )
    # 核对主PS公式字段不存在残余缺失。
    if (anyNA(completed[, all.vars(ps_formula_main), drop = FALSE])) {
      # 报告试验和插补编号。
      stop(trial, "第", imputation_index, "份完成数据仍有PS字段缺失。")
    }
    # 拟合主CBPS-ATE稳定化权重。
    cbps_fit <- fit_bootstrap_weights(
      completed,
      ps_formula_main,
      "cbps_ate"
    )
    # 拟合br.logit ATO重叠权重敏感性模型。
    ato_fit <- fit_bootstrap_weights(
      completed,
      ps_formula_main,
      "ato"
    )
    # 构造当前插补的患者追溯键。
    key_map <- data.frame(
      patient_id = saved_mice$raw_data$patient_id,
      hospital_admission_id = saved_mice$raw_data$hospital_admission_id,
      icu_stay_id = saved_mice$raw_data$icu_stay_id,
      stringsAsFactors = FALSE
    )
    # 将完成数据、两套权重和公式组织为单个可审计对象。
    weight_object <- list(
      trial = trial,
      imputation = imputation_index,
      key_map = key_map,
      completed_data = completed,
      ps_formula = ps_formula_main,
      cbps_ate = cbps_fit,
      ato_brlogit = ato_fit
    )
    # 构造当前插补权重对象输出路径。
    object_path <- file.path(
      output_root,
      "03_权重",
      paste0(
        "weights_",
        trial,
        "_imp",
        sprintf("%02d", imputation_index),
        ".rds"
      )
    )
    # 原子保存当前插补权重对象。
    atomic_save_rds(weight_object, object_path)
    # 对主CBPS和ATO分别登记诊断结果。
    for (model_name in c("cbps_ate", "ato_brlogit")) {
      # 选择当前模型的诊断对象。
      current_fit <- if (model_name == "cbps_ate") {
# 历史锁定输出使用CBPS-ATE对象；当前投稿主分析使用ATO，勿以本脚本旧输出覆盖ATO成品。
        cbps_fit
      } else {
        # 敏感性分析使用br.logit ATO对象。
        ato_fit
      }
      # 追加一行结构化诊断记录。
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
        trial = trial,
        imputation = imputation_index,
        model = model_name,
        max_absolute_smd = current_fit$max_absolute_smd,
        max_absolute_smd_term = current_fit$max_absolute_smd_term,
        ess_continuation = current_fit$ess_continuation,
        ess_deescalation = current_fit$ess_deescalation,
        max_weight = current_fit$max_weight,
        warning_n = current_fit$warning_n,
        stringsAsFactors = FALSE
      )
    }
    # 计算并登记当前权重对象的SHA-256。
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      trial = trial,
      imputation = imputation_index,
      relative_path = file.path(
        "09_输出结果",
        "03_权重",
        basename(object_path)
      ),
      sha256 = toupper(digest::digest(
        file = object_path,
        algo = "sha256",
        serialize = FALSE
      )),
      n = nrow(completed),
      stringsAsFactors = FALSE
    )
  }
}

# 合并全部40份插补的两模型诊断。
diagnostic_table <- do.call(rbind, diagnostic_rows)
# 合并全部40个权重对象的文件清单。
manifest_table <- do.call(rbind, manifest_rows)
# 原子写出权重诊断表。
atomic_write_csv(
  diagnostic_table,
  file.path(output_root, "03_权重", "weight_diagnostics.csv")
)
# 原子写出权重对象身份清单。
atomic_write_csv(
  manifest_table,
  file.path(output_root, "03_权重", "weight_object_manifest.csv")
)
# 向控制台报告拟合完成数量。
message("权重拟合完成：40份插补数据×2种模型，警告均由硬门槛处理。")
