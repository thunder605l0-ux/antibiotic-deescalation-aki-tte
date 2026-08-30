# ==============================================================================
# 独立运行包静态和输入完整性审计
# Method: 解析全部R脚本、核对数据哈希、检查方法来源块和绝对项目路径
# Purpose: 在正式运行前发现语法、文件身份或不可移植依赖问题
# Source type: self_written
# Source name or URL: 本项目临床框架可复现性要求
# Citation or commit/version: D-4.7N-ANALYSIS-LOCK-V2.0
# Access date: 2026-07-29
# Adaptation notes: 只读检查，不运行MICE、权重、Bootstrap或结局模型。
# Verification command: Rscript 10_方法来源与文件校验/validate_package.R
# ==============================================================================

# 读取统一配置。
source(
  file.path(Sys.getenv("TTE_RUNTIME_ROOT", unset = "."), "02_配置与变量字典", "analysis_config.R"),
  encoding = "UTF-8"
)
setwd(package_root)

# 核对SHA-256计算依赖。
require_packages("digest")

# 递归列出全部R脚本。
r_files <- list.files(
  package_root,
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

# 初始化脚本QA。
script_qa <- lapply(r_files, function(script_path) {
  # 读取脚本文本。
  script_text <- readLines(
    script_path,
    warn = FALSE,
    encoding = "UTF-8"
  )
  # 尝试解析脚本。
  parse_ok <- tryCatch(
    {
      parse(file = script_path)
      TRUE
    },
    error = function(error_condition) FALSE
  )
  # 检查主要方法来源字段是否齐全。
  source_fields <- c(
    "# Method:",
    "# Purpose:",
    "# Source type:",
    "# Source name or URL:",
    "# Citation or commit/version:",
    "# Access date:",
    "# Adaptation notes:",
    "# Verification command:"
  )
  # 返回单脚本QA。
  data.frame(
    file = substring(
      normalizePath(script_path, winslash = "/"),
      nchar(normalizePath(package_root, winslash = "/")) + 2L
    ),
    parse_ok = parse_ok,
    source_block_complete = all(
      vapply(source_fields, function(field) any(grepl(field, script_text, fixed = TRUE)), logical(1))
    ),
    contains_original_project_absolute_path = any(
      grepl(
        "source\\([^\\n]*([A-Za-z]:[/\\\\]|/Users/|/home/)",
        script_text
      )
    ),
    stringsAsFactors = FALSE
  )
})

# 合并脚本QA。
script_qa <- do.call(rbind, script_qa)

# 读取预先登记的数据哈希。
manifest_path <- file.path(
  package_root,
  "10_方法来源与文件校验",
  "file_manifest_sha256.csv"
)

# 要求哈希清单存在。
if (!file.exists(manifest_path)) {
  stop("缺少file_manifest_sha256.csv。")
}

# 读取哈希清单。
manifest <- utils::read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 仅核对不属于输出目录的固定输入和代码。
fixed_manifest <- manifest[
  !grepl("^09_输出结果/", manifest$relative_path),
  ,
  drop = FALSE
]

# 逐项复算哈希。
hash_qa <- lapply(seq_len(nrow(fixed_manifest)), function(row_index) {
  # 构造实际文件路径。
  current_path <- file.path(
    package_root,
    fixed_manifest$relative_path[row_index]
  )
  # 判断文件是否存在。
  exists_now <- file.exists(current_path)
  # 存在时复算SHA-256。
  current_hash <- if (exists_now) {
    toupper(digest::digest(
      object = current_path,
      algo = "sha256",
      file = TRUE,
      serialize = FALSE
    ))
  } else {
    NA_character_
  }
  # 比较当前SHA-256与登记值。
  hash_matches <- exists_now &&
    identical(
      current_hash,
      toupper(fixed_manifest$sha256[row_index])
    )
  data.frame(
    relative_path = fixed_manifest$relative_path[row_index],
    exists = exists_now,
    registered_sha256 = fixed_manifest$sha256[row_index],
    current_sha256 = current_hash,
    sha256_matches = hash_matches,
    stringsAsFactors = FALSE
  )
})

# 合并文件存在性QA。
hash_qa <- do.call(rbind, hash_qa)

# 定义QA输出目录。
qa_dir <- file.path(output_root, "01_输入QA")

# 创建QA输出目录。
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

# 保存脚本QA。
atomic_write_csv(
  script_qa,
  file.path(qa_dir, "package_script_static_QA.csv")
)

# 保存清单存在性QA。
atomic_write_csv(
  hash_qa,
  file.path(qa_dir, "package_manifest_existence_QA.csv")
)

# 设置硬门槛。
all_pass <- all(script_qa$parse_ok) &&
  all(script_qa$source_block_complete) &&
  !any(script_qa$contains_original_project_absolute_path) &&
  all(hash_qa$exists) &&
  all(hash_qa$sha256_matches)

# 任一硬门槛失败时停止。
if (!all_pass) {
  stop("运行包静态审计失败；请检查09_输出结果/01_输入QA。")
}

# 报告审计通过。
message("运行包静态审计通过。")
