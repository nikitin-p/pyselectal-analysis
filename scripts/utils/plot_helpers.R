# Shared plotting helpers.
#
# Bars in this project almost always carry a percentage, so the underlying
# number of alignments (or of unique positions) is invisible: a category holding
# 0.01% of the data looks exactly as trustworthy as one holding 40%. Every
# histogram therefore prints that number next to the bar, formatted compactly by
# fmt_n() so the annotation stays readable at small figure sizes.
#
# Source it from a script with:
#   .f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
#   source(if (length(.f))
#            file.path(dirname(sub("^--file=", "", .f[1])), "..", "utils", "plot_helpers.R")
#          else "scripts/utils/plot_helpers.R")

# Compact count label: 1234567 -> "1.2M", 45678 -> "46k", 1234 -> "1.2k", 987 -> "987".
fmt_n <- function(n) {
  n <- as.numeric(n)
  out <- ifelse(
    is.na(n), "",
    ifelse(n >= 1e6, sprintf("%.1fM", n / 1e6),
    ifelse(n >= 1e4, sprintf("%.0fk", n / 1e3),
    ifelse(n >= 1e3, sprintf("%.1fk", n / 1e3),
           formatC(n, format = "d", big.mark = "")))))
  out
}
