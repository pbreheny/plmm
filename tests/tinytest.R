if (requireNamespace("tinytest", quietly=TRUE) ){
  tinytest::test_package("penalizedLMM", pattern="^[^_].*\\.[rR]$")
}
