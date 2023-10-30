if (requireNamespace("tinytest", quietly=TRUE) ){
  tinytest::test_package("plmm", pattern="^[^_].*\\.[rR]$")
}
