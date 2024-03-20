if (requireNamespace("tinytest", quietly=TRUE) ){
  tinytest::test_package("plmmr", pattern="^[^_].*\\.[rR]$")
}
