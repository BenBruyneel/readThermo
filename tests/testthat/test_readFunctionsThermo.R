test_that("readChromatogram.Thermo works", {
  demoRaw <- fs::path_package(
    "extdata",
    "reserpine07.RAW",
    package = "readThermo"
  )
  result <- readChromatogram.Thermo(
    filename = demoRaw,
    type = "tic",
    filter = "ms2"
  )
  expect_equal(class(result), "function")
  result <- result()
  expect_equal(class(result), "list")
  expect_length(result, 1)
  expect_equal(class(result[[1]]), c("dataElement","R6"))
  expect_length(result[[1]]$info, 6)
  expect_equal(
    names(result[[1]]$info),
    c("source", "filename", "mz", "tolerance", "filter", "type")
  )
  expect_equal(result[[1]]$info$source, "thermo")
  demoRaw <- fs::path_package(
    "extdata",
    "reserpine07.RAW",
    package = "readThermo"
  )
  result <- readChromatogram.Thermo(
    filename = demoRaw,
    type = "xic",
    filter = "ms2",
    mz = c(397.16, 448.13),
    tolerance = 500
  )
  expect_equal(class(result), "function")
  result <- result()
  expect_equal(class(result), "list")
  expect_length(result, 2)
  expect_equal(class(result[[1]]), c("dataElement","R6"))
  expect_length(result[[1]]$info, 6)
  expect_equal(
    names(result[[1]]$info),
    c("source", "filename", "mz", "tolerance", "filter", "type")
  )
  expect_equal(result[[1]]$info$source, "thermo")
  expect_equal(class(result), "list")
  expect_length(result, 2)
  expect_length(result[[1]][["info"]], 6)
  expect_equal(
    names(result[[1]]$info),
    c("source", "filename", "mz", "tolerance", "filter", "type")
  )
  expect_equal(as.character(result[[1]]$info$mz), "397.16")
  expect_equal(names(result[[1]]$data), c("rt", "intensity"))
  expect_equal(nrow(result[[1]]$data), 403)
  expect_equal(as.character(min(result[[1]]$data[, 1])), "0.009")
  expect_equal(as.character(max(result[[1]]$data[, 1])), "7.00066666666667")
  expect_equal(as.character(min(result[[1]]$data[, 2])), "0")
  expect_equal(as.character(max(result[[1]]$data[, 2])), "5968176")
  expect_length(result[[2]][["info"]], 6)
  expect_equal(
    names(result[[2]]$info),
    c("source", "filename", "mz", "tolerance", "filter", "type")
  )
  expect_equal(as.character(result[[2]]$info$mz), "448.13")
  expect_equal(names(result[[2]]$data), c("rt", "intensity"))
  expect_equal(nrow(result[[2]]$data), 403)
  expect_equal(as.character(min(result[[2]]$data[, 1])), "0.009")
  expect_equal(as.character(max(result[[2]]$data[, 1])), "7.00066666666667")
  expect_equal(as.character(min(result[[2]]$data[, 2])), "0")
  expect_equal(as.character(max(result[[2]]$data[, 2])), "4614730")
})

test_that("readSpectrum.Thermo works", {
  demoRaw <- fs::path_package(
    "extdata",
    "reserpine07.RAW",
    package = "readThermo"
  )
  result <- readSpectrum.Thermo(
    filename = demoRaw,
    scan = 245,
    centroided = TRUE
  )
  expect_equal(class(result), "function")
  result <- result()
  expect_equal(class(result), "list")
  expect_length(result, 1)
  expect_equal(class(result[[1]]), c("dataElement","R6"))
  expect_length(result[[1]]$info, 6)
  expect_equal(
    names(result[[1]]$info),
    c("source", "filename", "rt", "scan", "scanType", "centroided")
  )
  expect_equal(result[[1]]$info$scan, 245)
  expect_equal(names(result[[1]]$data), c("mz", "intensity"))
  expect_equal(nrow(result[[1]]$data), 47)
  expect_equal(as.character(min(result[[1]]$data[, 1])), "174.356506347656")
  expect_equal(as.character(max(result[[1]]$data[, 1])), "610.101440429688")
  expect_equal(as.character(min(result[[1]]$data[, 2])), "2")
  expect_equal(as.character(max(result[[1]]$data[, 2])), "5968176")
  result <- readSpectrum.Thermo(
    filename = demoRaw,
    scan = 245,
    centroided = TRUE
  )
  expect_equal(class(result), "function")
  result <- result(parameters = "extended")
  expect_equal(class(result), "list")
  expect_length(result, 1)
  expect_equal(class(result[[1]]), c("dataElement","R6"))
  expect_length(result[[1]]$info, 9)
  expect_equal(
    names(result[[1]]$info),
    c(
      'source',
      'filename',
      'rt',
      'scan',
      'scanType',
      'centroided',
      'microScanCount',
      'ionInjectionTime',
      'elapsedTime'
    )
  )
  expect_equal(as.character(result[[1]]$info$ionInjectionTime), "24.33")
  result <- readSpectrum.Thermo(
    filename = demoRaw,
    scan = 245,
    centroided = TRUE
  )
  expect_equal(class(result), "function")
  result <- result(parameters = "full")
  expect_equal(class(result), "list")
  expect_length(result, 1)
  expect_equal(class(result[[1]]), c("dataElement","R6"))
  expect_length(result[[1]]$info, 25)
  result <- readSpectrum.Thermo(
    filename = demoRaw,
    scan = c(245, 246),
    centroided = TRUE
  )
  expect_equal(class(result), "function")
  result <- result()
  expect_equal(class(result), "list")
  expect_length(result, 2)
  expect_equal(class(result[[1]]), c("dataElement","R6"))
  expect_length(result[[1]]$info, 6)
  expect_equal(result[[1]]$info$scan, 245)
  expect_equal(class(result[[2]]), c("dataElement","R6"))
  expect_length(result[[2]]$info, 6)
  expect_equal(result[[2]]$info$scan, 246)
  expect_equal(nrow(result[[1]]$data), 47)
  expect_equal(as.character(min(result[[1]]$data[, 1])), "174.356506347656")
  expect_equal(as.character(max(result[[1]]$data[, 1])), "610.101440429688")
  expect_equal(as.character(min(result[[1]]$data[, 2])), "2")
  expect_equal(as.character(max(result[[1]]$data[, 2])), "5968176")
  expect_equal(nrow(result[[1]]$data), 47)
  expect_equal(as.character(min(result[[2]]$data[, 1])), "173.728240966797")
  expect_equal(as.character(max(result[[2]]$data[, 1])), "609.599975585938")
  expect_equal(as.character(min(result[[2]]$data[, 2])), "1")
  expect_equal(as.character(max(result[[2]]$data[, 2])), "4878015")
})

test_that("fileInfo.Thermo works", {
  demoRaw <- fs::path_package("extdata", "reserpine07.RAW", package = "readThermo")
  result <- fileInfo.Thermo(demoRaw, readIndex = TRUE)()
  expect_length(result[[1]]$info, 38)
  expect_equal(result[[1]]$info$filename,demoRaw)
  expect_equal(result[[1]]$info$RAWfile, "reserpine07.RAW")
  expect_equal(result[[1]]$info$Numberofscans, 403)
  expect_equal(ncol(result[[1]]$data), 9)
  expect_equal(colnames(result[[1]]$data),
                        c('scan','scanType','StartTime', 'precursorMass',
                          'MSOrder','charge','masterScan','dependencyType',
                          'monoisotopicMz'))
  expect_equal(nrow(result[[1]]$data), 403)
  expect_equal(as.character(unique(result[[1]]$data$precursorMass)), "609.2")
  result <- fileInfo.Thermo(demoRaw, readIndex = FALSE)()
  expect_length(result[[1]]$info, 38)
  expect_equal(result[[1]]$info$filename,demoRaw)
  expect_equal(result[[1]]$info$RAWfile, "reserpine07.RAW")
  expect_equal(result[[1]]$info$Numberofscans, 403)
  expect_equal(ncol(result[[1]]$data), 1)
  expect_equal(colnames(result[[1]]$data),
               "data")
  expect_equal(nrow(result[[1]]$data), 1)
  expect_equal(result[[1]]$data$data, "No Data")
})

test_that("readProcessing.Thermo works", {
  demoXLS <- fs::path_package(
    "extdata",
    "ExcelExp_Long.XLS",
    package = "readThermo"
  )
  result <- readProcessing.Thermo(filename = demoXLS)
  expect_length(result, 3)
  expect_equal(dim(result[[1]]), c(40, 47))
  expect_equal(dim(result[[2]]), c(40, 47))
  expect_equal(dim(result[[3]]), c(40, 47))
  result <- readProcessing.Thermo(filename = demoXLS, leaveIn = "One")
  expect_length(result, 1)
  expect_equal(dim(result[[1]]), c(40, 47))
  demoXLS <- fs::path_package(
    "extdata",
    "ExcelExp_Short.XLS",
    package = "readThermo"
  )
  result <- readProcessing.Thermo(filename = demoXLS)
  expect_length(result, 3)
  expect_equal(dim(result[[1]]), c(40, 16))
  expect_equal(dim(result[[2]]), c(40, 16))
  expect_equal(dim(result[[3]]), c(40, 16))
})

test_that("processingThermo.componentName", {
  demoXLS <- fs::path_package(
    "extdata",
    "ExcelExp_Long.XLS",
    package = "readThermo"
  )
  result <- readProcessing.Thermo(filename = demoXLS)
  expect_equal(processingThermo.componentName(result), c("One", "Two", "Three"))
  expect_equal(processingThermo.componentName(result[[1]]), "One")
})

test_that("processingThermo.creatorInfo", {
  demoXLS <- fs::path_package(
    "extdata",
    "ExcelExp_Long.XLS",
    package = "readThermo"
  )
  result <- readProcessing.Thermo(filename = demoXLS)
  theInfo <- processingThermo.creatorInfo(result)
  expect_equal(dim(theInfo), c(3, 3))
  expect_equal(colnames(theInfo), c("User Name", "Full Name", "Date"))
  expect_equal(unique(theInfo[, 1]), "Ben")
  expect_equal(unique(theInfo[, 2]), "Ben")
  theInfo <- processingThermo.creatorInfo(result[[1]])
  expect_equal(dim(theInfo), c(1, 3))
  expect_equal(colnames(theInfo), c("User Name", "Full Name", "Date"))
  expect_equal(unique(theInfo[, 1]), "Ben")
  expect_equal(unique(theInfo[, 2]), "Ben")
})

test_that("processingThermo.componentData", {
  demoXLS <- fs::path_package(
    "extdata",
    "ExcelExp_Long.XLS",
    package = "readThermo"
  )
  result <- readProcessing.Thermo(filename = demoXLS)
  theInfo <- processingThermo.componentData(result)
  expect_equal(dim(theInfo), c(42, 48))
  expect_equal(unique(theInfo$Component), c("One", "Two", "Three"))
  theInfo <- processingThermo.componentData(result[[1]])
  expect_equal(dim(theInfo), c(14, 48))
  expect_equal(unique(theInfo$Component), "One")
  demoXLS <- fs::path_package(
    "extdata",
    "ExcelExp_Short.XLS",
    package = "readThermo"
  )
  result <- readProcessing.Thermo(filename = demoXLS)
  theInfo <- processingThermo.componentData(result)
  expect_equal(dim(theInfo), c(42, 17))
  expect_equal(unique(theInfo$Component), c("One", "Two", "Three"))
  theInfo <- processingThermo.componentData(result[[1]])
  expect_equal(dim(theInfo), c(14, 17))
  expect_equal(unique(theInfo$Component), "One")
})

# No test for
# - getScans
#
# If the example(s) work it should be good
