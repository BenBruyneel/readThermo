% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readFunctionsThermo.R
\name{readProcessing.Thermo}
\alias{readProcessing.Thermo}
\title{readProcessing.Thermo}
\usage{
readProcessing.Thermo(
  filename,
  leaveIn = NA,
  leaveOut = c("Component", "mdlCalcs")
)
}
\arguments{
\item{filename}{name of the file from which the data is to be read. Must be
an export file in excel format (.XLS), see Xcalibur software manual.}

\item{leaveIn}{character vector (or NA). Defines which sheets (components) to
extract. Default is NA, which causes all sheets to be read (except those
defined im the 'leaveout' argument).}

\item{leaveOut}{character vector which defines which sheets to leave out of
the result. The internal functions (packages readxl & XLConnect) see two
extra sheets which do not contain info on the components defined in the
processing method. They have to do with formatting & macros. This package
doesn't do anything with that, so this argument can be used to remove them
from the result. If set to NA, then every sheet in the file is read.}
}
\value{
NA (in case it fails reading the sheets from the file) or a list of
data.frame's containing the integration data (etc) in the .XLS file exported
from Xcalibur
}
\description{
reads in the excel files generated by the processing methods in
the Thermo Scientific XCalibur software. The export file is an excel file
in which the integration data is stored per sheet
}
\examples{
demoXLS <- fs::path_package("extdata", "ExcelExp_Long.XLS", package = "readThermo")
result <- readProcessing.Thermo(filename = demoXLS)
result[[1]]

}
