# Silence R CMD check NOTEs about non-standard evaluation against data.table
# columns referenced by bare name in the wrapper code.
utils::globalVariables(c("e", "lat", "lon", "unit", "time", "keep"))
