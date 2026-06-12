# Build the pkgdown site, then scrub CLAUDE.md (internal dev instructions)
# from the published output. pkgdown renders every top-level .md file with
# no exclude option, so the removal has to happen post-build.
#
# Usage: Rscript pkgdown/build-site.R

pkgdown::build_site(preview = FALSE)

docs <- "docs"

unlink(file.path(docs, c("CLAUDE.html", "CLAUDE.md")))

sitemap <- file.path(docs, "sitemap.xml")
lines <- readLines(sitemap)
writeLines(lines[!grepl("CLAUDE", lines)], sitemap)

search_json <- file.path(docs, "search.json")
entries <- jsonlite::fromJSON(readLines(search_json), simplifyVector = FALSE)
keep <- vapply(entries, function(e) {
  # some entries (e.g. the lunr settings stub) carry no path field
  is.null(e$path) || !grepl("CLAUDE", e$path, fixed = TRUE)
}, logical(1))
writeLines(jsonlite::toJSON(entries[keep], auto_unbox = TRUE), search_json)

stopifnot(!any(grepl("CLAUDE", list.files(docs))),
          !any(grepl("CLAUDE", readLines(search_json))))
cat("Site built and CLAUDE pages scrubbed.\n")
