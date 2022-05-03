.onAttach <- function(libname, pkgname) {
    # packageStartupMessage("Welcome Zhuoer Dong's little package")
}

.onLoad <- function(libname, pkgname) {
    # rGEO.data.options <- list(
    #   rGEO.data.foo = "bar"
    # )

    # to_set <- !(names(rGEO.data.options) %in% names(options()))
    # if(any(to_set)) options(rGEO.data.options[to_set])

    invisible()
}
