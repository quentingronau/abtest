require(devtools)
options(error = NULL)

load_all()
roxygen2::roxygenize()

test()
