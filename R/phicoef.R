### Based on http://en.wikipedia.org/wiki/Phi_coefficient

phicoef <- function(x, y=NULL)
{
    if (is.null(y)) {
        if (!is.integer(x) || length(x) != 4L)
            stop("when 'y' is not supplied, 'x' must be ",
                 "a 2x2 integer matrix or an integer vector of length 4")
        a <- x[1L]
        c <- x[2L]
        b <- x[3L]
        d <- x[4L]
    } else {
        if (!is.logical(x) || !is.logical(y) || length(x) != length(y))
            stop("when 'y' is supplied, 'x' and 'y' must be ",
                 "2 logical vectors of the same length")
        a <- sum(x  &  y)
        b <- sum(x  & !y)
        c <- sum(!x &  y)
        d <- sum(!x & !y)
    }
    a <- as.double(a)
    b <- as.double(b)
    c <- as.double(c)
    d <- as.double(d)
    div <- sqrt((a + b) * (c + d) * (a + c) * (b + d))
    (a * d - b * c) / div
}

