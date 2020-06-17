# test on slc24a5 germline unviability

# values from excel file

inx <- c(0.034482759, 0.019230769, 0.23566879)
outx <- c(0.022222222, 0.028985507, 0)

t.test(inx, outx)