# Text plots?

# Load packages
require(txtplot)

# Load test data
datobj <- read.csv("deseq-test.csv", header = TRUE, row.names = 1)



# MA
txtplot(
	x = log10(datobj[, 1]), 
	y = datobj[, 2], 
	xlab = "log 2 fold change",
	ylab = "normalized mean",
	ylim = c(-4, 4)
)


# Volcano
txtplot(
	x = datobj[, 2],
	y = -log10(datobj[, 4]),
	xlab = "log2 fold change",
	ylab = "-log10(p-value)",
	xlim = c(-10, 10),
	ylim = c(0, 250)
)


# Output text
writeLines(
	capture.output(
		txtplot(
			x = datobj[, 2],
			y = -log10(datobj[, 4]),
			xlab = "log2 fold change",
			ylab = "-log10(p-value)",
			xlim = c(-10, 10),
			ylim = c(0, 250)
		)
	),
	con = "output.txt",
	sep = "\n"
)

tmp <- read.csv("ia.csv", header = TRUE)

plot(tmp$lng, tmp$lat)


library(d3scatter)
library(crosstalk)


