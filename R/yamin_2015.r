data <- readRDS("yamin_2015_data.rds")

days <- seq(as.Date("2014-5-1"), as.Date("2014-12-1"), by = "day")

viralLoad <- 1

generateR <- function(viralLoadIn_)
{
  RcoefOUT <- 0
  while (RcoefOUT <= 1)
  {
    if (viralLoadIN_ == 1)
    {
      meal <- max(1, exp(rnorm(1, log(1.2), 0.415)))
      RRmeal <- max(meal, exp(rnorm(1, log(2.2), 0.31))) / meal
      RcoefOUT <-
        RRmeal**(1 / ((52 * sum(data[["fatalMean"]][5:10]) / 6 +
                       3 * sum(data[["nonFatalMean"]][5:10]) / 6) / 55 -
                      (52 * sum(data[["fatalMean"]][1:4]) / 4 +
                       3 * sum(data[["nonFatalMean"]][1:4]) / 4) / 55))
    } else if (viralLoadIN_ == 2)
    {
      talk <- max(1, exp(rnorm(1, log(0.7), 0.525)))
      RRtalk <- max(talk, exp(rnorm(1, log(3.9), 0.58))) / talk
      RcoefOUT <-
        RRtalk**(1 / ((52 * sum(data[["fatalMean"]][5:10]) / 6 +
                       3 * sum(data[["nonFatalMean"]][5:10]) / 6) / 55 -
                      (52 * sum(data[["fatalMean"]][1:4]) / 4 +
                       3 * sum(data[["nonFatalMean"]][1:4]) / 4) / 55))
    } else if (viralLoadIN_ == 3)
    {
      bed <- max(1, exp(rnorm(1, log(1.3), 0.34)))
      RRbed <- max(bed, exp(rnorm(1, log(2.2), 0.33)))
      RcoefOUT <-
        RRbed**(1 / ((52 * sum(data[["fatalMean"]][5:10]) / 6 +
                      3 * sum(data[["nonFatalMean"]][5:10]) / 6) / 55 -
                     (52 * sum(data[["fatalMean"]][1:4]) / 4 +
                      3 * sum(data[["nonFatalMean"]][1:4]) / 4) / 55))
    }
    return(RcoefOUT)
  }
}

define <- function(iterationsIn_)
{
  numOfIterations <- iterationsIn_
  maxSickPeriod <- 14
  minSickPeriod <- 5
  maxIncubationPeriod <- 15
  minIncubationPeriod <- 5

  incubationDist <-
    sapply(seq_len(numOfIterations),
           function (x) {
             round(rtriangle(length(data[["id"]]),
                             minIncubationPeriod - 0.5,
                             maxIncubationPeriod + 0.5,
                             8))
           })

  sickPeriodDist <-
    sapply(seq_len(numOfIterations),
           function (x) {
             round(rtriangle(length(data[["id"]]),
                             minSickPeriod - 0.5,
                             maxSickPeriod + 0.5,
                             8))
           })

  contactFirstStageDist <-
    sapply(seq_len(numOfIterations),
           function (x) {
             sample(data[["numContacts"]], length(id), replace = TRUE)
           })

  contactSecondStageDist <-
    sapply(seq_len(numOfIterations),
           function (x) {
             sample(data[["numContacts"]][data[["numContacts"]] <= 5],
                    length(id), replace = TRUE)
           })

  stageCutoffDist <- 
    sapply(seq_len(numOfIterations),
           function (x) {
             round(runif(length(data[["id"]]), min = 1, max = 5))
           })

  Rcoef <- 
    sapply(seq_len(numOfIterations), function(x) generateR(viralLod))

  Rcoef <- Table[generateR[viralLoad],{i,1,numOfIterations}]
  Rcoef=Table[generateR[viralLoad],{i,1,numOfIterations}]
}
