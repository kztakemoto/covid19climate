library(incidence)
data <- read.csv("travel_restrictions_covid19_wiki", encoding = "UTF-8", sep="\t", header = T)
d_sub <- unique(data$Country)
sub <- data[d_sub,]
sub$Date <- as.Date(sub$Date, tryFormats = c("%m/%d/%y"))
sub_level <- na.omit(sub[c("Country","Date","Level")],)
i <- incidence(sub_level$Date, groups = sub_level$Level)
plot(i, border="white")
