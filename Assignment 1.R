##Read file in .

## Training dataset
King = read.csv(file.choose())


##Test dataset
KingTest = read.csv(file.choose())


library(car)



basemodel <- lm(price ~. , data = King)

pairs(basemodel)


residualPlot(basemodel)

Wpairs(King)

basemodel

summary(King)
