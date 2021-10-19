#' ---
#' title: "Class 05 Data Visualization"
#' author: "Dominique Lie - A15470100"
#' date: "October 12, 2021"
#' ---


# Class 05 Data Visualization

# Lets start with a scatterplot
# Before we can use ggplot2, we need to load it up!
# install.pacakges(ggplot2) before library(ggplot2)

library(ggplot2)

# Every ggplot has a data + aes + geoms
ggplot(data = cars) + aes (x = speed, y = dist) + geom_point() + geom_smooth()

# Change to linear model
# Save plot as p
p <- ggplot(data = cars) + aes (x = speed, y = dist) + geom_point() + geom_smooth(method = "lm")

p + labs(title = "My nice plot", x = "Speed (mph)", y = "distance (feet)")

#Base graphics is shorter
plot(cars)

#Adding more plot aesthetics (size, color, alpha)
url <-"https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)

# Function to find out how many genes are in dataset
nrow(genes)

# colnames() for column names and ncol() for number of columns
colnames(genes)
ncol(genes)

# table() for "state"upregulated" genes - just to access list genes["State"] or genes$State
table(genes$State)

# Q. What % are up/down? round function for significant figures (,2) means two sigfig
prec <- table(genes$State) / nrow(genes) * 100
round(prec, 2)

# Make scatterplot with new data genes
ggplot( data = genes) + aes( x = Condition1, y = Condition2, ) + geom_point()

# Add state column and color
ggplot( data = genes) + aes( x = Condition1, y = Condition2, col = State ) + geom_point()
p <- ggplot( data = genes) + aes( x = Condition1, y = Condition2, col = State ) + geom_point()

# Change default colors by adding another layer
p + scale_colour_manual( values = c("blue", "gray", "red"))
p <- p + scale_colour_manual( values = c("blue", "gray", "red"))

# Add plot notations to change name
p + labs( title = "Gene Expression Changes Upon Drug Treatment", x = "Control (no drug)", y = "Drug Treatment")

# Optional gapminder set install.packages("gapminder")
library(gapminder)

# To focus on a single year install.packages("dplyr")
# Filter data frame
library(dplyr)
gapminder_2007 <-  gapminder %>% filter(year==2007)

# basic scatterplot
ggplot( data = gapminder_2007) + aes( x= gdpPercap, y = lifeExp) + geom_point()

# make plot points transparent
ggplot( data = gapminder_2007) + aes( x= gdpPercap, y = lifeExp) + geom_point(alpha = 0.4)

ggplot( data = gapminder_2007) + aes( x= gdpPercap, y = lifeExp, color = continent, size = pop) + geom_point(alpha = 0.4)
ggplot(gapminder_2007) + aes( x = gdpPercap, y = lifeExp, color = pop) + geom_point(alpha = 0.8)
ggplot(gapminder_2007) + aes( x = gdpPercap, y = lifeExp, color = pop, size = pop) + geom_point(alpha = 0.8)
ggplot(gapminder_2007) + aes( x = gdpPercap, y = lifeExp, color = pop, size = pop) + geom_point(alpha = 0.8) + scale_size_area()



