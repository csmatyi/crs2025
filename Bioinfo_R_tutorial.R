# Creation Research Society conference July 22, 2025, St. Louis

### BASICS

# get the working directory that you are in
getwd()
# set this to your own directory
setwd("C:\\Users\\csmat\\OneDrive\\Desktop\\CreationScience\\Talks\\CRSconferences\\CRS2025\\workshop")

# math operators
11 + 9
12 - 4
2 * 6
7/3
15 %% 4

# NA types values
NA + 2
w <- NA
is.na(w)
1/0
-1/0
is.nan(Inf)
is.nan(0/0)
is.nan(1/0)
as.integer(9.9)
as.double(11)

# TRUE and FALSE
crs_is_cool <- TRUE
crs_is_cool
as.integer(T)
as.integer(F)
as.integer(TRUE)
as.integer(FALSE)
typeof(7)
typeof(2.4)
typeof(FALSE)
typeof("ACGTGGC")

# use built-in math functions
sqrt(9)
sqrt(25)
sin(pi/4) # here pi = 3.14159, pi/4is 45 degrees
log10(100)
log2(16)
log(121,base=11)
abs(-4)
factorial(5)
# n!/(k! x (n-k)!)
choose(5,3)

# print text
print("Hello world!")
print(7/9)

# change no. digits
options(digits=20)
print(7/9)
print(7/9,digits=2)
options(digits=3)

# concatenate
paste("The","rain","in","Spain","stays","mainly","in","the","plain")
paste("The","rain","in","Spain","stays","mainly","in","the","plain",sep="*")
paste("A","C","G","T",sep="")
# with paste0 the default separator is ""
paste0("A","C","G","T")

### Exercise 1 on DNA sequence manipulation ###



# look for help for a given function
help(sqrt)
?sqrt # same thing
help.search("log")
??"log"
# get an example of how to run a function
example(log)
# assign values to variables
X <- 1
1 -> X
X = 1

# create a list
v <- c(4, "ACCGT", FALSE)
x <- c(56, 95.3, 0.4)
z <- seq(2,5)
z2 <- seq(2,5,0.01)

# use sequence method
a <- 3
a + x
a - x
a * x
a / x
sqrt(x)

# reference list element
x[3]
vec <- c(a=3.4, b=5.4, c=0.4) # named list

# look at list name
names(vec)

# change list name
names(vec)[3] <- 'd'
vec

# get subset of list
x[c(2,3)]
nums <- seq(1:10)
nums[c(6,2,9)]
nums[3:6]
nums[5:1]

# order list in decreasing order
order(nums,decreasing=TRUE)

# random sample
sample(10)

# raise list elements by exponent
nums2 <- nums ** 2
nums2

# summarize
summary(nums2)

# sample list
nums2[sample(10)]

# sample list with replacement
nums2[sample(10, replace=TRUE)]

# check which elements in x are > 10
x > 10

# get those elements
x[x > 10]

# get elements via T/F
nums[c(T, T, F, F, T, F, F, F, F, T)]

# create matrix
mx <- matrix(c(11.5, -0.9, 5.6, 7.7, 13.5, 1.8), nrow = 3, ncol = 2)
mx

# matrix dimensions
dim(mx)
dim(mx)[1]
mx[1,2]

# out of bounds index
mx[1:2,3]

# matrix of zeroes
# 0, num rows, num cols
matrix(0,3,4)

# multiply 2 matrices together
mx1 <- matrix(c(11.5, -0.9, 5.6, 7.7, 13.5, 1.8, 0, -9.1), nrow = 4, ncol = 2)
mx1
mx2 <- matrix(c(1,-2,4,1,0,-9), nrow=2, ncol=3)
mx2
# n x k  X  k x m = n X m matrix
mx1 %*% mx2

# create a data frame of students and their grades
course <- data.frame(
  names = c("John", "Mary", "Scott"),
  scores = c(98,95,89),
  letter = c("A","A","B")
)
course

# get the column names and row names
colnames(course)
rownames(course)

# summarize the data frame
summary(course)

# new row
course <- rbind(course, c("Jan",72,"C"))
course

# new column
course <- cbind(course, major=c("compsci","compsci","art","history"))
course

# remove last row
course <- course[-c(4),]
course

# remove last col
course <- course[,-c(4)]
course

# functions, loops
# declare method that squares a number
square <- function(n) {
  m = n*n
  return(m)
}

# run it
square(5)

# show method code
square

### Exercise #2: Hardy-Weinberg function ###



# loop
for (i in 1:5) { 
  print(i) 
}

# if statement
x <- 10
if (x > 5) {
  print("x is greater than 5")
}

# if else statement
if (x > 5) {
  print("x is greater than 5")
} else {
  print("x is less than 5")
}

### Exercise 3 on dummy RNA-seq data


### Review slides on applying functions to data ###

# apply, lapply, sapply, mapply
# declare matrix
mxx <- matrix(1:12,nrow=3,ncol=4)
mxx

# apply mean function to rows
apply(mxx,1,mean)
# first row: (1+4+7+10)/4 = 5.5

# apply mean function to cols
apply(mxx,2,mean)
# (1+2+3)/3 = 2

# list apply: apply function to all elements in nums list
# new list
nums <- c(2,4,7)
lnums <- lapply(nums, sqrt)
lnums
lnums[[2]]

# new function
cubit <- function(n) {return(n*n*n)}

# lapply produces a list
lapply(nums, cubit)

# sapply (simple apply) produces a one-line list, a vector, similar to lapply
# except that formatting is different
sapply(nums, cubit)

# using SIMPLIFY in mapply to write out results on one line
mapply(cubit,nums,SIMPLIFY = FALSE)
mapply(cubit,nums,SIMPLIFY = TRUE)
mapply(sqrt,nums,SIMPLIFY = TRUE)

# use mapply with function created on the fly :-)
# the function is 3x + 1
mapply(function(x) 3*x+1,c(1,2,3,4,5))
