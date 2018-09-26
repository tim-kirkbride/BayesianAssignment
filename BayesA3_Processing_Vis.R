library(readr)
library(ggplot2)
library(data.table)

#set working directory to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read data
bf<-as.data.table(read.csv("BlackFriday.csv"))

#=========== Preprocessing ===============

#drop unnecessary/ambiguous columns
bf[ ,c("Product_ID","Occupation","City_Category",
       "Product_Category_1","Product_Category_2","Product_Category_3") := NULL]

#aggregate data to ensure unique user IDs
bf[, `:=` (Count = .N, Purch_Avg = mean(Purchase), Purch_Sum = sum(Purchase)), by = User_ID]

#remove original purchase column as well as duplicates
bf[ ,c("Purchase") := NULL]
bf<-subset(unique(bf), select= -User_ID)

#saving dataset for visualisation + making marital stuts levels more interpretable
bf_vis<-as.data.frame(bf)
bf_vis$Marital_Status<-ifelse(bf_vis$Marital_Status == "1","Married","Not Married")

#converting factors to binary
setDT(bf)[, c(levels(bf$Gender), "Gender") := 
                c(lapply(levels(Gender), function(x) as.integer(x == Gender)), .(NULL))]
setDT(bf)[, c(levels(bf$Age), "Age") := 
            c(lapply(levels(Age), function(x) as.integer(x == Age)), .(NULL))]
setDT(bf)[, c(levels(bf$Stay_In_Current_City_Years), "Stay_In_Current_City_Years") := 
            c(lapply(levels(Stay_In_Current_City_Years), function(x) as.integer(x == Stay_In_Current_City_Years)), .(NULL))]


#removing base level column from each factor
bf[ ,c("F","0-17","0") := NULL]

#============== Visualisation ================

#univariate

#avg purchase per cust
ggplot(bf_vis, aes(Purch_Avg))+geom_histogram(fill="#74c476")+
  labs(title="Histogram of Average Purchase Price per Customer",
       y="Frequency",
       x="Average Purchase Price ($)")+theme_minimal()

#total purchase price per cust
ggplot(bf_vis, aes(Purch_Sum))+geom_histogram(fill="#74c476")+
  labs(title="Histogram of Total Purchase Price per Customer",
       y="Frequency",
       x="Total Purchase Price ($)")+theme_minimal()

#avg purch vs gender
ggplot(bf_vis, aes(x=Gender, y=Purch_Avg))+
  geom_boxplot(fill="#74c476")+
  labs(title="Average Purchase Price per Customer by Gender",
       y="Average Purchase Price ($)",
       x="Gender")+theme_minimal()

#avg purch vs age
ggplot(bf_vis, aes(x=Age, y=Purch_Avg))+
  geom_boxplot(fill="#74c476")+
  labs(title="Average Purchase Price per Customer by Age",
       y="Average Purchase Price ($)",
       x="Age")+theme_minimal()

#avg purch vs city stay
ggplot(bf_vis, aes(x=Stay_In_Current_City_Years, y=Purch_Avg))+
  geom_boxplot(fill="#74c476")+
  labs(title="Average Purchase Price per Customer by \nYears Stayed in City",
       y="Average Purchase Price ($)",
       x="Years Stayed in City")+theme_minimal()

#avg purch vs marital
ggplot(bf_vis, aes(x=Marital_Status, y=Purch_Avg))+
  geom_boxplot(fill="#74c476")+
  labs(title="Average Purchase Price per Customer by Marital Status",
       y="Average Purchase Price ($)",
       x="Marital Status")+theme_minimal()

#avg purch vs count
ggplot(bf_vis, aes(x=Count, y=Purch_Avg))+
  geom_point(colour="#74c476")+
  labs(title="Average Purchase Price per Customer by \nNumber of Purchases",
       y="Average Purchase Price ($)",
       x="Number of Purchases")+theme_minimal()

#multivariate

ggplot(bf_vis, aes(x=Age, y=Purch_Avg))+
  geom_boxplot(aes(fill=Gender))

ggplot(bf_vis, aes(x=Count, y=Purch_Avg))+
  geom_point(colour="#74c476")+facet_wrap(~Age)

ggplot(bf_vis,aes(x=Purch_Avg,fill=Gender))+
  geom_density(alpha=.5)+facet_grid(Age~Stay_In_Current_City_Years)
