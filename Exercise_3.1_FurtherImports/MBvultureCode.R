library(move)
library(RCurl)
library(circular)

login <- movebankLogin(username="wdwalter", password="xxxxxx")

getMovebankAnimals(study="Turkey Vulture South Carolina USA",login=login)

turkey <- getMovebankData(study="Turkey Vulture South Carolina USA", animalName=c("Bird51","Bird52","Bird54","Bird55a","Bird56","Bird59a","Bird60a"),login=login, moveObject=TRUE)

n.locs(turkey)
# Bird51  Bird52  Bird54 Bird55a  Bird56 Bird59a Bird60a 
#   9655   11456    2378    2228    5876    8311    8594 
   